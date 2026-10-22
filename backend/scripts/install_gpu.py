"""Resume, verify and install the pinned official Windows CUDA wheel.

Does not change drivers, global Python, or system security policy. Stops before
package changes if the local API is still running, avoiding locked torch DLLs.
"""
import hashlib
import os
from pathlib import Path
import socket
import subprocess
import sys
import msvcrt

root = Path(__file__).resolve().parents[2]
wheel = root / '.tools' / 'torch-2.11.0+cu128-cp312-cp312-win_amd64.whl'
expected_size = 2753189216
expected_sha = '7c78215c3af4f62e63f2b2e360f1722fc719b0853c7ac22666483d9810613a4c'
url = 'https://download-r2.pytorch.org/whl/cu128/torch-2.11.0%2Bcu128-cp312-cp312-win_amd64.whl'
wheel.parent.mkdir(parents=True, exist_ok=True)
# Lock held until process exit, including during installation. This prevents a
# second setup invocation from writing over an active resumable download.
download_lock = (wheel.parent / 'gpu-install.lock').open('a+b')
download_lock.seek(0, os.SEEK_END)
if download_lock.tell() == 0:
    download_lock.write(b'1')
    download_lock.flush()
download_lock.seek(0)
try:
    msvcrt.locking(download_lock.fileno(), msvcrt.LK_NBLCK, 1)
except OSError:
    raise SystemExit('A GPU download/install is already running. Wait for it to finish before retrying.')
if not wheel.exists() or wheel.stat().st_size < expected_size:
    print('Resuming official CUDA wheel download. Ctrl+C safely pauses it.', flush=True)
    subprocess.run(['curl.exe', '--fail', '--location', '--continue-at', '-', '--connect-timeout', '20',
                    '--silent', '--show-error', '--speed-time', '120', '--speed-limit', '1024',
                    '--retry', '5', '--retry-delay', '5', '--output', str(wheel), url], check=True)
with wheel.open('rb') as stream:
    actual_sha = hashlib.file_digest(stream, 'sha256').hexdigest()
if wheel.stat().st_size != expected_size or actual_sha != expected_sha:
    raise RuntimeError('Wheel checksum mismatch. Installation was not attempted.')
print('Official wheel SHA-256 verified.', flush=True)
if '--download-only' in sys.argv:
    print('Download complete. Close the API and run setup_gpu.cmd to install and verify all GPU engines.', flush=True)
    raise SystemExit(0)
with socket.socket() as probe:
    probe.settimeout(1)
    if probe.connect_ex(('127.0.0.1', 8000)) == 0:
        raise SystemExit('Download is ready. Stop the API on port 8000, then run setup_gpu.cmd again to install.')
subprocess.run([sys.executable, '-m', 'pip', 'install', str(wheel)], check=True, cwd=root)
subprocess.run([sys.executable, '-m', 'pip', 'install', '-e', f'{root}[test]'], check=True, cwd=root)
subprocess.run([sys.executable, '-m', 'pip', 'check'], check=True, cwd=root)
env = dict(os.environ, PYTHONPATH=str(root / 'backend'))
subprocess.run([sys.executable, str(root / 'backend' / 'scripts' / 'verify_engines.py'), 'cuda'], check=True, cwd=root, env=env)
