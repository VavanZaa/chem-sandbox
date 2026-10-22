"""Resumable bounded-range download for this host's unreliable large transfers.

Only downloads HTTPS wheel artifacts into the project .tools directory. Requires
an expected byte count and SHA-256 from the package publisher; never installs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from pathlib import Path
import shutil
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('url')
parser.add_argument('filename')
parser.add_argument('size', type=int)
parser.add_argument('sha256')
parser.add_argument('--workers', type=int, default=6)
args = parser.parse_args()
if not args.url.startswith('https://') or Path(args.filename).name != args.filename or not args.filename.endswith('.whl'):
    raise ValueError('An HTTPS URL and plain wheel filename are required')
root = Path(__file__).resolve().parents[2] / '.tools'
target = root / args.filename
parts = root / (args.filename + '.parts')
parts.mkdir(parents=True, exist_ok=True)
chunk = 1024 * 1024
count = (args.size + chunk - 1) // chunk

def checksum(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def fetch(i):
    start, stop = i * chunk, min((i + 1) * chunk, args.size) - 1
    path = parts / f'{i:05d}.part'
    if path.exists() and path.stat().st_size == stop - start + 1:
        return i
    for attempt in range(5):
        result = subprocess.run(['curl.exe', '-4', '--fail', '--silent', '--show-error', '--connect-timeout', '10',
                                 '--max-time', '25', '--range', f'{start}-{stop}', '--output', str(path), args.url],
                                capture_output=True)
        if result.returncode == 0 and path.stat().st_size == stop - start + 1:
            return i
        time.sleep(min(attempt + 1, 3))
    raise RuntimeError(f'Chunk {i} failed: {result.stderr.decode(errors="replace")}')

if target.exists() and target.stat().st_size == args.size and checksum(target) == args.sha256:
    print(f'Already verified: {target}', flush=True)
else:
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=min(max(args.workers, 1), 8)) as pool:
        for done, task in enumerate(as_completed([pool.submit(fetch, i) for i in range(count)]), 1):
            task.result()
            if done % 50 == 0 or done == count:
                print(f'{args.filename}: {done}/{count} chunks ({100 * done / count:.1f}%), {time.monotonic() - started:.0f}s', flush=True)
    with target.open('wb') as destination:
        for i in range(count):
            with (parts / f'{i:05d}.part').open('rb') as source:
                shutil.copyfileobj(source, destination)
    if target.stat().st_size != args.size or checksum(target) != args.sha256:
        raise RuntimeError('Publisher checksum mismatch; do not install this file')
    print(f'SHA-256 verified: {target}', flush=True)
