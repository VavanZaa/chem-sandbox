"""Read the live size of an open download (directory listings can be stale)."""
import os
import sys
from pathlib import Path
cpu = len(sys.argv) > 1 and sys.argv[1] == 'cpu'
filename = 'torch-2.11.0+cpu-cp312-cp312-win_amd64.whl' if cpu else 'torch-2.11.0+cu128-cp312-cp312-win_amd64.whl'
total = 114469197 if cpu else 2753189216
target = Path(__file__).resolve().parents[2] / '.tools' / filename
with target.open('rb') as stream:
    size = os.fstat(stream.fileno()).st_size
print(f'{size:,} / {total:,} bytes ({100 * size / total:.2f}%)')
