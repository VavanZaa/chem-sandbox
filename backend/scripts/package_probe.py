import json
import urllib.request

with urllib.request.urlopen('https://pypi.org/pypi/rdkit/2026.3.6/json', timeout=15) as response:
    metadata = json.load(response)
for item in metadata['urls']:
    if 'cp312-cp312-win_amd64' in item['filename']:
        print(item['url'])
        print(item['digests']['sha256'])
