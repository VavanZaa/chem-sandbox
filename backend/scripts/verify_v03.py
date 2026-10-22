"""Small real-compute integration check. Does not claim chemical accuracy."""
import json
from pathlib import Path
import sys
import time
import httpx

client = httpx.Client(base_url='http://127.0.0.1:8000/api', timeout=30)
print(json.dumps(client.get('/compute').json(), indent=2), flush=True)
payload = {
    'experiment_id': 'general_integration_v03', 'mode': 'general',
    'molecules': [{'name': 'Methanol', 'smiles': 'CO', 'amount_mol': .05}, {'name': 'Water', 'smiles': 'O', 'amount_mol': .1}],
    'compute': {'engine': sys.argv[1] if len(sys.argv) > 1 else 'mace_mp', 'device': sys.argv[2] if len(sys.argv) > 2 else 'auto'},
    'search': {'runs': 1, 'max_evaluations_per_run': 6, 'max_candidates_per_state': 6,
               'max_depth': 2, 'optimization_steps': 5, 'candidate_optimization_steps': 3},
}
response = client.post('/experiments', json=payload)
response.raise_for_status()
job = response.json()
print('SEARCH JOB', job['id'], flush=True)
for _ in range(240):
    job = client.get('/experiments/' + job['id']).json()
    if job['status'] in ['completed', 'failed']:
        break
    time.sleep(1)
if job['status'] != 'completed':
    raise RuntimeError(str(job.get('error') or job.get('message')))
result = job['result']
assert result['statistics']['accounted_runs'] == 1
assert result['stoichiometry'] is None
assert result['model']['device'] == payload['compute']['device'] or payload['compute']['device'] == 'auto'
print('SEARCH', json.dumps({'model': result['model'], 'statistics': result['statistics'], 'runtime_seconds': result['runtime_seconds']}), flush=True)
if payload['compute']['engine'] == 'all_compatible':
    assert result['discovery_comparison']['successful_engine_count'] >= 1
    assert result['discovery_comparison']['calibrated_confidence'] is None
    print('DISCOVERY COMPARISON', json.dumps(result['discovery_comparison']['engine_runs']), flush=True)
for node in result['graph']['nodes']:
    filename = Path(node['structure_reference']).name
    geometry = client.get(f'/experiments/{job["id"]}/structures/{filename}')
    geometry.raise_for_status()
    assert len(geometry.json()['atoms']) == 9
if result['graph']['edges']:
    response = client.post('/experiments/' + job['id'] + '/refinement', json={
        'engines': ([row['engine'] for row in result['discovery_comparison']['engine_runs']]
                    if len(sys.argv) > 3 and sys.argv[3] == 'all'
                    else ['mace_mp', 'ani2x'] if len(sys.argv) > 3 and sys.argv[3] == 'both'
                    else [payload['compute']['engine']]), 'device': payload['compute']['device'],
        'images': 3, 'steps': 2, 'endpoint_steps': 2, 'max_edges': 1})
    response.raise_for_status()
    refinement = response.json()
    print('REFINEMENT JOB', refinement['id'], flush=True)
    for _ in range(240):
        refinement = client.get('/experiments/' + refinement['id']).json()
        if refinement['status'] in ['completed', 'failed']:
            break
        time.sleep(1)
    assert refinement['status'] == 'completed', refinement
    print('REFINEMENT', json.dumps(refinement['result']['records']), flush=True)
    assert refinement['result']['calibrated_product_probabilities'] is None
