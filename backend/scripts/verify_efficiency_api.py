"""Check actual HTTP efficiency requests using explicitly hypothetical evidence."""
import json
import math
import httpx

def barrier(k):
    return -8.31446261815324 * 298.15 * math.log(k * 6.62607015e-34 / (1.380649e-23 * 298.15)) / 1000

payload = {
    'species': [{'id': 'A', 'name': 'Ethanol', 'smiles': 'CCO', 'initial_mol': .05},
                {'id': 'B', 'name': 'Dimethyl ether', 'smiles': 'COC', 'initial_mol': 0},
                {'id': 'C', 'name': 'Ethene', 'smiles': 'C=C', 'initial_mol': 0},
                {'id': 'W', 'name': 'Water', 'smiles': 'O', 'initial_mol': 0}],
    'models': [{'id': 'Hypothetical analytic check', 'basis': 'hypothetical', 'channels': [
        {'id': 'ab', 'reactants': {'A': 1}, 'products': {'B': 1}, 'activation_free_energy_kj_mol': barrier(2), 'source': 'Invented k=2/s for analytic test, not real chemistry'},
        {'id': 'acw', 'reactants': {'A': 1}, 'products': {'C': 1, 'W': 1}, 'activation_free_energy_kj_mol': barrier(1), 'source': 'Invented k=1/s for analytic test, not real chemistry'}]}],
    'reference_reactant': 'A', 'target_product': 'B', 'product_pool': ['B', 'C'], 'duration_s': 1,
    'temperature_K': 298.15, 'sensitivity_kj_mol': 5, 'acknowledge_assumptions': True,
}
with httpx.Client(base_url='http://127.0.0.1:8000/api', timeout=30) as client:
    previous = None
    for temperature in [298.15, 330]:
        payload['temperature_K'] = temperature
        response = client.post('/efficiency', json=payload)
        response.raise_for_status()
        result = response.json()
        nominal = result['models'][0]['nominal']
        assert result['calibrated_product_probabilities'] is None
        if previous is not None:
            assert nominal['conversion'] > previous
        else:
            assert math.isclose(nominal['conversion'], 1 - math.exp(-3), rel_tol=1e-6)
            assert math.isclose(nominal['target_product_share'], 2 / 3, rel_tol=1e-6)
        previous = nominal['conversion']
        print(json.dumps({'temperature_K': temperature, **nominal}), flush=True)
    payload['models'][0]['channels'][1]['products'].pop('W')
    assert client.post('/efficiency', json=payload).status_code == 422
    print('Unbalanced channel rejected with HTTP 422', flush=True)
