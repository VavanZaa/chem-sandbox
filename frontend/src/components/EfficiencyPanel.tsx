import { useState } from 'react'
import { requestJson } from '../lib/api'
import type { EfficiencyResult, ExperimentDefinition } from '../lib/types'
import { downloadJson } from './AnalysisPanel'

const demonstration = {
  species: [
    { id: 'A', name: 'Ethanol', smiles: 'CCO', initial_mol: .05 },
    { id: 'B', name: 'Dimethyl ether', smiles: 'COC', initial_mol: 0 },
    { id: 'C', name: 'Ethene', smiles: 'C=C', initial_mol: 0 },
    { id: 'W', name: 'Water', smiles: 'O', initial_mol: 0 },
  ], reference_reactant: 'A', target_product: 'B', product_pool: ['B', 'C'],
  models: [{ id: 'Hypothetical competing channels', basis: 'hypothetical', channels: [
    { id: 'A_to_B', reactants: { A: 1 }, products: { B: 1 }, activation_free_energy_kj_mol: 80,
      reference_temperature_K: 298.15, activation_entropy_j_mol_K: 0, source: 'Invented demonstration barrier, NOT a prediction for ethanol chemistry' },
    { id: 'A_to_C_W', reactants: { A: 1 }, products: { C: 1, W: 1 }, activation_free_energy_kj_mol: 84,
      reference_temperature_K: 298.15, activation_entropy_j_mol_K: 0, source: 'Invented demonstration barrier, NOT a prediction for ethanol chemistry' },
  ] }],
}
const percentage = (value: number | null) => value === null ? 'Unavailable' : `${(100 * value).toFixed(3)}%`

export function EfficiencyPanel({ experiment }: { experiment: ExperimentDefinition }) {
  const [network, setNetwork] = useState(JSON.stringify(demonstration, null, 2))
  const [temperature, setTemperature] = useState(experiment.environment.temperature_K)
  const [pressure, setPressure] = useState(experiment.environment.pressure_atm)
  const [time, setTime] = useState(60)
  const [volume, setVolume] = useState(.05)
  const [phase, setPhase] = useState('solution')
  const [solvent, setSolvent] = useState<string>(experiment.environment.solvent)
  const [sensitivity, setSensitivity] = useState(5)
  const [ack, setAck] = useState(false)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const [output, setOutput] = useState<EfficiencyResult | null>(null)
  const clear = () => { setOutput(null); setError('') }
  async function simulate() {
    setBusy(true); clear()
    try {
      const data = JSON.parse(network)
      setOutput(await requestJson<EfficiencyResult>('/efficiency', { ...data, temperature_K: temperature, pressure_atm: pressure,
        duration_s: time, volume_L: volume, phase, solvent, sensitivity_kj_mol: sensitivity, acknowledge_assumptions: true }))
    } catch (caught) { setError(caught instanceof Error ? caught.message : 'Calculation failed') }
    finally { setBusy(false) }
  }
  function workspaceTemplate() {
    const species = experiment.mode === 'general' ? experiment.molecules.map((item, index) => ({ id: `R${index + 1}`, name: item.name, smiles: item.smiles, initial_mol: item.amount_mol })) : [
      { id: 'R1', name: 'Baking soda', smiles: '[Na+].O=C(O)[O-]', initial_mol: .05952 },
      { id: 'R2', name: 'Acetic acid', smiles: 'CC(=O)O', initial_mol: .05 },
      { id: 'R3', name: 'Water', smiles: 'O', initial_mol: 2.767 },
    ]
    // The bicarbonate demo has editable units: do not silently replace its inputs with default amounts.
    if (experiment.mode === 'bicarbonate') experiment.substances.forEach((item, i) => {
      const mass = [84.0066, 60.052, 18.01528][i]
      species[i].initial_mol = item.unit === 'mol' ? item.amount : item.amount * (item.unit === 'mL' ? .997 : 1) / mass
    })
    setNetwork(JSON.stringify({ species, reference_reactant: 'R1', target_product: 'ADD_PRODUCT_ID', product_pool: ['ADD_PRODUCT_ID'],
      models: [{ id: 'Your condition-specific evidence', basis: 'external_calculation', channels: [] }] }, null, 2))
    setTemperature(experiment.environment.temperature_K); setPressure(experiment.environment.pressure_atm); setSolvent(experiment.environment.solvent); setAck(false); clear()
  }
  return <div className="analysis-panel efficiency-panel"><span className="eyebrow">BULK REACTION EFFICIENCY · EVIDENCE REQUIRED</span><h2>Conversion is a kinetic question.</h2>
    <p>Calculate amounts over time from a balanced, competing-channel network. Temperature changes Eyring rates; concentration changes mass-action rates. Product yield is distinct from reactant conversion and listed-product composition.</p>
    <div className="analysis-warning">The initial example uses invented barriers to demonstrate competing products. It is not your experiment and is not a validated prediction of ethanol chemistry. Replace it with condition-appropriate evidence. Engine votes cannot supply product probabilities.</div>
    <fieldset disabled={busy}>
      <div className="kinetics-controls">
        <label>Temperature / K<input aria-label="Efficiency temperature" type="number" min="200" max="1500" value={temperature} onChange={(e) => { setTemperature(Number(e.target.value)); clear() }} /></label>
        <label>Reaction time / s<input aria-label="Reaction time" type="number" min="0.000001" max="10000000" value={time} onChange={(e) => { setTime(Number(e.target.value)); clear() }} /></label>
        <label>Phase / vessel model<select value={phase} onChange={(e) => { setPhase(e.target.value); clear() }}><option value="solution">Solution · fixed volume</option><option value="ideal_gas">Ideal gas · constant pressure</option></select></label>
        {phase === 'solution' ? <label>Solution volume / L<input type="number" min="0.000001" step="any" value={volume} onChange={(e) => { setVolume(Number(e.target.value)); clear() }} /></label> : <label>Pressure / atm<input type="number" min="0.002" max="1000" step="any" value={pressure} onChange={(e) => { setPressure(Number(e.target.value)); clear() }} /></label>}
        <label>Solvent / evidence context<input value={solvent} onChange={(e) => { setSolvent(e.target.value); clear() }} /></label>
        <label>Barrier sensitivity / ± kJ mol⁻¹<input type="number" min="0" max="30" value={sensitivity} onChange={(e) => { setSensitivity(Number(e.target.value)); clear() }} /></label>
      </div>
      <p>Pressure affects concentrations only in the ideal-gas model. In solution, solvent, pH, catalysts and mixing must be represented by your barrier evidence; this tool does not infer those effects.</p>
      <div className="analysis-title"><h3>Species & kinetic evidence</h3><div className="button-row"><button className="secondary-button" onClick={workspaceTemplate}>Use workspace amounts</button><button className="secondary-button" onClick={() => { setNetwork(JSON.stringify(demonstration, null, 2)); setAck(false); clear() }}>Load hypothetical example</button></div></div>
      <p>Add all products (initial_mol: 0), competing and reverse channels. Each channel must conserve atoms and charge. Supply standard-state ΔG‡ in kJ/mol, its reference temperature and ΔS‡ in J/(mol·K). ΔS‡ = 0 assumes constant ΔG‡. Up to four evidence models can be compared on identical channels.</p>
      <textarea className="network-editor" aria-label="Balanced kinetic network JSON" spellCheck={false} value={network} onChange={(e) => { setNetwork(e.target.value); setAck(false); clear() }} />
      <label className="acknowledge"><input type="checkbox" checked={ack} onChange={(e) => { setAck(e.target.checked); clear() }} />I have checked the network, elementary mass-action assumption, evidence and standard states (1 M solution / 1 atm gas). I understand hypothetical inputs give hypothetical results.</label>
      <button className="primary-action" disabled={!ack || busy} onClick={simulate}>{busy ? 'Integrating kinetic models…' : 'Predict conditional efficiency'}</button>
    </fieldset>
    {error && <p className="analysis-warning" role="alert">{error}</p>}
    {output && <section><div className="analysis-title"><h3>Conditional prediction · not calibrated probability</h3><button className="secondary-button" onClick={() => downloadJson(output, 'reaction-efficiency.json')}>Export evidence & results</button></div>
      {output.models.map((model) => <section className="yield-section" key={model.model_id}><h3>{model.model_id} · {model.basis.replaceAll('_', ' ')}</h3>
        <div className="summary-cards"><div><span>Reference conversion</span><strong>{percentage(model.nominal.conversion)}</strong></div><div><span>Target yield / ceiling</span><strong>{percentage(model.nominal.target_yield)}</strong><small>{model.nominal.target_formed_mol.toPrecision(5)} / {model.stoichiometric_ceiling_mol.toPrecision(5)} mol</small></div><div><span>Target share of listed products</span><strong>{percentage(model.nominal.target_product_share)}</strong><small>Molar composition, not probability</small></div></div>
        <p>Barrier-scenario yield range: {model.yield_sensitivity_range?.map(percentage).join(' – ') ?? 'Unavailable'}. This is sensitivity, not a confidence interval.</p>
        <div className="analysis-table-wrap"><table><thead><tr><th>Barrier scenario</th><th>Conversion</th><th>Target yield</th><th>Target product share</th></tr></thead><tbody>{model.scenarios.map((s) => <tr key={s.scenario}><td>{s.scenario.replaceAll('_', ' ')}</td><td>{percentage(s.conversion)}</td><td>{percentage(s.target_yield)}</td><td>{percentage(s.target_product_share)}</td></tr>)}</tbody></table></div>
        <h4>Final amounts / mol</h4><div className="summary-cards">{Object.entries(model.nominal.final_moles).map(([id, value]) => <div key={id}><span>{id}</span><strong>{value.toPrecision(5)}</strong></div>)}</div>
      </section>)}
      <ul>{output.warnings.map((warning) => <li key={warning}>{warning}</li>)}</ul>
    </section>}
  </div>
}
