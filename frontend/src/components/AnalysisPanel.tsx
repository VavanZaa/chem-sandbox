import { useState } from 'react'
import { runKinetics } from '../lib/api'
import type { ExperimentResult, KineticsRequest, KineticsResult, RateEvidence, Selection } from '../lib/types'

export function downloadJson(value: unknown, filename: string) {
  const url = URL.createObjectURL(new Blob([JSON.stringify(value, null, 2)], { type: 'application/json' }))
  const anchor = document.createElement('a')
  anchor.href = url; anchor.download = filename; anchor.click()
  window.setTimeout(() => URL.revokeObjectURL(url), 1000)
}

const readable = (text: string) => text.replaceAll('_', ' ')
const number = (value: number) => Number(value.toPrecision(5)).toLocaleString(undefined, { maximumSignificantDigits: 5 })

export function OutcomesPanel({ result, onSelect }: { result: ExperimentResult; onSelect: (selection: Selection) => void }) {
  const yieldData = result.stoichiometry
  return <div className="analysis-panel">
    <div className="analysis-title"><div><span className="eyebrow">SEARCH ACCOUNTING</span><h2>Every start, every stopping reason</h2></div>
      <button className="secondary-button" onClick={() => downloadJson(result, 'chem-sandbox-result.json')}>Export full result</button></div>
    <p>{result.statistics.accounted_runs}/{result.statistics.search_runs} starts accounted for. Each start explores multiple branches. Several states may be endpoints in the same start; their percentages must not be added as a product distribution.</p>
    <div className="summary-cards">{Object.entries(result.statistics.run_status_counts).map(([status, count]) => <div key={status}><span>{readable(status)}</span><strong>{count} starts</strong></div>)}</div>
    <p>“Explored” means the generated proposals were exhausted within the configured operator set—not that all possible chemistry was explored. No product is forced.</p>
    <div className="analysis-table-wrap"><table><thead><tr><th>Start</th><th>Status</th><th>Evaluations*</th><th>Expanded</th><th>Branch stops</th></tr></thead><tbody>
      {result.run_outcomes.map((run) => <tr key={run.run}><td>{run.run + 1}</td><td>{readable(run.status)}</td><td>{run.evaluations} + {run.initialization_evaluations}</td><td>{run.expanded_states}</td><td>
        {run.endpoints.map((endpoint, index) => { const node = result.graph.nodes.find((item) => item.id === endpoint.state_id)
          return <div className="endpoint-row" key={index}><button disabled={!node} onClick={() => node && onSelect({ kind: 'node', value: node })}>{node?.label ?? 'No valid state'}</button><span>{readable(endpoint.reason)} · depth {endpoint.depth}</span></div>
        })}
      </td></tr>)}
    </tbody></table></div>
    <small>*Candidate evaluations + initialization calculations. Rejected geometries remain saved in the local structure logs and full export.</small>
    <h3>Multiple-state discovery</h3>
    <div className="analysis-table-wrap"><table><thead><tr><th>State</th><th>Discovered in</th><th>Branch stopped in</th><th>Geometry converged?</th></tr></thead><tbody>
      {result.graph.nodes.map((node) => <tr key={node.id}><td><button className="text-button" onClick={() => onSelect({ kind: 'node', value: node })}>{node.label}</button></td><td>{node.runs_seen}/{result.statistics.search_runs}</td><td>{node.terminal_runs}/{result.statistics.search_runs}</td><td>{node.converged ? 'Yes — local geometry only' : 'No — step limit'}</td></tr>)}
    </tbody></table></div>
    {yieldData ? <section className="yield-section"><span className="eyebrow">SEPARATE BALANCED-REACTION CALCULATOR</span><h2>Theoretical yield ceiling</h2>
      <p className="equation">{yieldData.equation}</p>
      <p>From the submitted amounts, <strong>{yieldData.limiting_reagent}</strong> is limiting. These values do not come from search frequencies and are not an actual-yield prediction.</p>
      <div className="summary-cards"><div><span>Maximum CO₂</span><strong>{number(yieldData.co2_mol)} mol</strong><small>{number(yieldData.co2_mass_g)} g</small></div>
        <div><span>Ideal dry gas</span><strong>{number(yieldData.co2_ideal_gas_L)} L</strong><small>{yieldData.temperature_K} K · {yieldData.pressure_atm} atm</small></div>
        <div><span>Sodium acetate ceiling</span><strong>{number(yieldData.sodium_acetate_mol)} mol</strong></div></div>
      <p>Excess baking soda: {number(yieldData.excess_bicarbonate_mol)} mol · Excess acid: {number(yieldData.excess_acid_mol)} mol · Water formed: {number(yieldData.water_formed_mol)} mol</p>
      <ul>{yieldData.assumptions.map((note) => <li key={note}>{note}</li>)}</ul>
    </section> : <p>General molecules: open Efficiency and supply a balanced reaction network with kinetic evidence to calculate bulk conversion and yield. Search counts are not yields.</p>}
    <details><summary>Scientific assumptions & diagnostics ({result.errors.length})</summary><ul>{result.scientific_notes.map((note, index) => <li key={index}>{note}</li>)}</ul>
      {result.errors.slice(0, 50).map((error, index) => <p className="diagnostic" key={index}>Start {Number(error.run) + 1} · {String(error.name)}: {String(error.error)} {error.structure_reference ? `· saved: ${String(error.structure_reference)}` : ''}</p>)}
      {result.errors.length > 50 && <p>Remaining diagnostics are in the full JSON export.</p>}
    </details>
  </div>
}

export function KineticsPanel({ result, jobId }: { result: ExperimentResult; jobId: string }) {
  const [values, setValues] = useState<Record<string, string>>({})
  const [mode, setMode] = useState<'rate' | 'barrier'>('rate')
  const [source, setSource] = useState('')
  const [basis, setBasis] = useState<RateEvidence['basis']>('hypothetical')
  const [temperature, setTemperature] = useState(String(result.experiment?.environment.temperature_K ?? result.stoichiometry?.temperature_K ?? 298.15))
  const [duration, setDuration] = useState('10')
  const [initial, setInitial] = useState(result.graph.nodes.find((node) => node.starting_runs)?.id ?? '')
  const [acknowledged, setAcknowledged] = useState(false)
  const [output, setOutput] = useState<KineticsResult | null>(null)
  const [submitted, setSubmitted] = useState<KineticsRequest | null>(null)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const clear = () => { setOutput(null); setError('') }
  const ready = result.graph.edges.length > 0 && acknowledged && source.trim().length >= 3 && result.graph.edges.every((edge) => values[edge.id]?.trim() && Number.isFinite(Number(values[edge.id])) && Number(values[edge.id]) >= 0)
  async function simulate() {
    setBusy(true); setError(''); setOutput(null)
    const request: KineticsRequest = { initial_state_id: initial, temperature_K: Number(temperature), duration_s: Number(duration), acknowledge_conditional_model: true,
      rates: result.graph.edges.map((edge) => ({ edge_id: edge.id, temperature_K: Number(temperature), source: source.trim(), basis,
        ...(mode === 'rate' ? { rate_s_inv: Number(values[edge.id]) } : { activation_free_energy_ev: Number(values[edge.id]) }) })) }
    try { const response = await runKinetics(jobId, request); setSubmitted(request); setOutput(response) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Kinetics failed') }
    finally { setBusy(false) }
  }
  const colors = ['#85efbb', '#71cde4', '#ebbd65', '#c3a0ed', '#ee938b', '#b9d579']
  return <div className="analysis-panel kinetics-panel"><span className="eyebrow">EVIDENCE-GATED KINETICS</span><h2>What would these pathways do over time?</h2>
    <p>Supply a rate or a Gibbs activation barrier for <strong>every directed edge</strong>. No values are inferred from MACE ΔE or how often a pathway was found. Reverse paths require their own evidence.</p>
    <div className="analysis-warning">This is a conditional, effective first-order state model—not a bulk multi-reactant yield calculator. Concentration, solvent, pH and gas escape are not modeled. Hypothetical inputs produce hypothetical populations.</div>
    {!result.graph.edges.length ? <p>No pathways were found. Kinetics is unavailable for this search.</p> : <>
    <fieldset disabled={busy}>
      <div className="kinetics-controls">
        <label>Input type<select value={mode} onChange={(event) => { setMode(event.target.value as typeof mode); setValues({}); clear() }}><option value="rate">Effective rate k (s⁻¹)</option><option value="barrier">Activation free energy ΔG‡ (eV)</option></select></label>
        <label>Initial state<select value={initial} onChange={(event) => { setInitial(event.target.value); clear() }}>{result.graph.nodes.map((node) => <option key={node.id} value={node.id}>{node.label}</option>)}</select></label>
        <label>Temperature (K)<input type="number" min="250" max="450" value={temperature} onChange={(event) => { setTemperature(event.target.value); clear() }} /></label>
        <label>Duration (s)<input type="number" min="0.000000001" max="100000000" step="any" value={duration} onChange={(event) => { setDuration(event.target.value); clear() }} /></label>
        <label>Evidence basis<select value={basis} onChange={(event) => { setBasis(event.target.value as typeof basis); clear() }}><option value="hypothetical">Hypothetical / what-if</option><option value="experimental">Experimental (user supplied)</option><option value="external_calculation">External calculation (user supplied)</option></select></label>
        <label className="full-width">Source, method & conditions<input value={source} maxLength={1000} placeholder="Citation / calculation method / explicit hypothetical assumption" onChange={(event) => { setSource(event.target.value); clear() }} /></label>
      </div>
      <div className="analysis-table-wrap"><table><thead><tr><th>Directed pathway</th><th>{mode === 'rate' ? 'k / s⁻¹' : 'ΔG‡ / eV'}</th></tr></thead><tbody>{result.graph.edges.map((edge) => <tr key={edge.id}><td>{result.graph.nodes.find((node) => node.id === edge.source)?.label}<br /><span className="path-arrow">→</span> {result.graph.nodes.find((node) => node.id === edge.target)?.label}</td><td><input aria-label={`${edge.id} ${mode}`} type="number" min="0" step="any" placeholder="Required" value={values[edge.id] ?? ''} onChange={(event) => { setValues({ ...values, [edge.id]: event.target.value }); clear() }} /></td></tr>)}</tbody></table></div>
      <label className="acknowledge"><input type="checkbox" checked={acknowledged} onChange={(event) => { setAcknowledged(event.target.checked); clear() }} />I assume this network is complete, rates are effective first-order values at this temperature, and understand these are conditional populations, not validated yields.</label>
      <button className="primary-action" disabled={!ready || busy} onClick={simulate}>{busy ? 'Calculating…' : 'Calculate conditional populations'}</button>
    </fieldset>
    {mode === 'barrier' && <p>Uses <a href="https://goldbook.iupac.org/terms/view/T06470" target="_blank" rel="noreferrer">Eyring transition-state theory</a> with transmission coefficient 1. A potential-energy difference is not a free-energy barrier.</p>}
    {error && <p role="alert" className="analysis-warning">{error}</p>}
    {output && <section className="kinetics-output"><div className="analysis-title"><h3>Conditional state populations · {submitted?.rates[0]?.basis.replaceAll('_', ' ')}</h3><button className="secondary-button" onClick={() => downloadJson({ job_id: jobId, input: submitted, output }, 'conditional-kinetics.json')}>Export model & results</button></div>
      <svg className="population-chart" viewBox="0 0 640 210" role="img" aria-label="Conditional state populations over time"><path d="M40 10 V180 H620" stroke="#496157" fill="none" /><text x="3" y="17">100%</text><text x="15" y="180">0%</text><text x="40" y="203">0 s</text><text x="530" y="203">{submitted?.duration_s} s</text>
        {result.graph.nodes.map((node, index) => <polyline key={node.id} fill="none" stroke={colors[index % colors.length]} strokeWidth="2" points={output.samples.map((sample, i) => `${40 + i / (output.samples.length - 1) * 580},${180 - sample.fractions[node.id] * 170}`).join(' ')} />)}
      </svg>
      <table><thead><tr><th>State</th><th>Final fraction</th></tr></thead><tbody>{result.graph.nodes.map((node, index) => <tr key={node.id}><td><span style={{ color: colors[index % colors.length] }}>●</span> {node.label}</td><td>{(output.samples.at(-1)!.fractions[node.id] * 100).toFixed(3)}%</td></tr>)}</tbody></table>
      <ul>{output.warnings.map((warning) => <li key={warning}>{warning}</li>)}</ul>
    </section>}
    </>}
  </div>
}
