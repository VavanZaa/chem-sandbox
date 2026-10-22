import { useEffect, useState } from 'react'
import { requestJson } from '../lib/api'
import type { ExperimentResult, PathPopulationResult, RefinementJob } from '../lib/types'
import { downloadJson } from './AnalysisPanel'

export function EnginePanel({ result, jobId }: { result: ExperimentResult; jobId: string }) {
  const compatible = result.discovery_comparison?.engine_runs.map((run) => run.engine)
    ?? (result.experiment?.mode === 'general' ? ['mace_mp', 'ani2x', 'aimnet2', 'aimnet2_rxn'] : ['mace_mp'])
  const [engines, setEngines] = useState<string[]>(compatible)
  const [device, setDevice] = useState('auto')
  const [images, setImages] = useState(5)
  const [steps, setSteps] = useState(80)
  const [endpointSteps, setEndpointSteps] = useState(120)
  const [maxEdges, setMaxEdges] = useState(8)
  const [job, setJob] = useState<RefinementJob | null>(null)
  const [submitting, setSubmitting] = useState(false)
  const [error, setError] = useState('')
  const [ack, setAck] = useState(false)
  const [time, setTime] = useState(.001)
  const [temperature, setTemperature] = useState(result.experiment?.environment.temperature_K ?? 298.15)
  const initial = result.graph.nodes.find((node) => node.starting_runs)?.id ?? result.graph.nodes[0]?.id
  const [target, setTarget] = useState(result.graph.nodes.at(-1)?.id ?? '')
  const [populations, setPopulations] = useState<PathPopulationResult | null>(null)
  const busy = submitting || job?.status === 'queued' || job?.status === 'running'
  useEffect(() => {
    if (!job || !['queued', 'running'].includes(job.status)) return
    const timer = window.setInterval(() => requestJson<RefinementJob>(`/experiments/${job.id}`).then(setJob).catch((e) => setError(String(e))), 1500)
    return () => window.clearInterval(timer)
  }, [job?.id, job?.status])
  async function refine() {
    setSubmitting(true); setError(''); setPopulations(null); setJob(null)
    try { setJob(await requestJson<RefinementJob>(`/experiments/${jobId}/refinement`, { engines, device, images, steps, endpoint_steps: endpointSteps, max_edges: maxEdges })) }
    catch (e) { setError(String(e)) } finally { setSubmitting(false) }
  }
  async function predict() {
    setError(''); setPopulations(null); setSubmitting(true)
    try { setPopulations(await requestJson<PathPopulationResult>(`/refinements/${job!.id}/populations`, { temperature_K: temperature, duration_s: time, initial_state_id: initial, target_state_id: target, acknowledge_electronic_barrier_proxy: true })) }
    catch (e) { setError(String(e)) } finally { setSubmitting(false) }
  }
  return <div className="analysis-panel"><span className="eyebrow">MULTI-ENGINE DISCOVERY + MATCHED-PATH COMPARISON</span><h2>Compare discovery first. Then compute barriers.</h2>
    <p>Compatible engines perform independent bounded searches. Path refinement then evaluates the displayed graph's identical saved atom mappings with the selected models using climbing-image nudged elastic band (CI-NEB).</p>
    <div className="analysis-warning">Electronic barriers are not Gibbs free-energy barriers. No solvent, entropy, zero-point correction or transition-state vibrational validation is included. Failed paths stay unknown; agreement between models is not a calibrated product probability.</div>
    {result.discovery_comparison && <><div className="analysis-title"><h3>Independent discovery comparison</h3><span>{result.discovery_comparison.successful_engine_count} successful models</span></div>
      <div className="analysis-table-wrap"><table><thead><tr><th>Engine</th><th>Device</th><th>States</th><th>Paths</th><th>Runtime</th></tr></thead><tbody>{result.discovery_comparison.engine_runs.map((run) => <tr key={run.engine}><td>{run.model.model}</td><td>{run.model.device.toUpperCase()}</td><td>{run.states}</td><td>{run.paths}</td><td>{run.runtime_seconds.toFixed(1)} s</td></tr>)}</tbody></table></div>
      {(result.discovery_comparison.failed_engines.length > 0 || result.discovery_comparison.excluded_engines.length > 0) && <ul>{result.discovery_comparison.failed_engines.map((row) => <li key={row.engine}>{row.engine}: failed — {row.reason}</li>)}{result.discovery_comparison.excluded_engines.map((row) => <li key={row.engine}>{row.engine}: excluded — {row.reason}</li>)}</ul>}
      <h3>Connectivity agreement</h3><p>Agreement is model discovery frequency under this finite search budget, not chemical confidence.</p>
      <div className="analysis-table-wrap"><table><thead><tr><th>State</th><th>Models</th><th>Agreement</th></tr></thead><tbody>{result.discovery_comparison.states.map((row) => <tr key={row.comparison_key}><td>{row.labels.join(' / ')}</td><td>{row.engines.join(', ')}</td><td>{row.model_count}/{result.discovery_comparison!.successful_engine_count} ({(row.agreement_fraction * 100).toFixed(0)}%)</td></tr>)}</tbody></table></div>
      <h3>Pathway agreement</h3><div className="analysis-table-wrap"><table><thead><tr><th>Path</th><th>Models</th><th>Agreement</th><th>Within-model ΔE / eV</th></tr></thead><tbody>{result.discovery_comparison.paths.map((row) => <tr key={row.comparison_key}><td>{row.labels.join(' / ')}</td><td>{row.engines.join(', ')}</td><td>{row.model_count}/{result.discovery_comparison!.successful_engine_count}</td><td>{Object.entries(row.delta_energy_ev_by_engine).map(([engine, energy]) => `${engine}: ${energy.toFixed(3)}`).join(' · ')}</td></tr>)}</tbody></table></div>
      <ul>{result.discovery_comparison.warnings.map((warning) => <li key={warning}>{warning}</li>)}</ul></>}
    <fieldset disabled={busy}><div className="kinetics-controls">
      <div><span>Path engines</span>{compatible.map((engine) => <label className="acknowledge" key={engine}><input type="checkbox" checked={engines.includes(engine)} onChange={(e) => setEngines(e.target.checked ? [...engines, engine] : engines.filter((item) => item !== engine))} />{engine}</label>)}</div>
      <label>Device<select value={device} onChange={(e) => setDevice(e.target.value)}><option value="auto">Auto · prefer GPU</option><option value="cuda">Require CUDA GPU</option><option value="cpu">CPU</option></select></label>
      <label>Path images<input type="number" min="3" max="9" value={images} onChange={(e) => setImages(Number(e.target.value))} /></label>
      <label>NEB step budget<input type="number" min="1" max="300" value={steps} onChange={(e) => setSteps(Number(e.target.value))} /></label>
      <label>Endpoint step budget<input type="number" min="1" max="500" value={endpointSteps} onChange={(e) => setEndpointSteps(Number(e.target.value))} /></label>
      <label>Maximum edges<input type="number" min="1" max="24" value={maxEdges} onChange={(e) => setMaxEdges(Number(e.target.value))} /></label>
    </div><p>Models run sequentially to limit GPU memory use. Domain compatibility is checked before discovery; distorted path images can still leave a model's training domain.</p>
      <button className="primary-action" disabled={!result.graph.edges.length || !engines.length || busy} onClick={refine}>{busy ? 'Calculation in progress…' : 'Refine paths & compare engines'}</button></fieldset>
    {!result.graph.edges.length && <p>No candidate edges were discovered; no barrier calculation is available.</p>}
    {job && <p role="status">{job.status} · {job.progress}% · {job.message}</p>}
    {(error || job?.error) && <p className="analysis-warning" role="alert">{error || job?.error}</p>}
    {job?.result && <><div className="analysis-title"><h3>Path results · {job.result.unrefined_edges} unrefined edges</h3><button className="secondary-button" onClick={() => downloadJson(job.result, 'engine-comparison.json')}>Export comparison</button></div>
      <div className="analysis-table-wrap"><table><thead><tr><th>Path / engine</th><th>Device</th><th>ΔE / eV</th><th>ΔE‡ / eV</th><th>Validation</th></tr></thead><tbody>{job.result.records.map((row) => <tr key={`${row.engine}_${row.edge_id}`}><td>{row.edge_id}<br />{row.engine}</td><td>{row.metadata?.device ?? '—'}</td><td>{row.delta_energy_ev?.toFixed(4) ?? '—'}</td><td>{row.electronic_barrier_ev?.toFixed(4) ?? 'Unknown'}</td><td>{row.status.replaceAll('_', ' ')}<small className="record-note">{row.error ?? row.note}</small></td></tr>)}</tbody></table></div>
      <ul>{job.result.warnings.map((warning) => <li key={warning}>{warning}</li>)}</ul>
      <h3>Exploratory cluster-state populations</h3><p>Only available if every directed path is converged for an engine. This approximation uses ΔG‡ ≈ ΔE‡; it is not bulk efficiency. For molar yields, use Efficiency with condition-appropriate evidence.</p>
      <fieldset disabled={busy}><div className="kinetics-controls"><label>Target state<select value={target} onChange={(e) => { setTarget(e.target.value); setPopulations(null) }}>{result.graph.nodes.map((node) => <option key={node.id} value={node.id}>{node.label}</option>)}</select></label>
        <label>Temperature / K<input type="number" min="250" max="450" value={temperature} onChange={(e) => { setTemperature(Number(e.target.value)); setPopulations(null) }} /></label>
        <label>Time / s<input type="number" step="any" min="0.000000001" value={time} onChange={(e) => { setTime(Number(e.target.value)); setPopulations(null) }} /></label></div>
        <label className="acknowledge"><input type="checkbox" checked={ack} onChange={(e) => { setAck(e.target.checked); setPopulations(null) }} />I accept the unvalidated electronic-barrier proxy and complete cluster-network assumption.</label>
        <button className="primary-action" disabled={!ack || busy} onClick={predict}>Compare conditional state populations</button></fieldset>
      {populations && <><p className="analysis-warning">{populations.warning}</p>{populations.predictions.map((row) => <p key={row.engine}>{row.engine}: initial-state depletion {(row.state_conversion * 100).toFixed(3)}% · target population {(row.target_state_population * 100).toFixed(3)}%</p>)}{populations.unavailable.map((row) => <p key={row.engine}>{row.engine}: unavailable — {row.reason}</p>)}</>}
    </>}
  </div>
}
