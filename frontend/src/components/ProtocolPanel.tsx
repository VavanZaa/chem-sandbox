import { useState } from 'react'
import { requestJson } from '../lib/api'
import type { ProtocolResult } from '../lib/types'
import { downloadJson } from './AnalysisPanel'

const example = {
  id: 'staged_hexose_cycle_hypothetical',
  species: [
    { id: 'glucose', name: 'Glucose stereoisomer', smiles: 'OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O' },
    { id: 'fructose', name: 'Open-chain fructose', smiles: 'OCC(=O)C(O)C(O)C(O)CO' },
    { id: 'mannitol', name: 'Mannitol', smiles: 'OCC(O)C(O)C(O)C(O)CO' },
    { id: 'sorbose', name: 'Sorbose identity proxy', smiles: 'OCC(=O)C(O)C(O)C(O)CO' },
    { id: 'h2', name: 'Hydrogen', smiles: '[H][H]' },
  ],
  reactions: [
    { id: 'isomerize_1', name: 'Declared glucose to fructose isomerization', stoichiometry: { glucose: -1, fructose: 1 },
      rate_law: { kind: 'mass_action', rate_constant: { value: 0.08, unit: '1/s' }, source: 'hypothetical UI demonstration', basis: 'hypothetical' } },
    { id: 'hydrogenate', name: 'Declared fructose hydrogenation', stoichiometry: { fructose: -1, h2: -1, mannitol: 1 },
      rate_law: { kind: 'mass_action', rate_constant: { value: 2, unit: 'L/mol/s' }, source: 'hypothetical UI demonstration', basis: 'hypothetical' } },
    { id: 'dehydrogenate', name: 'Declared mannitol dehydrogenation', stoichiometry: { mannitol: -1, sorbose: 1, h2: 1 },
      rate_law: { kind: 'mass_action', rate_constant: { value: 0.04, unit: '1/s' }, source: 'hypothetical UI demonstration', basis: 'hypothetical' } },
    { id: 'isomerize_2', name: 'Declared sorbose to glucose isomerization', stoichiometry: { sorbose: -1, glucose: 1 },
      rate_law: { kind: 'mass_action', rate_constant: { value: 0.02, unit: '1/s' }, source: 'hypothetical UI demonstration', basis: 'hypothetical' } },
  ],
  initial_amounts: { glucose: { value: 0.01, unit: 'mol' } },
  initial_volume: { value: 0.1, unit: 'L' }, initial_temperature: { value: 298.15, unit: 'K' },
  stages: [
    { id: 'convert_1', name: 'First conversion', operation: 'reaction', duration: { value: 8, unit: 's' }, active_reactions: ['isomerize_1'] },
    { id: 'mid_addition', name: 'Add hydrogen after stage one', operation: 'addition', additions: [{ species_id: 'h2', amount: { value: 0.02, unit: 'mol' } }] },
    { id: 'convert_2', name: 'Hydrogenation stage', operation: 'reaction', duration: { value: 12, unit: 's' }, active_reactions: ['hydrogenate'] },
    { id: 'sample_1', name: 'Remove analytical sample', operation: 'sample', withdraw_fraction: 0.05 },
    { id: 'convert_3', name: 'Dehydrogenation stage', operation: 'reaction', duration: { value: 6, unit: 's' }, active_reactions: ['dehydrogenate'] },
    { id: 'hold_1', name: 'Timed hold', operation: 'hold', duration: { value: 3, unit: 's' } },
    { id: 'convert_4', name: 'Close declared cycle', operation: 'reaction', duration: { value: 5, unit: 's' }, active_reactions: ['isomerize_2'] },
  ],
  points_per_reaction_stage: 21, include_sbml: true, cross_check_roadrunner: true,
  acknowledge_evidence_limits: true,
}

function nonzero(values: Record<string, number>) {
  return Object.entries(values).filter(([, value]) => Math.abs(value) > 1e-12)
    .map(([key, value]) => `${key} ${value.toPrecision(3)} mol`).join(' · ') || 'none'
}

function downloadText(value: string, filename: string) {
  const url = URL.createObjectURL(new Blob([value], { type: 'application/sbml+xml' }))
  const anchor = document.createElement('a'); anchor.href = url; anchor.download = filename; anchor.click()
  URL.revokeObjectURL(url)
}

export function ProtocolPanel() {
  const [definition, setDefinition] = useState(JSON.stringify(example, null, 2))
  const [output, setOutput] = useState<ProtocolResult | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)

  async function run() {
    setBusy(true); setError(null); setOutput(null)
    try { setOutput(await requestJson<ProtocolResult>('/protocols/simulate', JSON.parse(definition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Protocol simulation failed') }
    finally { setBusy(false) }
  }

  return <div className="analysis-panel protocol-panel">
    <span className="eyebrow">SEQUENTIAL LAB PROTOCOL · USER-DECLARED EVIDENCE</span>
    <h2>Run additions and reactions in a controlled order.</h2>
    <div className="analysis-warning">The loaded hexose cycle is a hypothetical execution test, not a validated glucose process. Replace every channel and rate with sourced evidence. A protocol never feeds its expected products back into blind discovery.</div>
    <textarea className="network-editor" aria-label="Sequential protocol JSON" spellCheck={false} value={definition}
      onChange={(event) => { setDefinition(event.target.value); setOutput(null); setError(null) }} />
    <div className="button-row"><button className="primary-action" disabled={busy} onClick={run}>{busy ? 'Integrating stages…' : 'Run full staged protocol'}</button>
      <button className="secondary-button" disabled={busy} onClick={() => { setDefinition(JSON.stringify(example, null, 2)); setOutput(null) }}>Reset example</button></div>
    {error && <p className="analysis-warning" role="alert">{error}</p>}
    {output && <section>
      <div className="analysis-title"><h3>{output.protocol_id} · {output.total_duration_s.toPrecision(4)} s</h3><div className="button-row">
        <button className="secondary-button" onClick={() => downloadJson(output, 'lab-protocol-result.json')}>Export result</button>
        {output.sbml.document && <button className="secondary-button" onClick={() => downloadText(output.sbml.document!, 'lab-protocol-network.xml')}>Export SBML network</button>}
      </div></div>
      <div className="summary-cards">
        <div><span>Material ledger</span><strong>{output.material_ledger.closed ? 'Closed' : 'Failed'}</strong></div>
        <div><span>SBML validation</span><strong>{output.sbml.validation.status}</strong><small>libSBML {output.sbml.validation.libsbml_version}</small></div>
        <div><span>Independent solver</span><strong>{output.roadrunner_cross_check.status}</strong><small>{output.roadrunner_cross_check.engine}</small></div>
        <div><span>Stages</span><strong>{output.stages.length}</strong></div>
      </div>
      <h3>Stage ledger</h3><div className="analysis-table-wrap"><table><thead><tr><th>Stage</th><th>Operation / time</th><th>Active reactions</th><th>Added</th><th>Removed</th></tr></thead><tbody>
        {output.stages.map((stage) => <tr key={stage.id}><td>{stage.name}<small className="record-note">{stage.id}</small></td><td>{stage.operation}<br />{stage.start_time_s.toPrecision(4)} → {stage.end_time_s.toPrecision(4)} s</td><td>{stage.active_reactions.join(', ') || 'none'}</td><td>{nonzero(stage.additions_mol)}</td><td>{nonzero(stage.removed_mol)}</td></tr>)}
      </tbody></table></div>
      <h3>Final inventory</h3><div className="summary-cards">{Object.entries(output.final_amounts_mol).map(([id, value]) => <div key={id}><span>{id}</span><strong>{value.toPrecision(5)}</strong><small>mol</small></div>)}</div>
      {output.species_identity.filter((item) => item.stereo_status !== 'achiral').map((item) => <p key={item.id}>{item.name}: {item.stereo_status} · {item.stereochemistry.join('; ') || item.isomeric_smiles}</p>)}
      {output.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
    </section>}
  </div>
}
