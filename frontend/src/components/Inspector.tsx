import { lazy, Suspense } from 'react'
import { Activity, ArrowRight, Atom, CircleAlert, FlaskConical, Info, Route } from 'lucide-react'
import type { ExperimentResult, Selection } from '../lib/types'

const MoleculeViewer = lazy(() =>
  import('./MoleculeViewer').then((module) => ({ default: module.MoleculeViewer })),
)

const formatEnergy = (value: number) => `${value >= 0 ? '+' : ''}${value.toFixed(3)} eV`

interface Props { selection: Selection | null; result: ExperimentResult | null; jobId: string | null }

export function Inspector({ selection, result, jobId }: Props) {
  if (!selection || !result) {
    return (
      <aside className="inspector empty-inspector">
        <div className="inspector-orbit"><Atom size={28} /></div>
        <h3>Inspect the network</h3>
        <p>Select a state or transformation to reveal the calculation evidence behind it.</p>
      </aside>
    )
  }

  if (selection.kind === 'node') {
    const node = selection.value
    return (
      <aside className="inspector">
        <span className={`inspector-tag ${node.role}`}>{node.role === 'endpoint' ? 'search endpoint — not proven stable' : node.role}</span>
        <div className="inspector-heading"><div className="formula-mark"><FlaskConical size={20} /></div><div><span>SPECIES STATE</span><h2>{node.label}</h2></div></div>
        <Suspense fallback={<div className="viewer-loading-shell">Loading 3D viewer…</div>}>
          <MoleculeViewer jobId={jobId} structureReference={node.structure_reference} label={node.label} />
        </Suspense>
        <div className="species-interpretation">
          <div className="species-title"><span>SPECIES INTERPRETATION</span><small>identity metadata*</small></div>
          {node.species.map((species, index) => (
            <div className={`species-row ${species.resolved ? '' : 'unresolved'}`} key={`${species.formula}-${index}`}>
              <strong>{species.display_formula}</strong>
              <span>{species.name}</span>
              {species.aromatic === true && <small>Aromatic input template retained</small>}
              {species.stereo_status === 'specified' && <small>Stereochemistry: {species.stereo_assignments.join('; ') || species.stereochemistry}</small>}
              {species.stereo_status === 'unassigned' && <small>Stereochemistry: unspecified in input</small>}
              {species.isomeric_smiles && <small>Isomeric SMILES: {species.isomeric_smiles}</small>}
              {species.note && <small>{species.note}</small>}
            </div>
          ))}
        </div>
        <div className="evidence-card primary"><span>Relative energy</span><strong>{node.relative_energy_ev.toFixed(3)} <small>eV</small></strong><p>Relative to the lowest sampled state</p></div>
        <div className="metric-triple">
          <div><span>State appeared</span><strong>{node.runs_seen}/{result.statistics.search_runs}</strong></div>
          <div><span>Branch stopped</span><strong>{node.terminal_runs}/{result.statistics.search_runs}</strong></div>
          <div><span>Observations</span><strong>{node.occurrences}</strong></div>
        </div>
        <div className="detail-list">
          {Object.entries(node.endpoint_reasons).map(([reason, count]) => <div key={reason}><span>{reason.replaceAll('_', ' ')}</span><b>{count} branch stops</b></div>)}
          <div><span>Raw model energy</span><b>{node.energy_ev.toFixed(3)} eV</b></div>
          <div><span>Optimization</span><b>{node.converged ? 'Converged' : 'Step limit reached'}</b></div>
          <div><span>Max residual force</span><b>{node.max_force == null ? 'Unavailable' : `${node.max_force.toFixed(3)} eV/Å`}</b></div>
          <div><span>Structure log</span><b>{node.structure_reference ? 'Saved' : 'Unavailable'}</b></div>
        </div>
        <div className="honesty-note"><Info size={15} /><p>*For intact input-derived molecules, aromatic bond order comes from the mapped SMILES template and R/S or E/Z is reassigned from 3D coordinates. Changed or unmapped fragments remain unknown. Charges are metadata, not neural-model outputs. A converged geometry does not prove product stability.</p></div>
      </aside>
    )
  }

  const edge = selection.value
  const source = result.graph.nodes.find((node) => node.id === edge.source)
  const target = result.graph.nodes.find((node) => node.id === edge.target)
  return (
    <aside className="inspector">
      <span className="inspector-tag pathway">candidate pathway</span>
      <div className="inspector-heading"><div className="formula-mark"><Route size={20} /></div><div><span>TRANSFORMATION</span><h3 className="reaction-title">{source?.label}<ArrowRight size={16} />{target?.label}</h3></div></div>
      <div className="evidence-card primary"><span>Observed pathway frequency</span><strong>{Math.round(edge.run_frequency * 100)}<small>% of runs</small></strong><p>Present in {edge.run_count}/{result.statistics.search_runs} runs · {edge.transition_count} recorded transitions</p></div>
      <div className="metric-pair"><div><span>Δ energy</span><strong>{formatEnergy(edge.delta_energy_ev)}</strong></div><div><span>Barrier</span><strong>Not calculated</strong></div></div>
      <div className="detail-list">
        <div><span>Potential-energy change only</span><b>{edge.delta_energy_ev < -0.05 ? 'Lower sampled energy' : edge.delta_energy_ev > 0.05 ? 'Higher sampled energy' : 'Approximately neutral'}</b></div>
        <div><span>Proposal method</span><b>{edge.method.replaceAll('_', ' ')}</b></div>
        <div><span>Bonds created</span><b>{edge.bonds_created.length}</b></div>
        <div><span>Bonds removed</span><b>{edge.bonds_removed.length}</b></div>
        <div><span>Outgoing search share</span><b>{Math.round(edge.search_frequency * 100)}%</b></div>
      </div>
      <div className="honesty-note warning"><CircleAlert size={15} /><p>Discovery does not calculate barriers. Open Engines to attempt separate NEB refinement. ΔE is the mean of paired observations, not the difference between node minima or an activation barrier. Frequency is not reaction probability.</p></div>
    </aside>
  )
}
