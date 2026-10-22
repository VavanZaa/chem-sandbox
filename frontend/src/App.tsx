import { useEffect, useMemo, useState } from 'react'
import { Activity, AlertTriangle, CircleDot, Cpu, FlaskConical, Gauge, Hexagon, Search, Sparkles } from 'lucide-react'
import { ExperimentPanel } from './components/ExperimentPanel'
import { Inspector } from './components/Inspector'
import { ReactionGraph } from './components/ReactionGraph'
import { OutcomesPanel, KineticsPanel } from './components/AnalysisPanel'
import { EfficiencyPanel } from './components/EfficiencyPanel'
import { EnginePanel } from './components/EnginePanel'
import { RunProgress } from './components/RunProgress'
import { ProtocolPanel } from './components/ProtocolPanel'
import { OriginsPanel } from './components/OriginsPanel'
import { checkHealth, getExperiment, requestJson, startExperiment } from './lib/api'
import type { ComputeStatus, ExperimentDefinition, ExperimentJob, ExperimentResult, Selection } from './lib/types'

const createDefaultExperiment = (): ExperimentDefinition => ({
  experiment_id: 'molecular_sandbox_v04',
  mode: 'general',
  molecules: [{ name: 'Ethanol', smiles: 'CCO', amount_mol: .05, copies: 1 }, { name: 'Water', smiles: 'O', amount_mol: 2.767, copies: 1 }],
  compute: { engine: 'all_compatible', device: 'auto' },
  substances: [
    { id: 'sodium_bicarbonate', amount: 0.05952, unit: 'mol' },
    { id: 'acetic_acid', amount: 0.05, unit: 'mol' },
    { id: 'water', amount: 2.767, unit: 'mol' },
  ],
  environment: { temperature_K: 298.15, pressure_atm: 1, solvent: 'water', mixing: 'standard' },
  search: { effort: 'medium', runs: 8, seed: 42, max_candidates_per_state: 8, optimization_steps: 10,
    candidate_optimization_steps: 6, max_depth: 3, max_evaluations_per_run: 24, max_states_per_run: 16, frontier_width: 6 },
})

function App() {
  const [experiment, setExperiment] = useState(createDefaultExperiment)
  const [job, setJob] = useState<ExperimentJob | null>(null)
  const [result, setResult] = useState<ExperimentResult | null>(null)
  const [selection, setSelection] = useState<Selection | null>(null)
  const [minimumFrequency, setMinimumFrequency] = useState(0)
  const [view, setView] = useState<'network' | 'outcomes' | 'kinetics' | 'engines' | 'efficiency' | 'protocol' | 'origins'>('network')
  const [hardware, setHardware] = useState<ComputeStatus | null>(null)
  const [apiOnline, setApiOnline] = useState<boolean | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => { checkHealth().then(setApiOnline); requestJson<ComputeStatus>('/compute').then(setHardware).catch(() => {}) }, [])

  useEffect(() => {
    if (!job || !['queued', 'running'].includes(job.status)) return
    const timer = window.setInterval(async () => {
      try {
        const next = await getExperiment(job.id)
        setJob(next)
        if (['completed', 'failed'].includes(next.status) && next.result) {
          setResult(next.result)
          const co2 = next.result.graph.nodes.find((node) => node.formula.split('+').map((part) => part.trim()).includes('CO2'))
          const initial = co2 ?? next.result.graph.nodes.at(-1)
          if (initial) setSelection({ kind: 'node', value: initial })
        }
        if (next.status === 'failed') {
          setError(next.error ?? next.message)
        }
      } catch (caught) {
        setError(caught instanceof Error ? caught.message : 'Could not read experiment status')
      }
    }, 800)
    return () => window.clearInterval(timer)
  }, [job?.id, job?.status])

  const running = job?.status === 'queued' || job?.status === 'running'
  const co2DiscoveryFrequency = useMemo(
    () => result ? Math.round(result.statistics.co2_run_frequency * 100) : 0,
    [result],
  )

  const run = async () => {
    setError(null)
    setResult(null)
    setSelection(null)
    setView('network')
    try {
      const started = await startExperiment(experiment)
      setJob(started)
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : 'Could not start experiment')
    }
  }

  return (
    <div className="app-shell">
      <header className="topbar">
        <div className="brand"><div className="brand-mark"><Hexagon size={26} /><CircleDot size={10} /></div><div><strong>CHEM<span>SANDBOX</span></strong><small>EXPERIMENTAL MATTER LAB</small></div></div>
        <nav aria-label="Result views"><button className={view === 'network' ? 'active' : ''} onClick={() => setView('network')}>NETWORK</button><button className={view === 'outcomes' ? 'active' : ''} disabled={!result || running} onClick={() => setView('outcomes')}>OUTCOMES</button><button className={view === 'engines' ? 'active' : ''} disabled={!result || running} onClick={() => setView('engines')}>ENGINES</button><button className={view === 'kinetics' ? 'active' : ''} disabled={!result || running} onClick={() => setView('kinetics')}>STATE KINETICS</button><button className={view === 'efficiency' ? 'active' : ''} onClick={() => setView('efficiency')}>EFFICIENCY</button><button className={view === 'protocol' ? 'active' : ''} onClick={() => setView('protocol')}>PROTOCOL</button><button className={view === 'origins' ? 'active' : ''} onClick={() => setView('origins')}>ORIGINS</button></nav>
        <div className="system-status"><div><span className={`status-dot ${apiOnline ? 'online' : apiOnline === false ? 'offline' : ''}`} />{apiOnline ? 'SYSTEM READY' : apiOnline === false ? 'API OFFLINE' : 'CHECKING'}</div><span className="device-status" title={hardware?.reason}><Cpu size={15} /> {result ? `${result.model.engine} · ${result.model.device.toUpperCase()}` : hardware?.cuda_available ? 'CUDA AVAILABLE' : 'CPU'}</span></div>
      </header>

      <main className="workspace">
        <ExperimentPanel experiment={experiment} setExperiment={setExperiment} onRun={run} onReset={() => setExperiment(createDefaultExperiment())} running={Boolean(running)} />

        <section className="canvas-panel">
          <div className="canvas-header">
            <div><span className="eyebrow">{view === 'efficiency' ? 'SEPARATE BULK KINETIC MODEL' : view === 'protocol' ? 'SEQUENTIAL OPERATIONS' : view === 'origins' ? 'PREBIOTIC REACTOR DEVELOPMENT' : result ? `LAST SUBMITTED RUN · SEED ${result.search_settings.seed}` : 'LIVE WORKSPACE / 04'}</span><h1>{ { network: 'Candidate reaction network', outcomes: 'Search outcomes & yield ceiling', kinetics: 'Conditional state kinetics', engines: 'Multi-engine comparison', efficiency: 'Condition-dependent efficiency', protocol: 'Full staged laboratory protocol', origins: 'Miller–Urey apparatus laboratory' }[view] }</h1></div>
            {result && view === 'network' && <div className="graph-toolbar"><label>MIN. DISCOVERY <b>{minimumFrequency}%</b><input aria-label="Minimum pathway discovery" type="range" min="0" max="100" step="5" value={minimumFrequency} onChange={(event) => setMinimumFrequency(Number(event.target.value))} /></label><span>DRAG STATES TO MOVE</span><span>{result.graph.nodes.length} states</span><span>{result.graph.edges.length} paths</span></div>}
          </div>

          <div className="canvas-body">
            <div className="view-surface" hidden={view !== 'efficiency'}><EfficiencyPanel experiment={experiment} /></div>
            <div className="view-surface" hidden={view !== 'protocol'}><ProtocolPanel /></div>
            <div className="view-surface" hidden={view !== 'origins'}><OriginsPanel /></div>
            {running && job ? <div className="view-surface" hidden={view === 'efficiency' || view === 'protocol' || view === 'origins'}><RunProgress job={job} /></div> : result ? (
              <>
                <div className="view-surface" hidden={view !== 'network'}>
                <ReactionGraph key={job?.id} nodes={result.graph.nodes} edges={result.graph.edges} minimumFrequency={minimumFrequency} totalRuns={result.statistics.search_runs} onSelect={setSelection} />
                {result.graph.edges.length === 0 && (
                  <div className="no-pathway-note">
                    <AlertTriangle size={15} />
                    <div>
                      <strong>No accepted transformation in this search</strong>
                      <span>
                        {result.statistics.rejected_candidate_states} structures were rejected for unsupported connectivity. Check Outcomes & Yield for stopping reasons. Change the seed or search budgets to explore further; a missing pathway does not mean the real reaction cannot occur.
                      </span>
                    </div>
                  </div>
                )}
                </div>
                <div className="view-surface" hidden={view !== 'outcomes'}><OutcomesPanel result={result} onSelect={setSelection} /></div>
                <div className="view-surface" hidden={view !== 'kinetics'}><KineticsPanel key={job?.id} result={result} jobId={job!.id} /></div>
                <div className="view-surface" hidden={view !== 'engines'}><EnginePanel key={job?.id} result={result} jobId={job!.id} /></div>
              </>
            ) : (
              <div className="empty-canvas" hidden={view === 'efficiency' || view === 'protocol' || view === 'origins'}>
                <div className="constellation">
                  <span className="orb orb-a">SMILES</span><span className="orb orb-b">MACE / ANI / AIMNet</span><span className="orb orb-c">GPU</span>
                  <i className="line line-a" /><i className="line line-b" />
                  <div className="center-orb"><FlaskConical size={34} /></div>
                </div>
                <span className="eyebrow">READY TO EXPLORE</span>
                <h2>Mix matter. Search possibility.</h2>
                <p>Run a bounded atomistic search to reveal candidate states and transformations—without looking up a reaction answer.</p>
                <button onClick={run} disabled={apiOnline === false}><Sparkles size={16} /> START FIRST EXPERIMENT</button>
                <div className="method-chips"><span><Activity size={13} /> Energy + forces</span><span><Search size={13} /> Guided search</span><span><Gauge size={13} /> Honest uncertainty</span></div>
              </div>
            )}
            {error && <div className="error-banner"><AlertTriangle size={18} /><div><strong>Model failure / result uncertain</strong><span>{error}</span></div><button onClick={() => setError(null)}>×</button></div>}
          </div>

          {result && view !== 'protocol' && view !== 'efficiency' && view !== 'origins' && <footer className="result-strip">
            <div><span>STARTS ACCOUNTED</span><strong>{result.statistics.accounted_runs}/{result.statistics.search_runs}</strong></div>
            <div><span>VALID STATES</span><strong>{result.statistics.valid_optimized_states}</strong><small>{result.statistics.rejected_candidate_states} structures rejected</small></div>
            <div><span>UNIQUE STATES</span><strong>{result.statistics.unique_states}</strong></div>
            <div><span>CANDIDATE PATHS</span><strong>{result.statistics.candidate_pathways}</strong></div>
            {result.stoichiometry ? <div className="co2-stat"><span>CO₂ DISCOVERED IN</span><strong>{co2DiscoveryFrequency}%</strong><small>{result.statistics.co2_containing_runs}/{result.statistics.search_runs} starts · branch stopped with CO₂ in {result.statistics.co2_terminal_runs}/{result.statistics.search_runs}</small></div> : <div className="co2-stat"><span>DISCOVERY MODELS</span><strong>{result.discovery_comparison ? `${result.discovery_comparison.successful_engine_count} compared` : result.model.engine}</strong><small>{result.discovery_comparison ? `display: ${result.display_engine} · independent searches` : `${result.model.device.toUpperCase()} · ${result.runtime_seconds.toFixed(1)} s`}</small></div>}
            <p><AlertTriangle size={13} /> Search discovery and branch stops are not product yields or reaction probabilities.</p>
          </footer>}
        </section>

        <Inspector selection={selection} result={result} jobId={job?.id ?? null} />
      </main>
    </div>
  )
}

export default App
