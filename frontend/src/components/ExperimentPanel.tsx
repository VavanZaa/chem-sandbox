import { Beaker, Droplets, Gauge, Play, RotateCcw, Thermometer, Wind } from 'lucide-react'
import { useEffect, useState } from 'react'
import type { Effort, ExperimentDefinition, SubstanceId } from '../lib/types'

type AmountUnit = 'g' | 'mol' | 'mL'
type SubstanceMeta = {
  name: string
  formula: string
  icon: string
  tone: string
  allowedUnits: AmountUnit[]
  molarMass: number
  density?: number
}

const substanceMeta: Record<SubstanceId, SubstanceMeta> = {
  sodium_bicarbonate: { name: 'Baking soda', formula: 'NaHCO₃', icon: '◇', tone: 'mint', allowedUnits: ['mol', 'g'], molarMass: 84.0066 },
  acetic_acid: { name: 'Acetic acid', formula: 'CH₃COOH', icon: '⬡', tone: 'amber', allowedUnits: ['mol', 'g'], molarMass: 60.052 },
  water: { name: 'Water', formula: 'H₂O', icon: '●', tone: 'blue', allowedUnits: ['mol', 'mL', 'g'], molarMass: 18.01528, density: 0.997 },
}

function toMoles(meta: SubstanceMeta, amount: number, unit: AmountUnit): number {
  if (unit === 'mol') return amount
  const mass = unit === 'mL' ? amount * (meta.density ?? 1) : amount
  return mass / meta.molarMass
}

function fromMoles(meta: SubstanceMeta, moles: number, unit: AmountUnit): number {
  if (unit === 'mol') return moles
  const mass = moles * meta.molarMass
  return unit === 'mL' ? mass / (meta.density ?? 1) : mass
}

function readableAmount(value: number): number {
  if (value >= 10) return Number(value.toFixed(2))
  if (value >= 1) return Number(value.toFixed(3))
  return Number(value.toPrecision(4))
}

function equivalentLabel(meta: SubstanceMeta, amount: number, unit: AmountUnit): string {
  const moles = toMoles(meta, amount, unit)
  if (unit !== 'mol') return `≈ ${moles.toPrecision(3)} mol`
  const comparisonUnit: AmountUnit = meta.density ? 'mL' : 'g'
  const comparisonAmount = readableAmount(fromMoles(meta, moles, comparisonUnit))
  return `≈ ${comparisonAmount} ${comparisonUnit}`
}

interface Props {
  experiment: ExperimentDefinition
  setExperiment: (value: ExperimentDefinition) => void
  onRun: () => void
  onReset: () => void
  running: boolean
}

function AmountInput({ value, onCommit, label, disabled }: { value: number; onCommit: (value: number) => void; label: string; disabled: boolean }) {
  const [draft, setDraft] = useState(String(value))
  useEffect(() => setDraft(String(value)), [value])
  return <input aria-label={label} type="number" min="0.000000001" max="10000" step="any" value={draft} disabled={disabled}
    onChange={(event) => setDraft(event.target.value)}
    onBlur={() => { const amount = Number(draft); if (amount > 0 && amount <= 10000 && Number.isFinite(amount)) onCommit(amount); else setDraft(String(value)) }} />
}

export function ExperimentPanel({ experiment, setExperiment, onRun, onReset, running }: Props) {
  const updateAmount = (id: SubstanceId, amount: number) => {
    setExperiment({
      ...experiment,
      substances: experiment.substances.map((item) => (item.id === id ? { ...item, amount } : item)),
    })
  }

  const updateUnit = (id: SubstanceId, unit: AmountUnit) => {
    setExperiment({
      ...experiment,
      substances: experiment.substances.map((item) => {
        if (item.id !== id) return item
        const meta = substanceMeta[id]
        const moles = toMoles(meta, item.amount, item.unit)
        return { ...item, unit, amount: Number(fromMoles(meta, moles, unit).toPrecision(12)) }
      }),
    })
  }

  const updateEnvironment = <K extends keyof ExperimentDefinition['environment']>(
    key: K,
    value: ExperimentDefinition['environment'][K],
  ) => setExperiment({ ...experiment, environment: { ...experiment.environment, [key]: value } })

  const updateEffort = (effort: Effort) => {
    const runs = { low: 4, medium: 8, high: 12 }[effort]
    const optimization_steps = { low: 6, medium: 10, high: 10 }[effort]
    setExperiment({ ...experiment, search: { ...experiment.search, effort, runs, optimization_steps,
      max_depth: { low: 2, medium: 3, high: 4 }[effort],
      max_evaluations_per_run: { low: 12, medium: 24, high: 40 }[effort],
      frontier_width: { low: 4, medium: 6, high: 8 }[effort] } })
  }

  return (
    <aside className="experiment-panel">
      <div className="panel-title-row">
        <div>
          <span className="eyebrow">EXPERIMENT 01</span>
          <h2>Mixture</h2>
        </div>
        <button className="icon-button" onClick={onReset} aria-label="Reset experiment" disabled={running}>
          <RotateCcw size={16} />
        </button>
      </div>

      <label className="control-field wide mode-control"><span>Chemistry workspace</span><select aria-label="Chemistry workspace" value={experiment.mode} disabled={running}
        onChange={(event) => setExperiment({ ...experiment, mode: event.target.value as ExperimentDefinition['mode'], compute: { ...experiment.compute, engine: event.target.value === 'bicarbonate' ? 'all_compatible' : experiment.compute.engine } })}>
        <option value="general">General molecules · SMILES</option><option value="bicarbonate">Bicarbonate demonstration</option></select></label>
      {experiment.mode === 'general' ? <div className="molecule-input-list">
        {experiment.molecules.map((item, index) => <fieldset className="molecule-input" key={index} disabled={running}>
          <label>Name<input aria-label={`Molecule ${index + 1} name`} value={item.name} maxLength={80} onChange={(event) => setExperiment({ ...experiment, molecules: experiment.molecules.map((m, i) => i === index ? { ...m, name: event.target.value } : m) })} /></label>
          <label>SMILES<input aria-label={`Molecule ${index + 1} SMILES`} value={item.smiles} maxLength={500} onChange={(event) => setExperiment({ ...experiment, molecules: experiment.molecules.map((m, i) => i === index ? { ...m, smiles: event.target.value } : m) })} /></label>
          <div className="control-grid"><label>Bulk amount / mol<AmountInput label={`Molecule ${index + 1} amount mol`} value={item.amount_mol} disabled={running} onCommit={(value) => setExperiment({ ...experiment, molecules: experiment.molecules.map((m, i) => i === index ? { ...m, amount_mol: value } : m) })} /></label>
          <label>Cluster copies<input aria-label={`Molecule ${index + 1} copies`} type="number" min="1" max="3" value={item.copies} onChange={(event) => setExperiment({ ...experiment, molecules: experiment.molecules.map((m, i) => i === index ? { ...m, copies: Number(event.target.value) } : m) })} /></label></div>
          <button className="text-button" disabled={experiment.molecules.length === 1} onClick={() => setExperiment({ ...experiment, molecules: experiment.molecules.filter((_, i) => i !== index) })}>Remove molecule</button>
        </fieldset>)}
        <button className="secondary-button" disabled={running || experiment.molecules.length >= 6} onClick={() => setExperiment({ ...experiment, molecules: [...experiment.molecules, { name: 'Water', smiles: 'O', amount_mol: .05, copies: 1 }] })}>+ Add molecule</button>
        <p className="panel-footnote">RDKit generates 3D geometry. Use @/@@ for tetrahedral chirality and / or \ for E/Z double bonds. Unspecified stereochemistry stays unresolved. Up to 96 explicit atoms; newly formed product identities remain exploratory.</p>
      </div> :
      <div className="substance-list">
        {experiment.substances.map((item, index) => {
          const meta = substanceMeta[item.id]
          return (
            <div className={`substance-card ${meta.tone}`} key={item.id}>
              <div className="substance-number">0{index + 1}</div>
              <div className="substance-symbol">{meta.icon}</div>
              <div className="substance-copy">
                <strong>{meta.name}</strong>
                <span>{meta.formula}</span>
                <small>{equivalentLabel(meta, item.amount, item.unit)}</small>
              </div>
              <div className="amount-control">
                <AmountInput label={`${meta.name} amount`} value={item.amount} onCommit={(amount) => updateAmount(item.id, amount)} disabled={running} />
                <select
                  aria-label={`${meta.name} unit`}
                  value={item.unit}
                  onChange={(event) => updateUnit(item.id, event.target.value as AmountUnit)}
                  disabled={running}
                >
                  {meta.allowedUnits.map((unit) => <option key={unit} value={unit}>{unit}</option>)}
                </select>
              </div>
            </div>
          )
        })}
      </div>}

      <div className="section-label">Compute</div><div className="control-grid">
        <label className="control-field wide"><span>Discovery engine</span><select aria-label="Discovery engine" disabled={running} value={experiment.compute.engine} onChange={(event) => setExperiment({ ...experiment, compute: { ...experiment.compute, engine: event.target.value as ExperimentDefinition['compute']['engine'] } })}>
          <option value="all_compatible">All compatible · independent comparison</option><option value="mace_mp">MACE-MP · broad potential</option><option value="ani2x" disabled={experiment.mode !== 'general'}>ANI-2x · neutral HCNOFSCl</option><option value="aimnet2" disabled={experiment.mode !== 'general'}>AIMNet2 · charged main-group</option><option value="aimnet2_rxn" disabled={experiment.mode !== 'general'}>AIMNet2-RXN · neutral HCNO</option></select></label>
        <label className="control-field wide"><span>Device</span><select aria-label="Compute device" disabled={running} value={experiment.compute.device} onChange={(event) => setExperiment({ ...experiment, compute: { ...experiment.compute, device: event.target.value as ExperimentDefinition['compute']['device'] } })}>
          <option value="auto">Auto · prefer GPU</option><option value="cuda">NVIDIA GPU · require CUDA</option><option value="cpu">CPU</option></select></label>
      </div>
      {experiment.compute.engine === 'all_compatible' && <p className="panel-footnote">Each compatible installed model performs its own full bounded search with the same seed and budgets, sequentially on the GPU. Runtime scales with model count. Agreement is reported as discovery consensus, not probability.</p>}

      <details className="search-settings">
        <summary>Search budgets & reproducibility</summary>
        <div className="control-grid">
          {([
            ['runs', 'Independent starts', 1, 50], ['max_depth', 'Max depth', 1, 12],
            ['max_evaluations_per_run', 'Evaluations / start', 1, 160], ['max_candidates_per_state', 'Proposals / state', 1, 24],
            ['frontier_width', 'Width / depth', 1, 24], ['max_states_per_run', 'Admitted states', 1, 80],
            ['optimization_steps', 'Initial relax steps', 1, 60], ['candidate_optimization_steps', 'Candidate relax steps', 1, 60],
            ['seed', 'Random seed', 0, 2147483647],
          ] as const).map(([key, label, min, max]) => <label className="control-field" key={key}><span>{label}</span>
            <input aria-label={label} type="number" min={min} max={max} value={experiment.search[key]} disabled={running}
              onChange={(event) => setExperiment({ ...experiment, search: { ...experiment.search, [key]: Math.min(max, Math.max(min, Math.trunc(Number(event.target.value)))) } })} />
          </label>)}
        </div>
        <button className="secondary-button" disabled={running} onClick={() => setExperiment({ ...experiment, search: { ...experiment.search, seed: crypto.getRandomValues(new Uint32Array(1))[0] % 2147483648 } })}>New random seed</button>
        <p>Same inputs and seed repeat the search. Budgets limit exploration, not reaction time. Two initialization evaluations are additional per start.</p>
      </details>

      <div className="section-label"><Beaker size={14} /> Environment</div>
      <div className="control-grid">
        <label className="control-field">
          <span><Thermometer size={14} /> Temperature</span>
          <div><input type="number" min="250" max="450" value={experiment.environment.temperature_K} onChange={(event) => updateEnvironment('temperature_K', Number(event.target.value))} disabled={running} /><b>K</b></div>
        </label>
        <label className="control-field">
          <span><Gauge size={14} /> Pressure</span>
          <div><input type="number" min="0.1" max="10" step="0.1" value={experiment.environment.pressure_atm} onChange={(event) => updateEnvironment('pressure_atm', Number(event.target.value))} disabled={running} /><b>atm</b></div>
        </label>
        <label className="control-field wide">
          <span><Droplets size={14} /> Solvent</span>
          <select value={experiment.environment.solvent} disabled={running} onChange={(event) => updateEnvironment('solvent', event.target.value as ExperimentDefinition['environment']['solvent'])}><option value="water">Water (context only)</option><option value="vacuum">Vacuum</option><option value="other">Other (context only)</option></select>
        </label>
        <label className="control-field wide">
          <span><Wind size={14} /> Mixing</span>
          <select value={experiment.environment.mixing} onChange={(event) => updateEnvironment('mixing', event.target.value as ExperimentDefinition['environment']['mixing'])} disabled={running}>
            <option value="gentle">Gentle</option>
            <option value="standard">Standard</option>
            <option value="vigorous">Vigorous</option>
          </select>
        </label>
      </div>

      <div className="section-label effort-label"><span>◫</span> Search effort</div>
      <div className="segmented-control">
        {(['low', 'medium', 'high'] as Effort[]).map((effort) => (
          <button key={effort} className={experiment.search.effort === effort ? 'active' : ''} onClick={() => updateEffort(effort)} disabled={running}>
            <span>{effort}</span><small>{{ low: '4', medium: '8', high: '12' }[effort]} runs</small>
          </button>
        ))}
      </div>

      <button className="run-button" onClick={onRun} disabled={running}>
        <span className="run-icon"><Play size={17} fill="currentColor" /></span>
        <span>{running ? 'EXPERIMENT RUNNING' : 'RUN EXPERIMENT'}</span>
        <kbd>↵</kbd>
      </button>
      <p className="panel-footnote">Cluster copies are not bulk concentrations. Temperature and mixing affect search proposals, not physical kinetics; solvent and pressure do not alter neural energies. Use Efficiency for a condition-dependent bulk model.</p>
    </aside>
  )
}
