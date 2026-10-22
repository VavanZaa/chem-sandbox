import { Check, LoaderCircle } from 'lucide-react'
import type { ExperimentJob } from '../lib/types'

const stages = [
  ['preparing', 'Preparing molecular structures'],
  ['searching', 'Searching configuration space'],
  ['detecting', 'Detecting transformations'],
  ['building', 'Building reaction network'],
]

export function RunProgress({ job }: { job: ExperimentJob }) {
  const activeIndex = Math.max(0, stages.findIndex(([key]) => key === job.stage))
  return (
    <div className="progress-shell">
      <div className="progress-visual"><div className="pulse-ring ring-one" /><div className="pulse-ring ring-two" /><div className="progress-core"><LoaderCircle size={26} /></div></div>
      <span className="eyebrow">COMPUTATIONAL EXPERIMENT</span>
      <h2>{job.message}</h2>
      <div className="progress-track"><div style={{ width: `${job.progress}%` }} /></div>
      <span className="progress-percent">{job.progress}%</span>
      <div className="stage-list">
        {stages.map(([key, label], index) => <div key={key} className={index < activeIndex ? 'done' : index === activeIndex ? 'active' : ''}>{index < activeIndex ? <Check size={13} /> : <span />}{label}</div>)}
      </div>
      <p>Each candidate is evaluated with the configured interatomic model. This can take several minutes on CPU.</p>
    </div>
  )
}
