import type { ExperimentDefinition, ExperimentJob, StructureData, KineticsRequest, KineticsResult } from './types'

async function readJson<T>(response: Response): Promise<T> {
  if (!response.ok) {
    const payload = await response.json().catch(() => ({ detail: response.statusText }))
    throw new Error(typeof payload.detail === 'string' ? payload.detail : JSON.stringify(payload.detail))
  }
  return response.json() as Promise<T>
}

export async function requestJson<T>(path: string, payload?: unknown): Promise<T> {
  return readJson<T>(await fetch(`/api${path}`, payload === undefined ? undefined : {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload),
  }))
}

export async function startExperiment(definition: ExperimentDefinition): Promise<ExperimentJob> {
  return readJson<ExperimentJob>(
    await fetch('/api/experiments', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(definition),
    }),
  )
}

export async function getExperiment(jobId: string): Promise<ExperimentJob> {
  return readJson<ExperimentJob>(await fetch(`/api/experiments/${jobId}`))
}

export async function runKinetics(jobId: string, request: KineticsRequest): Promise<KineticsResult> {
  return readJson<KineticsResult>(await fetch(`/api/experiments/${encodeURIComponent(jobId)}/kinetics`, {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(request),
  }))
}

export async function checkHealth(): Promise<boolean> {
  try {
    return (await fetch('/api/health')).ok
  } catch {
    return false
  }
}

export async function getStructure(jobId: string, structureReference: string): Promise<StructureData> {
  const filename = structureReference.split('/').at(-1)
  if (!filename) throw new Error('This state has no saved structure filename')
  return readJson<StructureData>(
    await fetch(`/api/experiments/${encodeURIComponent(jobId)}/structures/${encodeURIComponent(filename)}`),
  )
}
