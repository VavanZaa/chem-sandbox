import { useEffect, useMemo, useRef, useState } from 'react'
import { requestJson } from '../lib/api'
import type { BoltzmannEEDFResult, CandidateRankingJob, CollisionRateResult, CrossSectionValidationResult, PrebioticCapabilities, PrebioticResult, RecordedPrebioticValidation } from '../lib/types'
import { downloadJson } from './AnalysisPanel'

const evidence = (value: number, unit: string) => ({
  value, unit,
  source: 'hypothetical apparatus dry-run value; replace before scientific use',
  method: 'declared scenario parameter', basis: 'hypothetical',
})

const example = {
  id: 'miller_urey_coupled_plasma_gate',
  species: [
    { id: 'H2O', name: 'Water', composition: { H: 2, O: 1 }, charge: 0, smiles: 'O',
      phase_behavior: 'condensable', evaporation_rate: evidence(0.02, '1/s'), condensation_rate: evidence(0.08, '1/s') },
    { id: 'CH4', name: 'Methane', composition: { C: 1, H: 4 }, charge: 0, smiles: 'C', phase_behavior: 'permanent_gas' },
    { id: 'NH3', name: 'Ammonia', composition: { N: 1, H: 3 }, charge: 0, smiles: 'N', phase_behavior: 'volatile',
      condensation_rate: evidence(0.01, '1/s') },
    { id: 'H2', name: 'Hydrogen', composition: { H: 2 }, charge: 0, smiles: '[H][H]', phase_behavior: 'permanent_gas' },
  ],
  reactions: [],
  initial_inventory: [
    { species_id: 'H2O', compartment: 'boiler_liquid', amount: { value: 2, unit: 'mol' } },
    { species_id: 'CH4', compartment: 'gas_loop', amount: { value: 0.08, unit: 'mol' } },
    { species_id: 'NH3', compartment: 'gas_loop', amount: { value: 0.04, unit: 'mol' } },
    { species_id: 'H2', compartment: 'gas_loop', amount: { value: 0.05, unit: 'mol' } },
  ],
  apparatus: {
    cycles: 20,
    boiler_liquid_volume: { value: 1, unit: 'L' }, gas_loop_volume: { value: 5, unit: 'L' },
    discharge_zone_volume: {
      value: 0.05, unit: 'L',
      source: 'illustrative 50 mL active spark zone; replace with measured apparatus geometry',
      method: 'declared active-volume scenario', basis: 'hypothetical',
    },
    condensate_trap_volume: { value: 0.25, unit: 'L' }, gas_circulation_rate: { value: 0.2, unit: 'L/s' },
    pressure: { value: 1, unit: 'atm' }, boiler_temperature: { value: 373.15, unit: 'K' },
    discharge_temperature: { value: 400, unit: 'K' }, condenser_temperature: { value: 280, unit: 'K' },
    boil_duration: { value: 30, unit: 's' }, discharge_duration: { value: 20, unit: 's' },
    condense_duration: { value: 30, unit: 's' }, condensate_return_fraction: 1,
  },
  spark: {
    mode: 'boltzmann_eedf', voltage: { value: 30000, unit: 'V' }, current: { value: 0.002, unit: 'A' },
    pulse_width: { value: 0.001, unit: 's' }, frequency: { value: 100, unit: 'Hz' }, coupling_efficiency: 0.3,
    boltzmann: {
      solver: 'bolos_two_term', collision_set: 'lxcat_morgan_2026-09-22', reduced_electric_field_Td: 200,
      electron_number_density: {
        value: 1e15, unit: '1/m^3', source: 'illustrative plasma density; replace with a measured or solved value',
        method: 'declared zero-dimensional discharge condition', basis: 'hypothetical',
      },
      energy_max_eV: 60, energy_cells: 400, initial_electron_temperature_eV: 2,
      maximum_iterations: 200, relative_tolerance: 1e-6, acknowledge_two_term_limit: true,
    },
    source: 'hypothetical coupled-plasma program; not the historical apparatus', basis: 'hypothetical',
  },
  interventions: [{
    cycle: 10, additions: [{ species_id: 'H2', compartment: 'gas_loop', amount: { value: 0.01, unit: 'mol' } }],
    note: 'Demonstration: add gas after nine complete circulation cycles',
  }],
  mechanism_origin: 'user_declared', knowledge_policy: 'blind_physics_only', acknowledge_reduced_order_limit: true,
}

const historicalMaterialExample = {
  ...example,
  id: 'miller_urey_1953_material_time_profile',
  initial_inventory: [
    { species_id: 'H2O', compartment: 'boiler_liquid', amount: { value: 11.101687012358397, unit: 'mol' } },
    { species_id: 'CH4', compartment: 'gas_loop', amount: { value: 0.053781637532012414, unit: 'mol' } },
    { species_id: 'NH3', compartment: 'gas_loop', amount: { value: 0.053781637532012414, unit: 'mol' } },
    { species_id: 'H2', compartment: 'gas_loop', amount: { value: 0.026890818766006207, unit: 'mol' } },
  ],
  apparatus: {
    ...example.apparatus,
    cycles: 168,
    boiler_liquid_volume: { value: 0.5, unit: 'L' },
    pressure: { value: 0.6894736842105263, unit: 'atm' },
    boil_duration: { value: 1200, unit: 's' },
    discharge_duration: { value: 1200, unit: 's' },
    condense_duration: { value: 1200, unit: 's' },
  },
  spark: {
    ...example.spark,
    voltage: { value: 60000, unit: 'V' },
    current: { value: 0.001, unit: 'A' },
    source: '1953 material and one-week time profile; discharge current, pulse structure, coupling, and reduced field remain hypothetical sensitivity values',
  },
  interventions: [],
}

const n2ReactiveExample = {
  ...example,
  id: 'n2_reactive_pathway_probe',
  species: [
    ...example.species.filter((species) => species.id !== 'NH3'),
    { id: 'N2', name: 'Nitrogen', composition: { N: 2 }, charge: 0, smiles: 'N#N', phase_behavior: 'permanent_gas' },
  ],
  initial_inventory: [
    { species_id: 'H2O', compartment: 'boiler_liquid', amount: { value: 2, unit: 'mol' } },
    { species_id: 'CH4', compartment: 'gas_loop', amount: { value: 0.08, unit: 'mol' } },
    { species_id: 'N2', compartment: 'gas_loop', amount: { value: 0.08, unit: 'mol' } },
    { species_id: 'H2', compartment: 'gas_loop', amount: { value: 0.08, unit: 'mol' } },
  ],
  apparatus: {
    ...example.apparatus,
    cycles: 20,
    operation_mode: 'continuous_overlap',
    integration_substeps_per_cycle: 2,
  },
  spark: {
    ...example.spark,
    boltzmann: {
      ...example.spark.boltzmann,
      collision_set: 'hybrid_morgan_phelps_tabata_2026-09-23',
    },
    source: 'reviewed N2 dissociation pathway probe; electrical controls remain hypothetical',
  },
  interventions: [],
}

const collisionExample = {
  id: 'abstract_ionization_channel', target_species: 'X', product_species: 'X+', kind: 'ionization',
  energy_levels_eV: [0, 5, 10, 15, 20], cross_sections_m2: [0, 0, 1e-22, 4e-21, 2e-21],
  source: 'abstract UI validation fixture; replace with a citable collision dataset',
  method: 'tabulated electron collision cross section', basis: 'experimental', mechanism_origin: 'user_declared',
}

const collisionRateExample = {
  dataset: collisionExample,
  eedf: {
    mode: 'maxwellian', mean_energy: { value: 5, unit: 'eV' },
    source: 'illustrative Maxwellian assumption; replace with a solved or measured EEDF',
    method: 'normalized Maxwell-Boltzmann electron-energy distribution', basis: 'hypothetical',
  },
  electron_number_density: {
    value: 1e15, unit: '1/m^3', source: 'illustrative density; replace with a measured or solved value',
    method: 'declared UI scenario parameter', basis: 'hypothetical',
  },
}

const boltzmannExample = {
  solver: 'bolos_two_term', collision_set: 'hybrid_morgan_phelps_tabata_2026-09-23',
  mole_fractions: { H2: 0.35, H2O: 0.25, CH4: 0.2, NH3: 0.1, N2: 0.1 },
  gas_temperature_K: 400, reduced_electric_field_Td: 200,
  energy_max_eV: 60, energy_cells: 400, initial_electron_temperature_eV: 2,
  maximum_iterations: 200, relative_tolerance: 1e-6, acknowledge_two_term_limit: true,
}

const communityCh4CollisionSet = 'hybrid_morgan_community_phelps_tabata_2026-09-23'
const istH2CommunityCh4CollisionSet = 'hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23'
const istCo2O2CollisionSet = 'hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26'
const carbonOxygenEedfExample = {
  ...boltzmannExample,
  collision_set: istCo2O2CollisionSet,
  mole_fractions: { H2: 0.01, CO2: 0.05, O2: 0.001, H2O: 0.02, N2: 0.919 },
  gas_temperature_K: 298.15,
}

type HCNAtomCeiling = {
  cyanide_family_total_ceiling_mol: number
  molecular_hcn_total_ceiling_mol: number
  maximum_net_new_cyanide_family_mol: number
  available_atom_mol: Record<string, number>
  interpretation: string
}

type ReducedHCNResult = {
  final_gas_hcn_mol: number
  final_aqueous_hcn_plus_cn_mol: number
  net_hcn_plus_cn_formed_mol: number
  source_extents_mol: Record<string, number>
  atom_ceiling: HCNAtomCeiling
  atom_ledger_closed: boolean
  experimental_yield_validated: boolean
  hypothetical_input_fields: string[]
  limitations: string[]
}

type ReducedHCNAudit = {
  numerical_checks_passed: boolean
  audit_outcome: 'checks_passed' | 'check_failed' | 'incomplete'
  solver_step_consistent: boolean
  coarse_maximum_solver_step_s: number
  refined_maximum_solver_step_s: number
  source_off_control: { final_total_hcn_plus_cn_mol: number; below_reporting_floor: boolean }
  atom_removal_controls: Record<string, { status: string; final_total_hcn_plus_cn_mol: number | null }>
  step_comparisons: Record<string, { within_declared_numerical_tolerance: boolean }>
  refined_result: ReducedHCNResult
  interpretation: string
}

// The one source hazard is deliberately zero. No N2 fragmentation is supplied:
// the reviewed total dissociation record does not identify ground-state atoms.
const reducedHCNTemplate = {
  initial_gas_mol: { CH4: 0.01, N2: 0.02, H2: 0.02, H2O: 0.01 },
  gas_volume: evidence(1, 'L'), water_volume: evidence(0.05, 'L'),
  temperature: evidence(298.15, 'K'), duration: evidence(60, 's'),
  fixed_pH: { value: 7, source: 'illustrative pH, replace with independent measurement',
    method: 'declared fixed-pH scenario', basis: 'hypothetical' },
  hcn_henry_hcp: { value: 0.12, unit: 'mol/(m^3*Pa)',
    source: 'NIST HCN water solubility review at 298 K: 12 mol/(kg bar)',
    method: 'converted with 1000 kg/m^3 water density', basis: 'literature' },
  hcn_liquid_side_kla: evidence(0.001, '1/s'),
  sources: [
    { id: 'methane_cleavage', parent: 'CH4', products: { CH3: 1, H: 1 },
      parent_conversion_hazard: evidence(0, '1/s'),
      product_state_source: 'ground-state fragment identity requires independent evidence',
      product_state_basis: 'hypothetical' },
  ],
  sample_count: 24, hcn_assay_used_to_set_inputs: false,
}

function nonzero(values: Record<string, number>) {
  const entries = Object.entries(values).filter(([, value]) => Math.abs(value) > 1e-12)
    .sort((left, right) => Math.abs(right[1]) - Math.abs(left[1]))
  const shown = entries.slice(0, 12).map(([key, value]) => `${key} ${value.toPrecision(3)} mol`).join(' · ')
  return shown ? `${shown}${entries.length > 12 ? ` · +${entries.length - 12} more` : ''}` : 'none'
}

export function OriginsPanel() {
  const resultRef = useRef<HTMLElement | null>(null)
  const [definition, setDefinition] = useState(JSON.stringify(example, null, 2))
  const [atomCeiling, setAtomCeiling] = useState<HCNAtomCeiling | null>(null)
  const [reducedDefinition, setReducedDefinition] = useState(JSON.stringify(reducedHCNTemplate, null, 2))
  const [reducedOutput, setReducedOutput] = useState<ReducedHCNResult | null>(null)
  const [reducedAudit, setReducedAudit] = useState<ReducedHCNAudit | null>(null)
  const [reducedError, setReducedError] = useState<string | null>(null)
  const [reducedBusy, setReducedBusy] = useState(false)
  const [output, setOutput] = useState<PrebioticResult | null>(null)
  const [capabilities, setCapabilities] = useState<PrebioticCapabilities | null>(null)
  const [recordedValidation, setRecordedValidation] = useState<RecordedPrebioticValidation | null>(null)
  const [collisionDefinition, setCollisionDefinition] = useState(JSON.stringify(collisionExample, null, 2))
  const [collisionResult, setCollisionResult] = useState<CrossSectionValidationResult | null>(null)
  const [rateDefinition, setRateDefinition] = useState(JSON.stringify(collisionRateExample, null, 2))
  const [rateResult, setRateResult] = useState<CollisionRateResult | null>(null)
  const [boltzmannDefinition, setBoltzmannDefinition] = useState(JSON.stringify(boltzmannExample, null, 2))
  const [boltzmannResult, setBoltzmannResult] = useState<BoltzmannEEDFResult | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)
  const [rankingJob, setRankingJob] = useState<CandidateRankingJob | null>(null)
  const [rankingError, setRankingError] = useState<string | null>(null)
  const [evaluationBusy, setEvaluationBusy] = useState(false)
  useEffect(() => { setAtomCeiling(null) }, [definition])
  useEffect(() => { requestJson<PrebioticCapabilities>('/prebiotic/capabilities').then(setCapabilities).catch(() => {}) }, [])
  useEffect(() => { requestJson<RecordedPrebioticValidation>('/prebiotic/validation/recorded').then(setRecordedValidation).catch(() => {}) }, [])
  useEffect(() => {
    if (!output) return
    window.requestAnimationFrame(() => resultRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' }))
  }, [output])
  const shownCycles = useMemo(() => output ? output.cycles.filter((cycle, index) =>
    index === 0 || index === output.cycles.length - 1 || cycle.intervention_notes.length > 0) : [], [output])
  const plasmaCycles = output?.plasma_model?.cycle_summaries ?? []
  const lastPlasma = plasmaCycles.length ? plasmaCycles[plasmaCycles.length - 1] : null
  const firstResolvedO1d = plasmaCycles.find((row) =>
    row.o1d_fate_diagnostic?.status === 'conditional_fate_resolved')?.o1d_fate_diagnostic
  const lastO1d = lastPlasma?.o1d_fate_diagnostic
  const heavyCycles = output?.plasma_model?.heavy_species_cycle_summaries ?? []
  const firstHeavy = heavyCycles.length ? heavyCycles[0] : null
  const firstPathwayAudit = firstHeavy?.integrated_pathway_diagnostics
  const lastHeavy = heavyCycles.length ? heavyCycles[heavyCycles.length - 1] : null
  const partitionCycles = output?.plasma_model?.phase_partition_cycle_summaries ?? []
  const lastPartition = partitionCycles.length ? partitionCycles[partitionCycles.length - 1] : null
  const reactiveSources = lastPlasma?.source_terms.filter((source) => source.reactive).slice(0, 8) ?? []
  const lastN2Dissociation = lastPlasma?.source_terms.find((source) =>
    source.target === 'N2' && source.product_label === 'N + N')
  const n2Dissociation = boltzmannResult?.rate_coefficients.find((rate) =>
    rate.target === 'N2' && rate.kind === 'dissociation' && rate.product_label === 'N + N')
  const emergedGasSpecies = lastHeavy ? Object.entries(lastHeavy.formed_mol)
    .filter(([species]) => !['H2', 'H2O', 'CH4', 'NH3', 'N2'].includes(species))
    .sort((left, right) => right[1] - left[1]).slice(0, 10) : []
  const shadowComparison = lastHeavy?.mechanism_comparison
  const chargedShadow = shadowComparison?.charged_shadow
  const propagatedCharged = output?.plasma_model.propagated_charged_branch
  const chargedHcnPathways = chargedShadow?.final_instantaneous_pathway_diagnostics?.HCN
  const literatureShadow = shadowComparison?.literature_hcn_shadow
  const shadowSpecies = shadowComparison?.species_comparison
    ? ['HCN', 'H2CN', 'CH2NH', 'NH3', 'CH2O'].filter((species) => shadowComparison.species_comparison?.[species])
    : []
  const capturedPrecursors = output ? ['HCN', 'CH2O', 'CH3CHO', 'CH3OH', 'NH3']
    .map((species) => [species, (output.final_compartments_mol.boiler_liquid?.[species] ?? 0)
      + (output.final_compartments_mol.condensate_trap?.[species] ?? 0)
      + (output.final_compartments_mol.collected_condensate?.[species] ?? 0)] as const)
    .filter(([, amount]) => amount > 1e-14) : []
  const chargedCapturedPrecursors = propagatedCharged ? ['HCN', 'CH2O', 'CH3CHO', 'CH3OH', 'NH3']
    .map((species) => [species, (propagatedCharged.final_compartments_mol.boiler_liquid?.[species] ?? 0)
      + (propagatedCharged.final_compartments_mol.condensate_trap?.[species] ?? 0)
      + (propagatedCharged.final_compartments_mol.collected_condensate?.[species] ?? 0)] as const)
    .filter(([, amount]) => amount > 1e-14) : []
  const aminoAcidMotifCandidates = output?.aqueous_discovery.formula_records
    .flatMap((record) => record.candidates
      .filter((candidate) => candidate.motifs.includes('alpha_amino_carboxylic_acid'))
      .map((candidate) => ({ ...candidate, sources: record.source_combinations }))) ?? []

  async function loadStriblingMillerProfile(coupledScenario = false, communityMethane = false, istHydrogen = false, quickRun = false) {
    setBusy(true); setError(null); setOutput(null)
    try {
      const profile = await requestJson<Record<string, unknown>>('/prebiotic/benchmarks/stribling-miller-1986/profile')
      if (coupledScenario) {
        const apparatus = profile.apparatus as Record<string, unknown>
        apparatus.primary_source_forcing = 'time_averaged'
        apparatus.operation_mode = 'continuous_overlap'
        apparatus.integration_substeps_per_cycle = 4
        apparatus.condensate_nonwater_return_fraction = 0
        apparatus.condenser_phase_transfer_model = 'coupled_flow_transport'
        apparatus.condenser_contact_water_holdup = {
          value: 0.1, unit: 'L', basis: 'hypothetical',
          source: 'illustrative contact-water hold-up, not measured historical geometry',
          method: 'declared coupled-flow sensitivity scenario',
        }
      }
      if (communityMethane || istHydrogen) {
        const spark = profile.spark as Record<string, unknown>
        const boltzmann = spark.boltzmann as Record<string, unknown>
        boltzmann.collision_set = istHydrogen ? istH2CommunityCh4CollisionSet : communityCh4CollisionSet
      }
      if (quickRun) {
        const apparatus = profile.apparatus as Record<string, unknown>
        apparatus.cycles = 1
        apparatus.integration_substeps_per_cycle = 2
      }
      setDefinition(JSON.stringify(profile, null, 2))
      if (quickRun) setOutput(await requestJson<PrebioticResult>('/prebiotic/apparatus/simulate', profile))
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : quickRun ? 'Quick HCN scenario failed' : 'Measurement-timing profile failed to load')
    } finally {
      setBusy(false)
    }
  }

  async function run() {
    setBusy(true); setError(null); setOutput(null); setRankingJob(null); setRankingError(null)
    try { setOutput(await requestJson<PrebioticResult>('/prebiotic/apparatus/simulate', JSON.parse(definition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Apparatus simulation failed') }
    finally { setBusy(false) }
  }

  async function calculateAtomCeiling() {
    setError(null); setAtomCeiling(null)
    try { setAtomCeiling(await requestJson<HCNAtomCeiling>('/prebiotic/hcn/atom-ceiling', JSON.parse(definition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Could not calculate atom ceiling') }
  }

  async function runReducedHCN() {
    setReducedBusy(true); setReducedError(null); setReducedOutput(null); setReducedAudit(null)
    try { setReducedOutput(await requestJson<ReducedHCNResult>('/prebiotic/hcn/reduced-mechanistic', JSON.parse(reducedDefinition))) }
    catch (caught) { setReducedError(caught instanceof Error ? caught.message : 'Reduced HCN calculation failed') }
    finally { setReducedBusy(false) }
  }

  async function auditReducedHCN() {
    setReducedBusy(true); setReducedError(null); setReducedOutput(null); setReducedAudit(null)
    try {
      const audit = await requestJson<ReducedHCNAudit>('/prebiotic/hcn/reduced-mechanistic/audit', JSON.parse(reducedDefinition))
      setReducedAudit(audit); setReducedOutput(audit.refined_result)
    } catch (caught) { setReducedError(caught instanceof Error ? caught.message : 'Reduced HCN audit failed') }
    finally { setReducedBusy(false) }
  }

  async function rankFrozenCandidates() {
    if (!output) return
    setRankingError(null)
    try {
      let job = await requestJson<CandidateRankingJob>('/prebiotic/post-discovery/rank', {
        discovery: output.aqueous_discovery,
        discovery_artifact_sha256: output.post_discovery_validation.discovery_artifact_sha256,
        engines: ['mace_mp', 'ani2x', 'aimnet2', 'aimnet2_rxn'],
        device: 'auto', optimization_steps: 60, force_threshold_eV_A: 0.15,
        acknowledge_vacuum_energy_limit: true,
      })
      setRankingJob(job)
      while (job.status === 'queued' || job.status === 'running') {
        await new Promise((resolve) => window.setTimeout(resolve, 1500))
        job = await requestJson<CandidateRankingJob>(`/prebiotic/post-discovery/rank/${encodeURIComponent(job.id)}`)
        setRankingJob(job)
      }
      if (job.status === 'failed') setRankingError(job.error ?? 'Candidate ranking failed')
    } catch (caught) {
      setRankingError(caught instanceof Error ? caught.message : 'Candidate ranking failed')
    }
  }

  async function revealHeldOutRanking() {
    if (!rankingJob?.result) return
    setEvaluationBusy(true); setRankingError(null)
    try {
      const evaluation = await requestJson<NonNullable<CandidateRankingJob['evaluation']>>(
        `/prebiotic/post-discovery/rank/${encodeURIComponent(rankingJob.id)}/held-out-evaluation`,
        {},
      )
      setRankingJob({ ...rankingJob, evaluation })
    } catch (caught) {
      setRankingError(caught instanceof Error ? caught.message : 'Held-out evaluation failed')
    } finally {
      setEvaluationBusy(false)
    }
  }

  async function validateCollision() {
    setBusy(true); setError(null); setCollisionResult(null)
    try { setCollisionResult(await requestJson<CrossSectionValidationResult>('/prebiotic/plasma/cross-sections/validate', JSON.parse(collisionDefinition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Cross-section validation failed') }
    finally { setBusy(false) }
  }

  async function calculateCollisionRate() {
    setBusy(true); setError(null); setRateResult(null)
    try { setRateResult(await requestJson<CollisionRateResult>('/prebiotic/plasma/rate-coefficient', JSON.parse(rateDefinition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Electron-collision rate calculation failed') }
    finally { setBusy(false) }
  }

  async function solveBoltzmannEEDF() {
    setBusy(true); setError(null); setBoltzmannResult(null)
    try { setBoltzmannResult(await requestJson<BoltzmannEEDFResult>('/prebiotic/plasma/eedf/solve', JSON.parse(boltzmannDefinition))) }
    catch (caught) { setError(caught instanceof Error ? caught.message : 'Boltzmann EEDF solve failed') }
    finally { setBusy(false) }
  }

  return <div className="analysis-panel protocol-panel">
    <span className="eyebrow">PREBIOTIC APPARATUS · BLIND-PRODUCT FIREWALL</span>
    <h2>Closed spark-discharge and condensation loop.</h2>
      <div className="analysis-warning">The gas inventory drives a self-consistent plasma calculation in every discharge cycle. Resolved, energy-limited neutral fragments feed a Cantera gas network, and reviewed Henry-law data move supported products into the water phase. Broad structure enumeration, post-discovery cycle replay, a frozen held-out rank evaluation, and a measured precursor energy-yield gate are reported separately. Transition paths, full aqueous kinetics, and validated yields remain open.</div>
    {capabilities && <div className="summary-cards">
      <div><span>Apparatus</span><strong>{capabilities.apparatus.available ? 'Ready' : 'Unavailable'}</strong><small>{capabilities.apparatus.compartments.length} compartments</small></div>
      <div><span>Cantera plasma</span><strong>{capabilities.cantera_cross_section_plasma.available ? 'Verified' : 'Unavailable'}</strong><small>v{capabilities.cantera_cross_section_plasma.version}</small></div>
      <div><span>Boltzmann EEDF</span><strong>{capabilities.boltzmann_eedf.available ? 'Ready' : 'Unavailable'}</strong><small>{capabilities.boltzmann_eedf.targets.join(' / ')}</small></div>
      <div><span>Pavan methane-plasma shadow</span><strong>{capabilities.cantera_cross_section_plasma.pavan_charge_balanced_shadow_connected ? 'Charge-balanced' : capabilities.cantera_cross_section_plasma.pavan_correlated_shadow_connected ? 'Neutral only' : 'Unavailable'}</strong><small>{capabilities.cantera_cross_section_plasma.pavan_methane_plasma_reactions ?? 0} packaged reactions · correlated</small></div>
      <div><span>Literature HCN network</span><strong>{capabilities.cantera_cross_section_plasma.kamarinopoulou_hcn_shadow_connected ? 'Connected' : 'Unavailable'}</strong><small>{capabilities.cantera_cross_section_plasma.kamarinopoulou_hcn_neutral_reactions ?? 0} neutral reactions · 62 fixed electron rates excluded</small></div>
      <div><span>Known products in discovery</span><strong>{capabilities.blind_validation.known_products_available_to_discovery ? 'Yes' : 'Blocked'}</strong></div>
    </div>}
    {recordedValidation && <section aria-label="Recorded Miller-Urey scenario result">
      <div className="analysis-title"><h3>Archived full-week development scenarios</h3></div>
      <p className="analysis-warning">These saved 168-hour runs predate the correction to the N₂ electron-impact fragment states: they assigned both fragments to ground-state N atoms. They are retained as a development record, not current-code predictions or validated yields. Their electrical and transport inputs were hypothetical.</p>
      <div className="analysis-table-wrap"><table><thead><tr><th>Source / condenser</th><th>Non-water return</th><th>Cycles / steps</th><th>Primary HCN</th><th>Primary formaldehyde</th><th>Charged HCN</th><th>Charged formaldehyde</th></tr></thead><tbody>
        {recordedValidation.runs.map((run) => <tr key={run.artifact_id}><td>{(run.primary_source_forcing ?? 'impulse').replaceAll('_', ' ')}<small> · {(run.condenser_phase_transfer_model ?? 'step_equilibrium').replaceAll('_', ' ')}{run.condenser_contact_water_holdup_L ? ` · ${run.condenser_contact_water_holdup_L} L hold-up` : ''} · {run.model_revision ? run.model_revision.replaceAll('_', ' ') : 'legacy density treatment'}</small></td><td>{((run.condensate_nonwater_return_fraction ?? 1) * 100).toFixed(0)}%</td><td>{run.cycles_completed} / {run.integration_steps}</td><td>{run.primary.hcn_zero_time_intercept_nmol_per_J.toExponential(4)}<small> · R² {run.primary.hcn_fit_r_squared.toFixed(3)}</small></td><td>{run.primary.formaldehyde_zero_time_intercept_nmol_per_J.toExponential(4)}<small> · R² {run.primary.formaldehyde_fit_r_squared.toFixed(3)}</small></td><td>{run.propagated_charged.hcn_zero_time_intercept_nmol_per_J.toExponential(4)}<small> · R² {run.propagated_charged.hcn_fit_r_squared.toFixed(3)}</small></td><td>{run.propagated_charged.formaldehyde_zero_time_intercept_nmol_per_J.toExponential(4)}<small> · R² {run.propagated_charged.formaldehyde_fit_r_squared.toFixed(3)}</small></td></tr>)}
      </tbody></table></div>
      <p className="analysis-warning">All four precursor columns are in nmol/J and use a recorded zero-time linear fit; R² near zero means that fit poorly describes the time series. Different source and trap treatments show model sensitivity, not independent votes or product probabilities.</p>
      {recordedValidation.convergence_studies.map((study) => <div key={study.artifact_id}>
        <div className="analysis-title"><h3>24-hour timestep study · {(study.primary_source_forcing ?? 'impulse').replaceAll('_', ' ')} · {(study.condenser_phase_transfer_model ?? 'step_equilibrium').replaceAll('_', ' ')} · {((study.condensate_nonwater_return_fraction ?? 1) * 100).toFixed(0)}% non-water return · {study.model_revision ? study.model_revision.replaceAll('_', ' ') : 'legacy density treatment'}</h3></div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Substeps/hour</th><th>Primary HCN</th><th>Primary formaldehyde</th><th>Charged HCN</th><th>Charged formaldehyde</th></tr></thead><tbody>
          {study.results.map((row) => <tr key={row.substeps}><td>{row.substeps}</td><td>{row.primary_water_hcn_mol.toExponential(4)}</td><td>{row.primary_water_formaldehyde_mol.toExponential(4)}</td><td>{row.charged_water_hcn_mol.toExponential(4)}</td><td>{row.charged_water_formaldehyde_mol.toExponential(4)}</td></tr>)}
        </tbody></table></div>
        <p className="analysis-warning">Water-phase amount in mol after 24 hours · {study.convergence_status.replaceAll('_', ' ')}</p>
      </div>)}
      {recordedValidation.sensitivity_studies?.map((study) => <div key={study.artifact_id}>
        <div className="analysis-title"><h3>24-hour contact-water sensitivity · coupled flow</h3></div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Contact hold-up / L</th><th>Primary HCN</th><th>Primary formaldehyde</th><th>Charged HCN</th><th>Charged formaldehyde</th></tr></thead><tbody>
          {study.results.map((row) => <tr key={row.condenser_contact_water_holdup_L}><td>{row.condenser_contact_water_holdup_L}</td><td>{row.primary_water_hcn_mol.toExponential(4)}</td><td>{row.primary_water_formaldehyde_mol.toExponential(4)}</td><td>{row.charged_water_hcn_mol.toExponential(4)}</td><td>{row.charged_water_formaldehyde_mol.toExponential(4)}</td></tr>)}
        </tbody></table></div>
        <p className="analysis-warning">Amounts in mol after 24 hours. Hold-up values are hypothetical sensitivity cases, not measurements or a statistical confidence interval. No benchmark yields were used to select them.</p>
      </div>)}
    </section>}
    <section aria-label="Reduced HCN study">
      <div className="analysis-title"><h3>Reduced HCN study · independent source inputs</h3></div>
      <p className="analysis-warning">This route skips the electron-distribution and spark-geometry calculation. Supply independently established fragment-production rates, HCN gas–water contact, and pH. The template's methane rate is zero and it supplies no nitrogen-fragment channel: the total N₂ collision data do not justify N₂ → 2 ground-state N. Never set source rates using an HCN assay intended for validation.</p>
      <details><summary>View or edit source-driven study inputs</summary>
        <textarea className="network-editor" aria-label="Reduced HCN study JSON" spellCheck={false} value={reducedDefinition}
          onChange={(event) => { setReducedDefinition(event.target.value); setReducedOutput(null); setReducedAudit(null); setReducedError(null) }} />
      </details>
      <div className="button-row"><button className="secondary-button" disabled={reducedBusy} onClick={runReducedHCN}>
        {reducedBusy ? 'Solving conditional study…' : 'Run reduced HCN study'}
      </button><button className="secondary-button" disabled={reducedBusy} onClick={auditReducedHCN}>
        Audit solver and controls
      </button></div>
      {reducedError && <p className="analysis-warning" role="alert">{reducedError}</p>}
      {reducedAudit && <div aria-live="polite">
        <p className="analysis-warning">Numerical audit: {reducedAudit.audit_outcome.replaceAll('_', ' ')} at {reducedAudit.coarse_maximum_solver_step_s.toPrecision(3)} s and {reducedAudit.refined_maximum_solver_step_s.toPrecision(3)} s maximum steps. Source-off background: {reducedAudit.source_off_control.final_total_hcn_plus_cn_mol.toExponential(3)} mol. Atom-removal controls: {Object.entries(reducedAudit.atom_removal_controls).map(([key, control]) => `${key.replaceAll('_', ' ')} ${control.status.replaceAll('_', ' ')}`).join(' · ')}. {reducedAudit.interpretation}</p>
        {!reducedAudit.solver_step_consistent && <p className="analysis-warning" role="alert">Step-sensitive outputs: {Object.entries(reducedAudit.step_comparisons).filter(([, check]) => !check.within_declared_numerical_tolerance).map(([name]) => name).join(', ')}</p>}
      </div>}
      {reducedOutput && <div aria-live="polite">
        <div className="summary-cards">
          <div><span>Gas HCN</span><strong>{reducedOutput.final_gas_hcn_mol.toExponential(3)}</strong><small>mol · conditional model</small></div>
          <div><span>Water HCN + CN⁻</span><strong>{reducedOutput.final_aqueous_hcn_plus_cn_mol.toExponential(3)}</strong><small>mol · fixed-pH HCN equivalents</small></div>
          <div><span>Atom ceiling</span><strong>{reducedOutput.atom_ceiling.cyanide_family_total_ceiling_mol.toExponential(3)}</strong><small>mol maximum from C and N, not an achievable yield</small></div>
          <div><span>Atom ledger</span><strong>{reducedOutput.atom_ledger_closed ? 'Closed' : 'Failed'}</strong><small>full heat and aqueous charge ledgers remain open</small></div>
          <div><span>Hypothetical inputs</span><strong>{reducedOutput.hypothetical_input_fields.length}</strong><small>not experimental measurements</small></div>
        </div>
        <p className="analysis-warning">Source conversion: {Object.entries(reducedOutput.source_extents_mol).map(([id, amount]) => `${id} ${amount.toExponential(3)} mol`).join(' · ') || 'none'}. This is not experimentally validated, and the GRI network already contains HCN reactions.</p>
        {reducedOutput.limitations.map((item) => <p className="analysis-warning" key={item}>{item}</p>)}
      </div>}
    </section>
    <textarea className="network-editor" aria-label="Prebiotic apparatus JSON" spellCheck={false} value={definition}
      onChange={(event) => { setDefinition(event.target.value); setOutput(null); setAtomCeiling(null); setError(null) }} />
    <div className="button-row"><button className="primary-action" disabled={busy} onClick={run}>{busy ? 'Solving cycles…' : 'Run coupled apparatus'}</button>
      <button className="secondary-button" disabled={busy} onClick={calculateAtomCeiling}>HCN atom ceiling</button>
      <button className="secondary-button" disabled={busy} onClick={() => { setDefinition(JSON.stringify(example, null, 2)); setOutput(null) }}>Reset coupled example</button>
      <button className="secondary-button" disabled={busy} onClick={() => { setDefinition(JSON.stringify(historicalMaterialExample, null, 2)); setOutput(null); setError(null) }}>Load 1953 material/time profile</button>
      <button className="secondary-button" disabled={busy} onClick={() => { setDefinition(JSON.stringify(n2ReactiveExample, null, 2)); setOutput(null); setError(null) }}>Load N2 reactive pathway</button>
      <button className="secondary-button" disabled={busy} onClick={() => loadStriblingMillerProfile()}>Load 1 h–7 d measurement profile</button>
      <button className="secondary-button" disabled={busy} onClick={() => loadStriblingMillerProfile(true)}>Load coupled-flow sensitivity</button>
      <button className="secondary-button" disabled={busy} onClick={() => loadStriblingMillerProfile(true, true)}>Load reviewed methane-set sensitivity</button>
      <button className="secondary-button" disabled={busy || !capabilities?.boltzmann_eedf.collision_sets?.some((item) => item.dataset_id === istH2CommunityCh4CollisionSet && item.available)} onClick={() => loadStriblingMillerProfile(true, true, true)}>Load opt-in IST H₂ + methane sensitivity</button>
      <button className="secondary-button" disabled={busy || !capabilities?.boltzmann_eedf.collision_sets?.some((item) => item.dataset_id === istH2CommunityCh4CollisionSet && item.available)} onClick={() => loadStriblingMillerProfile(true, true, true, true)}>Run quick HCN scenario</button></div>
    {atomCeiling && <div className="summary-cards" aria-live="polite">
      <div><span>HCN + CN⁻ atom ceiling</span><strong>{atomCeiling.cyanide_family_total_ceiling_mol.toExponential(3)}</strong><small>mol · min(C, N) across feed and scheduled additions</small></div>
      <div><span>Molecular HCN atom ceiling</span><strong>{atomCeiling.molecular_hcn_total_ceiling_mol.toExponential(3)}</strong><small>mol · min(C, N, H)</small></div>
      <div><span>Predicted yield</span><strong>No</strong><small>Atoms alone cannot predict conversion or collection</small></div>
    </div>}
    <p className="analysis-warning">The methane-set sensitivity uses a private, hash-verified LXCat Community/Song-Bouwman complete CH4 set. Its four neutral fragmentation channels are contributor-defined but not a validated wet-apparatus yield. H2 and NH3 neutral fragment sources remain unresolved.</p>
    <p className="analysis-warning">The opt-in IST-Lisbon H₂ replacement changes the electron distribution, but its dissociation labels include excited H and cannot yet supply ground-state H chemistry. Charged outcomes are not neutralized in the primary gas branch.</p>
    <p className="analysis-warning">The 1953 profile matches the reported water amount, initial gas partial pressures, and one-week duration. Its discharge current, pulse structure, coupling efficiency, and reduced electric field are explicitly hypothetical because the benchmark sources do not provide a complete deposited-power history.</p>
    <p className="analysis-warning">The Stribling–Miller profile simulates 168 hourly cycles across the published one-hour-to-seven-day window; these are not 168 experimental samples. It matches the CH4/N2/H2O/H2 feed family with H2/CH4 = 1. Absolute feed amounts, plasma controls, and deposited energy remain scenarios; independent pressure, temperature, and resistance-heater calorimetry evidence are needed before comparison.</p>
    <p className="analysis-warning">The coupled-flow sensitivity preset adds a hypothetical 0.1 L condenser contact-water hold-up, generic 0% non-water return, and four numerical substeps per hour. These are declared model assumptions, not reconstructed historical settings or a validated yield.</p>
    {error && <p className="analysis-warning" role="alert">{error}</p>}
    {output && <section ref={resultRef} aria-label="Coupled apparatus result" aria-live="polite">
      <div className="analysis-title"><h3>{output.apparatus_id} · {output.cycles_completed} cycles</h3>
        <button className="secondary-button" onClick={() => downloadJson(output, 'prebiotic-apparatus-result.json')}>Export result</button></div>
      <div className="summary-cards">
        <div><span>Material ledger</span><strong>{output.material_ledger.closed ? 'Closed' : 'Failed'}</strong></div>
        <div><span>Deposited energy</span><strong>{output.energy_ledger.deposited_energy_used_J.toPrecision(5)}</strong><small>J · {output.energy_ledger.deposited_energy_basis.replaceAll('_', ' ')}</small></div>
        <div><span>Yield claim</span><strong>{output.scientific_firewall.quantitative_yield_claim_allowed ? 'Allowed' : 'Blocked'}</strong></div>
        <div><span>Stable gas network</span><strong>{output.plasma_model.heavy_species_inventory_coupled ? 'Active' : 'Pending'}</strong><small>{output.plasma_model.heavy_species_engine ?? 'not connected'}</small></div>
        <div><span>Apparatus integration</span><strong>{output.operating_conditions.operation_mode.replaceAll('_', ' ')}</strong><small>{output.operating_conditions.integration_substeps_per_cycle} substep{output.operating_conditions.integration_substeps_per_cycle === 1 ? '' : 's'} per cycle · GRI source {output.operating_conditions.primary_source_forcing.replaceAll('_', ' ')} · convergence {output.operating_conditions.operator_splitting_convergence_assessed_in_this_run ? 'assessed' : 'pending'}</small></div>
        <div><span>Cold-trap return</span><strong>{(100 * output.operating_conditions.condensate_water_return_fraction).toFixed(0)}% water</strong><small>{(100 * output.operating_conditions.condensate_nonwater_return_fraction).toFixed(0)}% of each non-water species returned · {output.operating_conditions.condenser_phase_transfer_model.replaceAll('_', ' ')}{output.operating_conditions.condenser_contact_water_holdup_L ? ` · ${output.operating_conditions.condenser_contact_water_holdup_L} L hypothetical contact hold-up` : ''}</small></div>
        {output.operating_conditions.simulated_gas_pressure_range_Pa && <div><span>Simulated gas pressure</span><strong>{(output.operating_conditions.simulated_gas_pressure_range_Pa[0] / 1e5).toPrecision(4)}–{(output.operating_conditions.simulated_gas_pressure_range_Pa[1] / 1e5).toPrecision(4)}</strong><small>bar from inventory · declared {(output.operating_conditions.declared_reference_pressure_Pa / 1e5).toPrecision(4)} bar is reference only</small></div>}
        <div><span>Spark-zone exposure</span><strong>{output.operating_conditions.effective_heavy_chemistry_residence_time_s_per_cycle.toPrecision(4)}</strong><small>s per cycle · {output.operating_conditions.effective_heavy_chemistry_residence_time_s_per_integration_step.toPrecision(4)} s/substep · {(100 * output.operating_conditions.discharge_zone_fraction).toPrecision(3)}% of gas loop{output.operating_conditions.discharge_zone_evidence ? ` · ${output.operating_conditions.discharge_zone_evidence.basis}` : ' · full-loop upper bound'}</small></div>
      </div>
      {output.hcn_formation && <>
        <div className="analysis-title"><h3>HCN formation and fate</h3></div>
        <div className="summary-cards">
          <div><span>Formation check</span><strong>{output.hcn_formation.status === 'unseeded_hcn_formed' ? 'Formed in simulation' : output.hcn_formation.status.replaceAll('_', ' ')}</strong><small>not an experimental detection</small></div>
          <div><span>New HCN</span><strong>{output.hcn_formation.net_hcn_generated_mol.toExponential(3)}</strong><small>mol across all compartments</small></div>
          <div><span>In water or collected</span><strong>{output.hcn_formation.final_hcn_plus_cyanide_aqueous_mol.toExponential(3)}</strong><small>mol HCN + CN⁻; gas and wall remain separate</small></div>
          <div><span>Feed carbon in HCN</span><strong>{output.hcn_formation.net_hcn_per_feed_carbon === null ? 'N/A' : `${(100 * output.hcn_formation.net_hcn_per_feed_carbon).toPrecision(3)}%`}</strong><small>scenario conversion, not a validated yield</small></div>
          {output.hcn_formation.stoichiometric_atom_ceiling && <div><span>HCN + CN⁻ atom ceiling</span><strong>{output.hcn_formation.stoichiometric_atom_ceiling.cyanide_family_total_ceiling_mol.toExponential(3)}</strong><small>mol · element-only upper bound, not predicted yield</small></div>}
        </div>
        {output.hcn_formation.aqueous_hydrolysis_screen && <p className="analysis-warning">
          Aqueous HCN hydrolysis screen: {output.hcn_formation.aqueous_hydrolysis_screen.status === 'read_only_step_endpoint_hydrolysis_screen'
            ? `${output.hcn_formation.aqueous_hydrolysis_screen.combined_screen_mol.nominal_mol.toExponential(3)} mol estimated from simulated water snapshots (rate-parameter range ${output.hcn_formation.aqueous_hydrolysis_screen.combined_screen_mol.parameter_low_mol.toExponential(3)}–${output.hcn_formation.aqueous_hydrolysis_screen.combined_screen_mol.parameter_high_mol.toExponential(3)} mol).`
            : `Incomplete: ${output.hcn_formation.aqueous_hydrolysis_screen.unavailable_exposure_count} cyanide-bearing water snapshots lack a supported pH or temperature.`}
          {' '}Read-only screening; this is not subtracted from the predicted amount or a validated yield.
        </p>}
        <p className="analysis-warning">{output.hcn_formation.interpretation}</p>
      </>}
      {lastPlasma && <>
        <div className="analysis-title"><h3>Coupled plasma result</h3></div>
        <div className="summary-cards">
          <div><span>Cycle mixture</span><strong>Connected</strong><small>{(lastPlasma.covered_gas_fraction * 100).toPrecision(4)}% collision coverage</small></div>
          <div><span>Mean electron energy</span><strong>{lastPlasma.mean_electron_energy_eV?.toPrecision(5) ?? 'unavailable'}</strong><small>eV · cycle {lastPlasma.cycle}</small></div>
          <div><span>Potential collision outcomes</span><strong>{lastPlasma.reactive_channel_count ?? 0}</strong><small>{lastPlasma.balanced_reactive_channel_count ?? 0} atom/charge balanced; electronic products may still be unresolved</small></div>
          {lastPlasma.drift_power_diagnostic && <div><span>Electron power check</span><strong>{lastPlasma.drift_power_diagnostic.drift_to_deposited_energy_ratio.toPrecision(4)}x</strong><small>declared-density drift energy / deposited energy; uniform-plasma diagnostic</small></div>}
          {lastPlasma.drift_power_diagnostic?.used_to_scale_source_terms && <div><span>Energy-bound electron density</span><strong>{lastPlasma.electron_number_density_m3?.toExponential(3) ?? 'unavailable'}</strong><small>electrons/m^3; declared {lastPlasma.declared_electron_number_density_m3?.toExponential(3)}; assumes at most all deposited energy drives electrons</small></div>}
          <div><span>Inventory coupling</span><strong>{output.plasma_model.heavy_species_inventory_coupled ? 'Active' : 'Next gate'}</strong><small>source terms are conditional</small></div>
          {lastN2Dissociation && <div><span>N₂ dissociation outcome</span><strong>Electronic states unresolved</strong><small>{lastN2Dissociation.conditional_source_rate_mol_s.toExponential(4)} mol/s conditional total-event rate before energy cap; not injected as two ground-state N atoms</small></div>}
        </div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Target</th><th>Primary channel</th><th>Role</th><th>Conditional source</th><th>Identity</th></tr></thead><tbody>
          {reactiveSources.map((source) => <tr key={source.id}><td>{source.target}</td><td>{source.product_label}</td><td>{source.role}</td><td>{source.conditional_source_rate_mol_s.toExponential(4)} mol/s</td><td>{source.product_state_unresolved ? 'electronic states unresolved' : source.atoms_and_product_charge_balanced ? 'balanced atoms/charge' : source.unresolved_dissociation ? 'fragments unresolved' : 'unresolved'}</td></tr>)}
        </tbody></table></div>
        <p className="analysis-warning">These electron-impact source opportunities use the energy-bound density when the declared density would exceed the deposited-energy budget. That bound assumes all deposited energy could power electron drift; the actual fraction and spark geometry are unknown. Only resolved, energy-limited neutral fragments feed the gas inventory. Source opportunities are not product amounts or yields.</p>
        {firstResolvedO1d && <p><strong>O(¹D) first-collision audit:</strong> if electron-produced O(¹D) has thermalized, at the earliest covered interval {(100 * (firstResolvedO1d.physical_ground_O_fraction ?? 0)).toFixed(1)}% of calculated removal is physical quenching to ground-state O and {(100 * (firstResolvedO1d.reactive_OH_H_fraction ?? 0)).toFixed(1)}% follows OH/H-forming channels; estimated mean lifetime {(firstResolvedO1d.mean_first_collision_lifetime_s ?? 0).toExponential(2)} s. These are conditional fractions, not injected products or HCN yield.</p>}
        {lastO1d?.status === 'conditional_fate_partially_resolved' && <p className="analysis-warning"><strong>Later O(¹D) fate is only partly known:</strong> measured CO/HCN removal rates, where present, contribute to the {(lastO1d.mean_first_collision_lifetime_s ?? 0).toExponential(2)} s lifetime estimate{lastO1d.lifetime_range_from_CO_rate_only_s ? ` (CO-rate parameter range ${lastO1d.lifetime_range_from_CO_rate_only_s[0].toExponential(2)}–${lastO1d.lifetime_range_from_CO_rate_only_s[1].toExponential(2)} s)` : ''}. Of first collisions, {(100 * (lastO1d.unresolved_CO_fate_fraction ?? 0)).toFixed(2)}% have unresolved CO products and {(100 * (lastO1d.unresolved_HCN_reactive_fraction ?? 0)).toFixed(2)}% have unresolved HCN reaction products; the measured HCN quenching-to-ground-O share is counted separately. These are conditional branches, not an HCN yield or full uncertainty range. No products are injected from this audit.</p>}
        {lastO1d?.status === 'unavailable' && <p className="analysis-warning">Later O(¹D) fate is unresolved: {lastO1d.reason}. {lastO1d.unreviewed_bath_species?.length ? `Additional bath species: ${lastO1d.unreviewed_bath_species.slice(0, 8).join(', ')} (${((lastO1d.unreviewed_bath_mole_fraction ?? 0) * 100).toPrecision(3)}% of gas molecules).` : ''} The first-interval estimate must not be reused for the later mixture.</p>}
        {lastPlasma.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
      </>}
      {lastHeavy && <>
        {firstPathwayAudit && <>
          <div className="analysis-title"><h3>First-discharge integrated pathway audit</h3></div>
          <p className="analysis-warning">These are signed reaction amounts integrated over the first modeled discharge, not final-instant rates or measured yields. The audit is diagnostic only; it does not alter the solver or use historical product amounts.</p>
          <div className="analysis-table-wrap"><table><thead><tr><th>Species</th><th>Integrated formation</th><th>Integrated consumption</th><th>Chemical net</th><th>Audit residual</th></tr></thead><tbody>
            {['N2', 'H2CN', 'HCN'].filter((species) => firstPathwayAudit.species[species]).map((species) => {
              const item = firstPathwayAudit.species[species]
              return <tr key={species}><td>{species}</td><td>{item.formation.slice(0, 2).map((reaction) => <div key={reaction.reaction_index}>{reaction.equation}<small>+{reaction.species_net_mol.toExponential(3)} mol</small></div>)}</td><td>{item.consumption.slice(0, 2).map((reaction) => <div key={reaction.reaction_index}>{reaction.equation}<small>{reaction.species_net_mol.toExponential(3)} mol</small></div>)}</td><td>{item.integrated_chemical_net_mol.toExponential(4)} mol</td><td>{item.quadrature_closure_residual_mol.toExponential(2)} mol</td></tr>
            })}
          </tbody></table></div>
          <p className="analysis-warning">Computed with {firstPathwayAudit.method} over {firstPathwayAudit.adaptive_interval_count} intervals. A small audit residual checks numerical integration only; it does not validate the reaction mechanism.</p>
        </>}
        <div className="analysis-title"><h3>Emergent neutral gas species · cycle {lastHeavy.cycle}</h3></div>
        {lastHeavy.collision_energy_allocation && <div className="summary-cards">
          <div><span>Collision energy scale</span><strong>{lastHeavy.collision_energy_allocation.global_collision_energy_scale.toExponential(4)}</strong><small>shared by reactive and nonreactive collision channels</small></div>
          <div><span>Fragment chemistry energy</span><strong>{lastHeavy.collision_energy_allocation.resolved_fragment_energy_used_J.toExponential(4)}</strong><small>J of {lastHeavy.collision_energy_allocation.allocated_collision_energy_J.toExponential(4)} J allocated collision energy</small></div>
          <div><span>Other collision demand</span><strong>{lastHeavy.collision_energy_allocation.nonfragment_or_inventory_limited_energy_J.toExponential(4)}</strong><small>J assigned to excitation, unresolved, or inventory-limited channels</small></div>
          {(lastHeavy.collision_energy_allocation.excluded_fragment_channels?.length ?? 0) > 0 && <div><span>Excluded fragment channels</span><strong>{lastHeavy.collision_energy_allocation.excluded_fragment_channels?.length}</strong><small>Excited, charged, or annotation-inconsistent products are not injected as neutral ground-state fragments</small></div>}
          {(lastHeavy.collision_energy_allocation.overlapping_outcome_raw_threshold_J ?? 0) > 0 && <div><span>Overlapping N2 outcome</span><strong>{lastHeavy.collision_energy_allocation.overlapping_outcome_raw_threshold_J?.toExponential(4)}</strong><small>J event threshold reported but not charged twice to the electron-energy ledger</small></div>}
        </div>}
        <div className="analysis-table-wrap"><table><thead><tr><th>Species</th><th>Net formed this cycle</th><th>Origin</th></tr></thead><tbody>
          {emergedGasSpecies.map(([species, amount]) => <tr key={species}><td>{species}</td><td>{amount.toExponential(5)} mol</td><td>{lastHeavy.mechanism}</td></tr>)}
        </tbody></table></div>
        {shadowComparison && <>
          <div className="analysis-title"><h3>Downstream mechanism shadow</h3></div>
          <div className="summary-cards">
            <div><span>Primary</span><strong>GRI-Mech</strong><small>inventory-driving branch</small></div>
            <div><span>Shadow</span><strong>Pavan 2023</strong><small>{shadowComparison.shadow.retained_heavy_reactions ?? 0} neutral/heavy reactions</small></div>
            <div><span>Charged Pavan state</span><strong>{chargedShadow?.status.endsWith('_complete') ? 'Solved' : 'Diagnostic failed'}</strong><small>{chargedShadow?.retained_two_temperature_reactions ?? 0} ion/electron recombination reactions · {chargedShadow?.source_extents?.length ?? 0} mapped source channels</small></div>
            <div><span>Charge residual</span><strong>{(chargedShadow?.final_charge_residual_mol_e ?? 0).toExponential(3)}</strong><small>mol charge equivalents after this cycle</small></div>
            <div><span>Mapped N2 excitations</span><strong>{shadowComparison.shadow.mapped_excitation_source_extents?.length ?? 0}</strong><small>{(shadowComparison.shadow.mapped_excitation_energy_used_J ?? 0).toExponential(4)} J this cycle</small></div>
            <div><span>Kamarinopoulou 2024</span><strong>{literatureShadow?.status.endsWith('_complete') ? 'Solved' : 'Diagnostic failed'}</strong><small>{literatureShadow?.retained_heavy_reactions ?? 0} reactions · {((literatureShadow?.integrated_duration_s ?? 0) * 1000).toPrecision(3)} ms published domain</small></div>
            <div><span>Independent confidence</span><strong>No</strong><small>Pavan extends GRI-Mech; this is a correlated sensitivity branch</small></div>
          </div>
          {shadowSpecies.length > 0 && <div className="analysis-table-wrap"><table><thead><tr><th>Species</th><th>GRI net change</th><th>Pavan neutral</th><th>Pavan charged</th><th>Kamarinopoulou net change</th><th>Literature / GRI</th></tr></thead><tbody>
            {shadowSpecies.map((species) => {
              const item = shadowComparison.species_comparison?.[species]
              return item ? <tr key={species}><td>{species}</td><td>{item.gri_net_change_mol.toExponential(5)} mol</td><td>{item.pavan_shadow_net_change_mol.toExponential(5)} mol</td><td>{item.pavan_charged_shadow_net_change_mol.toExponential(5)} mol</td><td>{item.literature_hcn_shadow_net_change_mol.toExponential(5)} mol</td><td>{item.literature_hcn_to_gri_ratio === null ? 'undefined' : item.literature_hcn_to_gri_ratio.toPrecision(4)}</td></tr> : null
            })}
          </tbody></table></div>}
          {shadowComparison.shadow.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
          {chargedShadow?.error && <p className="analysis-warning">Charged Pavan solver: {chargedShadow.error}</p>}
          {chargedShadow?.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
          {literatureShadow?.error && <p className="analysis-warning">Literature shadow solver: {literatureShadow.error}</p>}
          {literatureShadow?.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
        </>}
        {propagatedCharged && <>
          <div className="analysis-title"><h3>Propagated charge-balanced branch</h3></div>
          <div className="summary-cards">
            <div><span>Cycle continuity</span><strong>{propagatedCharged.cycles_completed}/{output.cycles_completed}</strong><small>independent Pavan state carried between cycles</small></div>
            <div><span>Branch ledger</span><strong>{propagatedCharged.material_ledger.closed ? 'Closed' : 'Failed'}</strong><small>atoms and electrical charge</small></div>
            <div><span>Branch EEDF</span><strong>{chargedShadow?.branch_plasma_snapshot?.status === 'apparatus_mixture_eedf_solved' ? 'Self-consistent' : 'Unavailable'}</strong><small>re-solved from the branch gas composition each cycle</small></div>
            <div><span>Yield authority</span><strong>Comparison only</strong><small>does not overwrite the primary inventory</small></div>
          </div>
          {output.hcn_formation.charged_branch_aqueous_hydrolysis_screen && <p className="analysis-warning">Charged-branch HCN hydrolysis screen: {output.hcn_formation.charged_branch_aqueous_hydrolysis_screen.status === 'read_only_step_endpoint_hydrolysis_screen'
            ? `${output.hcn_formation.charged_branch_aqueous_hydrolysis_screen.combined_screen_mol.nominal_mol.toExponential(3)} mol conditional loss`
            : `${output.hcn_formation.charged_branch_aqueous_hydrolysis_screen.unavailable_exposure_count} unsupported cyanide-bearing snapshots`}.
            {' '}Read-only comparison; not subtracted from branch inventory.</p>}
          {chargedCapturedPrecursors.length > 0 && <div className="analysis-table-wrap"><table><thead><tr><th>Species</th><th>Charged-branch final water amount</th><th>Primary final water amount</th></tr></thead><tbody>
            {chargedCapturedPrecursors.map(([species, amount]) => {
              const primary = capturedPrecursors.find(([name]) => name === species)?.[1] ?? 0
              return <tr key={species}><td>{species}</td><td>{amount.toExponential(5)} mol</td><td>{primary.toExponential(5)} mol</td></tr>
            })}
          </tbody></table></div>}
          <div className="summary-cards">
            {Object.entries(propagatedCharged.precursor_energy_benchmark.species_comparison).map(([species, item]) => {
              const regression = propagatedCharged.precursor_energy_benchmark.zero_time_extrapolation[species]
              return <div key={species}><span>{species} charged-branch diagnostic</span><strong>{item.simulated_nmol_per_J === null ? 'unavailable' : item.simulated_nmol_per_J.toExponential(4)}</strong><small>nmol/J · reference {item.reference_nmol_per_J}{regression?.available ? ' · zero-time extrapolated' : ' · final inventory'}</small></div>
            })}
          </div>
          {chargedHcnPathways && <div className="analysis-table-wrap"><table><thead><tr><th>Final-cycle HCN direction</th><th>Pavan reaction</th><th>Instantaneous contribution</th></tr></thead><tbody>
            {chargedHcnPathways.formation.map((item) => <tr key={`form:${item.reaction_index}`}><td>Formation</td><td>{item.equation}</td><td>+{item.species_net_rate_mol_s.toExponential(4)} mol/s</td></tr>)}
            {chargedHcnPathways.consumption.map((item) => <tr key={`consume:${item.reaction_index}`}><td>Consumption</td><td>{item.equation}</td><td>{item.species_net_rate_mol_s.toExponential(4)} mol/s</td></tr>)}
          </tbody></table></div>}
          {propagatedCharged.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
        </>}
        {lastHeavy.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
      </>}
      {lastPartition && <>
        <div className="analysis-title"><h3>Cold-trap capture and aqueous inventory · cycle {lastPartition.cycle}</h3></div>
        {lastPartition.water_vapor_capacity && <div className="summary-cards">
          <div><span>Water vapor cap</span><strong>{(lastPartition.water_vapor_capacity.maximum_water_vapor_mol ?? 0).toExponential(4)}</strong><small>mol at pure-water saturation</small></div>
          <div><span>Returned before discharge</span><strong>{lastPartition.water_vapor_capacity.recondensed_water_mol.toExponential(4)}</strong><small>mol supersaturated water</small></div>
          <div><span>Boiler gas pressure</span><strong>{((lastPartition.water_vapor_capacity.ideal_total_pressure_Pa ?? 0) / 1e5).toPrecision(4)}</strong><small>bar ideal-gas inventory · water capped at saturation</small></div>
          {lastPartition.condenser_water_vapor_capacity && <div><span>Cold-end water returned</span><strong>{lastPartition.condenser_water_vapor_capacity.recondensed_water_mol.toExponential(4)}</strong><small>mol before Henry partition</small></div>}
        </div>}
        {capturedPrecursors.length ? <div className="analysis-table-wrap"><table><thead><tr><th>Species</th><th>Final aqueous + isolated collection</th><th>Isolated collection</th></tr></thead><tbody>
          {capturedPrecursors.map(([species, amount]) => <tr key={species}><td>{species}</td><td>{amount.toExponential(5)} mol</td><td>{(output.final_compartments_mol.collected_condensate?.[species] ?? 0).toExponential(5)} mol</td></tr>)}
        </tbody></table></div> : <p>No supported soluble precursor was captured in this run.</p>}
        <p className="analysis-warning">Partitioning is evidence-labelled and mass-conserving, but it is not yet an aqueous reaction or amino-acid result. pH, activity coefficients, hydration equilibria, and competing liquid chemistry are not yet solved.</p>
      </>}
      {output.aqueous_discovery && <>
        <div className="analysis-title"><h3>Blind aqueous structure candidates</h3></div>
        <div className="summary-cards">
          <div><span>Candidate formulas</span><strong>{output.aqueous_discovery.formula_count}</strong></div>
          <div><span>Neutral structures</span><strong>{output.aqueous_discovery.candidate_structure_count}</strong><small>acyclic CHNO gate</small></div>
          <div><span>Reaction templates</span><strong>{output.aqueous_discovery.scientific_firewall.reaction_templates_accessible_to_generator ? 'Used' : 'Blocked'}</strong></div>
          <div><span>Post-generation amino-acid motif</span><strong>{aminoAcidMotifCandidates.length}</strong><small>not a yield</small></div>
        </div>
        {aminoAcidMotifCandidates.length > 0 && <div className="analysis-table-wrap"><table><thead><tr><th>Formula</th><th>Generated structure</th><th>Element sources</th><th>Status</th></tr></thead><tbody>
          {aminoAcidMotifCandidates.slice(0, 12).map((candidate) => <tr key={`${candidate.formula}-${candidate.smiles}`}><td>{candidate.formula}</td><td>{candidate.smiles}</td><td>{candidate.sources.slice(0, 2).map((source) => source.join(' + ')).join('; ')}</td><td>candidate only</td></tr>)}
        </tbody></table></div>}
        {output.aqueous_discovery.limitations.map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
      </>}
      {output.post_discovery_validation && <>
        <div className="analysis-title"><h3>Post-discovery evidence check</h3></div>
        <p className="analysis-warning">{output.post_discovery_validation.hydrolysis_sensitivity.cycles_replayed} apparatus cycles replayed after the blind candidate artifact was frozen. This sidecar cannot alter gas production or candidate discovery.</p>
        <div className="summary-cards">
          <div><span>Frozen discovery</span><strong>{output.post_discovery_validation.discovery_artifact_sha256.slice(0, 12)}</strong><small>SHA-256 prefix</small></div>
          <div><span>Candidate identity</span><strong>{output.post_discovery_validation.candidate_product.name ?? 'not identified'}</strong><small>{output.post_discovery_validation.candidate_product.smiles}</small></div>
            <div><span>Reference aminonitrile</span><strong>{output.post_discovery_validation.reference_equilibrium.estimated_equilibrium_mol_at_298_15_K.toExponential(4)}</strong><small>mol · {output.post_discovery_validation.reference_equilibrium.aqueous_volume_L.toExponential(3)} L water basis</small></div>
          <div><span>Kinetic applicability</span><strong>{output.post_discovery_validation.hydrolysis_sensitivity.applicability.startsWith('within') ? 'Matched' : 'Sensitivity only'}</strong><small>yield claim {output.post_discovery_validation.scientific_firewall.quantitative_yield_claim_allowed ? 'allowed' : 'blocked'}</small></div>
        </div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Hydrolysis scenario</th><th>Rate</th><th>Converted fraction</th><th>Candidate amount</th></tr></thead><tbody>
          {Object.entries(output.post_discovery_validation.hydrolysis_sensitivity.scenarios).map(([label, scenario]) => <tr key={label}><td>{label.replaceAll('_', ' ')}</td><td>{scenario.rate_s_inv.toExponential(4)} s⁻¹</td><td>{(scenario.converted_fraction * 100).toPrecision(4)}%</td><td>{scenario.candidate_product_mol.toExponential(4)} mol</td></tr>)}
        </tbody></table></div>
        {output.post_discovery_validation.limitations.map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
      </>}
      {output.historical_benchmark && <>
        <div className="analysis-title"><h3>Held-out 1953 benchmark</h3></div>
        <div className="summary-cards">
          <div><span>Reproduction status</span><strong>{output.historical_benchmark.validated_reproduction ? 'Validated' : 'Not yet reproduced'}</strong><small>{output.historical_benchmark.status.replaceAll('_', ' ')}</small></div>
          <div><span>Historical conditions</span><strong>{output.historical_benchmark.historical_conditions_matched ? 'Matched' : 'Different'}</strong><small>water, pressure, and one-week duration</small></div>
          <div><span>Enumeration coverage</span><strong>{(100 * output.historical_benchmark.enumeration_coverage).toPrecision(3)}%</strong><small>present somewhere in broad list · not a prediction score</small></div>
          <div><span>Benchmark leakage</span><strong>{output.historical_benchmark.scientific_firewall.benchmark_observations_available_to_discovery ? 'Detected' : 'Blocked'}</strong><small>no fitting to historical yields</small></div>
        </div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Historical observation</th><th>Blind candidate</th><th>Captured precursors</th><th>Scientific status</th></tr></thead><tbody>
          {output.historical_benchmark.observation_comparison.map((item) => <tr key={item.name}><td>{item.name}</td><td>{item.candidate_recovered ? item.candidate_smiles.join(' / ') : 'not recovered'}</td><td>{Object.entries(item.captured_precursors_mol).map(([species, amount]) => `${species} ${amount.toExponential(3)} mol`).join(' · ') || 'pathway not implemented'}</td><td>{item.pathway_status.replaceAll('_', ' ')}</td></tr>)}
        </tbody></table></div>
        {output.historical_benchmark.remaining_blockers.map((blocker) => <p className="analysis-warning" key={blocker}>{blocker}</p>)}
      </>}
      {output.precursor_energy_benchmark && <>
        <div className="analysis-title"><h3>Measured precursor energy-yield gate</h3></div>
        <div className="summary-cards">
          <div><span>Comparison status</span><strong>{output.precursor_energy_benchmark.comparable ? 'Comparable' : 'Diagnostic only'}</strong><small>{output.precursor_energy_benchmark.status.replaceAll('_', ' ')}</small></div>
          <div><span>Feed composition</span><strong>{output.precursor_energy_benchmark.composition_matched ? 'Matched' : 'Different'}</strong><small>published CH4/N2/H2O/H2 condition</small></div>
          <div><span>Zero-time basis</span><strong>{output.precursor_energy_benchmark.feed_checks.zero_time_yield_basis_matched ? 'Extrapolated' : 'Incomplete'}</strong><small>{output.precursor_energy_benchmark.zero_time_extrapolation.HCN?.sample_count ?? 0} samples in window</small></div>
          <div><span>Energy basis</span><strong>{output.precursor_energy_benchmark.feed_checks.calorimetric_energy_basis_matched ? 'Calorimetric' : 'Electrical scenario'}</strong><small>{output.precursor_energy_benchmark.energy_basis.replaceAll('_', ' ')}</small></div>
          <div><span>Pressure closure</span><strong>{output.precursor_energy_benchmark.pressure_self_consistent === true ? 'Consistent' : output.precursor_energy_benchmark.pressure_self_consistent === false ? 'Mismatch' : 'Unassessed'}</strong><small>inventory-derived gas pressure vs declared reference; a mismatch blocks comparison</small></div>
          {Object.entries(output.precursor_energy_benchmark.species_comparison).map(([species, item]) => {
            const regression = output.precursor_energy_benchmark.zero_time_extrapolation[species]
            return <div key={species}><span>{species} energy yield</span><strong>{item.simulated_nmol_per_J === null ? 'unavailable' : item.simulated_nmol_per_J.toExponential(4)}</strong><small>nmol/J · reference {item.reference_nmol_per_J} · {item.observable_basis.replaceAll('_', ' ')}{regression?.available && regression.intercept_standard_error_nmol_per_J !== null ? ` · intercept SE ${regression.intercept_standard_error_nmol_per_J.toExponential(3)} · R² ${regression.r_squared?.toPrecision(4) ?? 'unavailable'}` : ''}</small></div>
          })}
        </div>
        {output.precursor_energy_benchmark.illustrative_sparse_sampling_sensitivity && <>
          <div className="analysis-table-wrap"><table><thead><tr><th>Precursor</th><th>Hourly fit</th><th>Illustrative sparse fit</th><th>Sparse fit quality</th></tr></thead><tbody>
            {Object.entries(output.precursor_energy_benchmark.zero_time_extrapolation).map(([species, hourly]) => {
              const sparse = output.precursor_energy_benchmark.illustrative_sparse_sampling_sensitivity?.zero_time_extrapolation[species]
              return <tr key={species}><td>{species}</td><td>{hourly.intercept_nmol_per_J?.toPrecision(5) ?? 'unavailable'} nmol/J</td><td>{sparse?.intercept_nmol_per_J?.toPrecision(5) ?? 'unavailable'} nmol/J</td><td>{sparse?.r_squared?.toPrecision(4) ?? 'unavailable'} R² · {sparse?.sample_count ?? 0} points</td></tr>
            })}
          </tbody></table></div>
          <p className="analysis-warning">Sparse resampling uses an illustrative 1-hour, 2/4/6-day, 7-day grid after simulation. It is not the original experimental schedule or an independent yield measurement; the exact sampling times and fit weighting are unavailable in the accessible abstracts.</p>
        </>}
        <p className="analysis-warning">{output.precursor_energy_benchmark.comparable
          ? 'Feed, zero-time observable, and calorimetric energy basis are comparable. A validation pass remains blocked until the published regression weighting and reference uncertainty are available.'
          : output.precursor_energy_benchmark.composition_matched
          ? 'The feed composition matches the published CH4/N2/H2O/H2 case, but the measurement basis is incomplete; the published values used samples from one hour to seven days, zero-time extrapolation, and calorimetry.'
          : 'This feed does not match the later N2 energy-yield experiment. The displayed simulated values are final-inventory diagnostics, while the published values were extrapolated to zero time.'}</p>
        {output.precursor_energy_benchmark.remaining_blockers.map((blocker) => <p className="analysis-warning" key={blocker}>{blocker}</p>)}
      </>}
      <div className="analysis-title"><h3>Frozen-candidate model comparison</h3></div>
      <p className="analysis-warning">This optional background calculation applies an amino-carboxylic functional-class filter after blind enumeration, then ranks those structures independently with MACE-MP, ANI-2x, AIMNet2, and AIMNet2-RXN. It is a post-discovery class comparison, not a global blind product predictor. Energies remain formula-local vacuum values—not solution probabilities or yields.</p>
      <button className="secondary-button" disabled={rankingJob?.status === 'queued' || rankingJob?.status === 'running'} onClick={rankFrozenCandidates}>{rankingJob?.status === 'queued' || rankingJob?.status === 'running' ? `${rankingJob.message} · ${rankingJob.progress}%` : 'Compare frozen candidates on all models'}</button>
      {rankingError && <p className="analysis-warning" role="alert">{rankingError}</p>}
      {rankingJob?.result && <>
        <div className="summary-cards">
          <div><span>Frozen candidates</span><strong>{rankingJob.result.candidate_count}</strong><small>{rankingJob.result.candidate_filter}</small></div>
          <div><span>Models completed</span><strong>{rankingJob.result.engines_succeeded.length}/{rankingJob.result.engines_requested.length}</strong><small>{rankingJob.result.engines_succeeded.join(' · ')}</small></div>
          <div><span>Discovery hash</span><strong>{rankingJob.result.discovery_artifact_sha256.slice(0, 12)}</strong><small>verified before ranking</small></div>
          <div><span>Yield claim</span><strong>{rankingJob.result.scientific_firewall.pathway_or_yield_claimed ? 'Claimed' : 'Blocked'}</strong><small>ranking is not probability</small></div>
        </div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Formula</th><th>Frozen structure</th><th>Mean local rank</th><th>Rank range</th><th>Relative-energy range</th></tr></thead><tbody>
          {rankingJob.result.formula_local_consensus.map((item) => <tr key={`${item.formula}:${item.smiles}`}><td>{item.formula}</td><td>{item.smiles}</td><td>{item.mean_formula_local_rank.toPrecision(3)}</td><td>{item.rank_range.join('–')}</td><td>{item.relative_energy_range_kJ_mol.map((value) => value.toPrecision(4)).join('–')} kJ/mol</td></tr>)}
        </tbody></table></div>
        <p className="analysis-warning">R/S mirror pairs use reflected starting coordinates and tied ranks because exact enantiomers are energy-degenerate in this achiral vacuum calculation.</p>
        <button className="secondary-button" disabled={evaluationBusy} onClick={revealHeldOutRanking}>{evaluationBusy ? 'Revealing held-out identities…' : rankingJob.evaluation ? 'Re-run frozen held-out score' : 'Freeze prediction, then reveal held-out identities'}</button>
        {rankingJob.evaluation && <>
          <div className="summary-cards">
            <div><span>Prediction artifact</span><strong>{rankingJob.evaluation.prediction_artifact_sha256.slice(0, 12)}</strong><small>frozen before benchmark read</small></div>
            <div><span>Class coverage</span><strong>{(100 * rankingJob.evaluation.class_filter_coverage).toPrecision(3)}%</strong><small>conditional on post-discovery class filter</small></div>
            <div><span>Top-{rankingJob.evaluation.top_k} recall</span><strong>{(100 * rankingJob.evaluation.top_k_formula_local_recall).toPrecision(3)}%</strong><small>formula-local rank · not yield</small></div>
            <div><span>Reproduction claim</span><strong>{rankingJob.evaluation.validated_reproduction ? 'Validated' : 'Not validated'}</strong><small>pathways and quantitative yields still required</small></div>
          </div>
          <div className="analysis-table-wrap"><table><thead><tr><th>Held-out observation</th><th>Frozen hit</th><th>Best local rank</th><th>Top-{rankingJob.evaluation.top_k}</th></tr></thead><tbody>
            {rankingJob.evaluation.observation_comparison.map((item) => <tr key={item.name}><td>{item.name}</td><td>{item.candidate_smiles.join(' / ') || 'not in class set'}</td><td>{item.best_mean_formula_local_rank?.toPrecision(3) ?? 'not ranked'}</td><td>{item.top_k_recovered ? 'yes' : 'no'}</td></tr>)}
          </tbody></table></div>
          {rankingJob.evaluation.limitations.map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
        </>}
        {rankingJob.result.engine_failures.map((failure) => <p className="analysis-warning" key={failure.engine}>{failure.engine}: {failure.error}</p>)}
        {rankingJob.result.limitations.map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
      </>}
      <h3>Selected circulation cycles</h3><div className="analysis-table-wrap"><table><thead><tr><th>Cycle</th><th>Intervention</th><th>Evaporated</th><th>Condensed</th><th>Electrical energy</th></tr></thead><tbody>
        {shownCycles.map((cycle) => <tr key={cycle.cycle}><td>{cycle.cycle}<small className="record-note">{cycle.start_time_s} → {cycle.end_time_s} s</small></td><td>{cycle.intervention_notes.join('; ') || 'none'}</td><td>{nonzero(cycle.evaporated_mol)}</td><td>{nonzero(cycle.condensed_mol)}</td><td>{cycle.deposited_energy_used_J.toPrecision(4)} J</td></tr>)}
      </tbody></table></div>
      <h3>Final inventory by compartment</h3>
      {Object.entries(output.final_compartments_mol).map(([compartment, values]) => <p key={compartment}><strong>{compartment}</strong>: {nonzero(values)}</p>)}
      {output.warnings.map((warning) => <p className="analysis-warning" key={warning}>{warning}</p>)}
    </section>}
    <section>
      <h3>Cantera electron-collision dataset gate</h3>
      <p>Validate the energy grid, nonnegative cross sections, provenance and Cantera 3.2 rate-object compatibility. Passing this gate does not calculate a rate without an electron-energy distribution.</p>
      <textarea className="network-editor compact" aria-label="Electron collision cross-section JSON" spellCheck={false}
        value={collisionDefinition} onChange={(event) => { setCollisionDefinition(event.target.value); setCollisionResult(null); setError(null) }} />
      <button className="secondary-button" disabled={busy} onClick={validateCollision}>Validate cross section with Cantera</button>
      {collisionResult && <div className="summary-cards">
        <div><span>Cantera object</span><strong>{collisionResult.cantera_rate_type}</strong><small>v{collisionResult.cantera_version}</small></div>
        <div><span>Data points</span><strong>{collisionResult.points}</strong><small>{collisionResult.energy_range_eV.join('–')} eV</small></div>
        <div><span>First nonzero</span><strong>{collisionResult.first_nonzero_energy_eV ?? 'none'}</strong><small>eV</small></div>
        <div><span>Rate-ready</span><strong>{collisionResult.scientific_firewall.usable_for_rate_without_eedf ? 'Yes' : 'Needs EEDF'}</strong></div>
      </div>}
      <h3>Electron-energy distribution → rate coefficient</h3>
      <p>Combine one collision curve with an explicit EEDF and electron number density. The result is a conditional ⟨σv⟩ and target lifetime, not a product yield.</p>
      <textarea className="network-editor compact" aria-label="Electron collision rate JSON" spellCheck={false}
        value={rateDefinition} onChange={(event) => { setRateDefinition(event.target.value); setRateResult(null); setError(null) }} />
      <button className="secondary-button" disabled={busy} onClick={calculateCollisionRate}>Calculate conditional collision rate</button>
      {rateResult && <>
        <div className="summary-cards">
          <div><span>⟨σv⟩</span><strong>{rateResult.rate_coefficient_m3_s.toExponential(4)}</strong><small>m³/s</small></div>
          <div><span>Mean electron energy</span><strong>{rateResult.eedf.mean_energy_eV.toPrecision(5)}</strong><small>eV · {rateResult.eedf.mode}</small></div>
          <div><span>Target lifetime</span><strong>{rateResult.target_lifetime_s === null ? '∞' : rateResult.target_lifetime_s.toExponential(4)}</strong><small>s at declared nₑ</small></div>
          <div><span>Yield claim</span><strong>{rateResult.scientific_firewall.quantitative_yield_claim_allowed ? 'Allowed' : 'Blocked'}</strong><small>needs full mechanism + calibration</small></div>
        </div>
        {rateResult.limitations.map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
      </>}
      <h3>Self-consistent mixture EEDF</h3>
      <p>Solve the steady-state two-term electron Boltzmann equation from the gas composition, gas temperature, and reduced electric field. The default hash-verified hybrid set combines Morgan/LXCat H2, H2O, CH4, and NH3 data with Cantera's packaged Phelps N2 data and the reviewed QST/Tabata neutral N2 dissociation outcome.</p>
      <textarea className="network-editor compact" aria-label="Boltzmann EEDF JSON" spellCheck={false}
        value={boltzmannDefinition} onChange={(event) => { setBoltzmannDefinition(event.target.value); setBoltzmannResult(null); setError(null) }} />
      <div className="button-row">
        <button className="primary-action" disabled={busy} onClick={solveBoltzmannEEDF}>{busy ? 'Solving...' : 'Solve mixture EEDF'}</button>
        <button className="secondary-button" disabled={busy} onClick={() => { setBoltzmannDefinition(JSON.stringify(boltzmannExample, null, 2)); setBoltzmannResult(null) }}>Reset five-gas example</button>
        <button className="secondary-button" disabled={busy} onClick={() => { setBoltzmannDefinition(JSON.stringify({ ...boltzmannExample, collision_set: communityCh4CollisionSet }, null, 2)); setBoltzmannResult(null); setError(null) }}>Use reviewed methane set</button>
        <button className="secondary-button" disabled={busy || !capabilities?.boltzmann_eedf.collision_sets?.some((item) => item.dataset_id === istH2CommunityCh4CollisionSet && item.available)} onClick={() => { setBoltzmannDefinition(JSON.stringify({ ...boltzmannExample, collision_set: istH2CommunityCh4CollisionSet }, null, 2)); setBoltzmannResult(null); setError(null) }}>Use opt-in IST H₂ set</button>
        <button className="secondary-button" disabled={busy || !capabilities?.boltzmann_eedf.collision_sets?.some((item) => item.dataset_id === istCo2O2CollisionSet && item.available)} onClick={() => { setBoltzmannDefinition(JSON.stringify(carbonOxygenEedfExample, null, 2)); setBoltzmannResult(null); setError(null) }}>Load CO₂/O₂ EEDF example</button>
      </div>
      <p className="analysis-warning">The IST-Lisbon H₂ contributor reports 1–200 Td coverage, but the default BOLOS grid can produce non-finite low-field rates and now rejects them. Its labelled H₂ fragments include excited hydrogen; without quenching kinetics they do not become ground-state H radicals or a validated product yield.</p>
      <p className="analysis-warning">The opt-in CO₂/O₂ example uses private, hash-verified LXCat complete sets at 200 Td. It solves electron collisions only. Excited CO and O products, ions, discharge exposure, and aqueous chemistry are not validated; this does not predict an HCN amount.</p>
      {boltzmannResult && <>
        <div className="summary-cards">
          <div><span>Mean electron energy</span><strong>{boltzmannResult.eedf.mean_energy_eV.toPrecision(5)}</strong><small>eV at {boltzmannResult.conditions.reduced_electric_field_Td} Td</small></div>
          <div><span>Electron temperature</span><strong>{boltzmannResult.eedf.electron_temperature_K.toExponential(3)}</strong><small>K equivalent</small></div>
          <div><span>Convergence</span><strong>{boltzmannResult.solver.converged ? 'Converged' : 'Failed'}</strong><small>{boltzmannResult.solver.iterations} iterations</small></div>
          <div><span>EEDF input coverage</span><strong>{boltzmannResult.collision_set.coverage.all_required_targets_complete ? `${Object.keys(boltzmannResult.collision_set.coverage.targets).length} / ${Object.keys(boltzmannResult.collision_set.coverage.targets).length} supported` : 'Incomplete'}</strong><small>{boltzmannResult.collision_set.coverage.process_count} imported processes; not chemical-fragment coverage</small></div>
          {n2Dissociation && <div><span>N₂ dissociation outcome</span><strong>Electronic states unresolved</strong><small>{n2Dissociation.rate_coefficient_m3_s.toExponential(4)} m3/s conditional total-outcome coefficient; not a ground-state N source</small></div>}
        </div>
        {boltzmannResult.collision_set.coverage.any_unresolved_dissociation_labels && <p className="analysis-warning">Some dissociation curves identify electron energy loss but not the resulting fragments. {Object.entries(boltzmannResult.collision_set.coverage.targets).filter(([, target]) => target.unresolved_dissociation_channel_count > 0).map(([name, target]) => `${name}: ${target.unresolved_dissociation_channel_count}`).join('; ')} unresolved channel(s) cannot directly seed the downstream chemistry.</p>}
        {!boltzmannResult.collision_set.coverage.all_required_targets_have_explicit_neutral_fragment_source && <p className="analysis-warning">No explicit neutral-fragment source is available for {Object.entries(boltzmannResult.collision_set.coverage.targets).filter(([, target]) => !target.has_explicit_neutral_fragment_source).map(([name]) => name).join(', ')}. Electron-distribution support does not make their radical chemistry complete.</p>}
        {boltzmannResult.collision_set.coverage.all_required_targets_have_resolved_ground_state_fragment_source === false && <p className="analysis-warning">No usable ground-state neutral-fragment source is available for {Object.entries(boltzmannResult.collision_set.coverage.targets).filter(([, target]) => !target.has_resolved_ground_state_fragment_source).map(([name]) => name).join(', ')}. A source label can name fragments but still require excited-state or ion chemistry before injection.</p>}
        <div className="analysis-title"><h3>Fastest calculated electron channels</h3>
          <button className="secondary-button" onClick={() => downloadJson(boltzmannResult, 'self-consistent-eedf-result.json')}>Export EEDF</button></div>
        <div className="analysis-table-wrap"><table><thead><tr><th>Target</th><th>Channel</th><th>Contributor product label</th><th>Rate coefficient</th></tr></thead><tbody>
          {[...boltzmannResult.rate_coefficients].sort((left, right) => right.rate_coefficient_m3_s - left.rate_coefficient_m3_s).slice(0, 8).map((rate) =>
            <tr key={rate.id}><td>{rate.target}</td><td>{rate.kind}</td><td>{rate.product_label}</td><td>{rate.rate_coefficient_m3_s.toExponential(4)} m3/s</td></tr>)}
        </tbody></table></div>
        <p>{boltzmannResult.conditions.note}</p>
        <p><strong>Independent engine:</strong> {boltzmannResult.independent_reference.available_for_this_mixture
          ? `${boltzmannResult.independent_reference.engine}, mean ${boltzmannResult.independent_reference.mean_energy_eV?.toPrecision(5)} eV`
          : `${boltzmannResult.independent_reference.engine}: ${boltzmannResult.independent_reference.reason ?? `unsupported targets ${boltzmannResult.independent_reference.unsupported_targets?.join(' / ') ?? 'unknown'}`}`}</p>
        {boltzmannResult.same_collision_set_reference?.available_for_this_mixture && <>
          <p><strong>Same-data engine check:</strong> {boltzmannResult.same_collision_set_reference.engine} gives {boltzmannResult.same_collision_set_reference.mean_energy_eV?.toPrecision(5)} eV on the same collision curves; BOLOS differs by {((boltzmannResult.same_collision_set_reference.relative_mean_energy_difference ?? 0) * 100).toFixed(2)}%. This checks numerical agreement, not experimental swarm data or yield.</p>
          {boltzmannResult.same_collision_set_reference.rate_comparison && <>
            <p><strong>Collision-rate check:</strong> {boltzmannResult.same_collision_set_reference.rate_comparison.material_within_factor_ten_count} of {boltzmannResult.same_collision_set_reference.rate_comparison.material_process_count} material channels agree within a factor of ten. Rates, not just mean electron energy, are compared; this is still not an experimental validation.</p>
            <div className="analysis-table-wrap"><table><thead><tr><th>Target</th><th>Reactive channel</th><th>BOLOS</th><th>Cantera EEDF</th><th>Difference</th></tr></thead><tbody>
              {boltzmannResult.same_collision_set_reference.rate_comparison.largest_reactive_channel_by_target.map((row) => <tr key={`${row.target}-${row.kind}-${row.product_label}`}><td>{row.target}</td><td>{row.product_label}</td><td>{row.bolos_rate_coefficient_m3_s.toExponential(3)}</td><td>{row.cantera_eedf_rate_coefficient_m3_s.toExponential(3)}</td><td>{row.factor_difference === null ? 'unresolved' : `${row.factor_difference.toFixed(2)}x`}</td></tr>)}
            </tbody></table></div>
          </>}
        </>}
        {boltzmannResult.same_collision_set_reference && !boltzmannResult.same_collision_set_reference.available_for_this_mixture && <p className="analysis-warning">Same-data engine check unavailable: {boltzmannResult.same_collision_set_reference.reason}</p>}
        {[...boltzmannResult.numerical_warnings, ...boltzmannResult.limitations].map((limitation) => <p className="analysis-warning" key={limitation}>{limitation}</p>)}
      </>}
    </section>
  </div>
}
