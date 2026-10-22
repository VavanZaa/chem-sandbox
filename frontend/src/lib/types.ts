export type SubstanceId = 'sodium_bicarbonate' | 'acetic_acid' | 'water'
export type Effort = 'low' | 'medium' | 'high'

export interface SubstanceInput {
  id: SubstanceId
  amount: number
  unit: 'g' | 'mol' | 'mL'
}

export interface ExperimentDefinition {
  experiment_id: string
  mode: 'bicarbonate' | 'general'
  molecules: Array<{ name: string; smiles: string; amount_mol: number; copies: number }>
  compute: { engine: 'all_compatible' | 'mace_mp' | 'ani2x' | 'aimnet2' | 'aimnet2_rxn'; device: 'auto' | 'cuda' | 'cpu' }
  substances: SubstanceInput[]
  environment: {
    temperature_K: number
    pressure_atm: number
    solvent: 'water' | 'vacuum' | 'other'
    mixing: 'gentle' | 'standard' | 'vigorous'
  }
  search: {
    effort: Effort
    runs: number
    seed: number
    max_candidates_per_state: number
    optimization_steps: number
    candidate_optimization_steps: number
    max_depth: number
    max_evaluations_per_run: number
    max_states_per_run: number
    frontier_width: number
  }
}

export interface GraphNode {
  id: string
  label: string
  formula: string
  species: SpeciesInfo[]
  interpretation_quality: 'resolved' | 'uncertain'
  energy_ev: number
  relative_energy_ev: number
  occurrences: number
  runs_seen: number
  run_frequency: number
  starting_runs: number
  terminal_runs: number
  terminal_run_frequency: number
  converged: boolean
  max_force: number | null
  structure_reference: string | null
  role: 'reactant' | 'intermediate' | 'endpoint'
  endpoint_reasons: Record<string, number>
  comparison_key: string
}

export interface SpeciesInfo {
  formula: string
  display_formula: string
  name: string
  inferred_charge: number | null
  resolved: boolean
  note: string | null
  atom_indices: number[]
  canonical_smiles: string | null
  isomeric_smiles: string | null
  stereochemistry: string | null
  stereo_assignments: string[]
  stereo_status: 'specified' | 'unassigned' | 'achiral' | 'indeterminate' | 'unknown'
  aromatic: boolean | null
  identity_basis: string
}

export interface StructureData {
  filename: string
  atoms: Array<{ index: number; element: string; position: [number, number, number] }>
  bonds: Array<{ source: number; target: number }>
  species: SpeciesInfo[]
  note: string
}

export interface GraphEdge {
  id: string
  source: string
  target: string
  transition_count: number
  run_count: number
  search_frequency: number
  run_frequency: number
  delta_energy_ev: number
  energy_assessment: 'favorable' | 'unfavorable' | 'approximately_neutral'
  barrier_ev: number | null
  method: string
  bonds_created: number[][]
  bonds_removed: number[][]
  comparison_key: string
}

export interface DiscoveryComparison {
  mode: 'independent_searches'
  successful_engine_count: number
  engine_runs: Array<{ engine: string; status: string; model: ExperimentResult['model']; result_job_id: string;
    states: number; paths: number; runtime_seconds: number; failed_optimizations: number }>
  failed_engines: Array<{ engine: string; status: string; reason: string }>
  excluded_engines: Array<{ engine: string; supported: boolean; reason: string }>
  states: Array<{ comparison_key: string; labels: string[]; formulas: string[]; engines: string[];
    model_count: number; agreement_fraction: number }>
  paths: Array<{ comparison_key: string; labels: string[]; engines: string[]; model_count: number;
    agreement_fraction: number; delta_energy_ev_by_engine: Record<string, number> }>
  calibrated_confidence: null
  warnings: string[]
}

export interface ExperimentResult {
  experiment_id: string
  status: 'completed' | 'model_failure'
  model: { engine: string; model: string; device: string; dtype: string; domain_note: string }
  graph: { nodes: GraphNode[]; edges: GraphEdge[] }
  statistics: {
    search_runs: number
    valid_optimized_states: number
    failed_optimizations: number
    rejected_candidate_states: number
    unique_states: number
    candidate_pathways: number
    co2_containing_runs: number
    co2_run_frequency: number
    co2_terminal_runs: number
    co2_terminal_frequency: number
    disclaimer: string
    accounted_runs: number
    candidate_evaluations: number
    initialization_evaluations: number
    run_status_counts: Record<string, number>
    endpoint_reason_counts: Record<string, number>
  }
  scientific_notes: string[]
  errors: Array<Record<string, unknown>>
  runtime_seconds: number
  barriers_calculated: boolean
  run_outcomes: RunOutcome[]
  stoichiometry: Stoichiometry | null
  experiment?: ExperimentDefinition
  search_settings: ExperimentDefinition['search']
  discovery_comparison?: DiscoveryComparison
  display_engine?: string
}

export interface RunOutcome {
  run: number
  status: string
  evaluations: number
  initialization_evaluations: number
  discovered_states: number
  expanded_states: number
  attempt_counts: Record<string, number>
  endpoints: Array<{ state_id: string | null; depth: number; reason: string }>
}

export interface Stoichiometry {
  equation: string
  limiting_reagent: string
  co2_mol: number
  co2_mass_g: number
  co2_ideal_gas_L: number
  sodium_acetate_mol: number
  water_formed_mol: number
  excess_bicarbonate_mol: number
  excess_acid_mol: number
  temperature_K: number
  pressure_atm: number
  assumptions: string[]
}

export interface RateEvidence {
  edge_id: string
  rate_s_inv?: number
  activation_free_energy_ev?: number
  temperature_K: number
  source: string
  basis: 'hypothetical' | 'experimental' | 'external_calculation'
}

export interface KineticsRequest {
  initial_state_id: string
  temperature_K: number
  duration_s: number
  rates: RateEvidence[]
  acknowledge_conditional_model: true
}

export interface KineticsResult {
  status: 'conditional'
  samples: Array<{ time_s: number; fractions: Record<string, number> }>
  warnings: string[]
  rates: Array<RateEvidence & { computed_rate_s_inv: number }>
}

export interface ExperimentJob {
  id: string
  status: 'queued' | 'running' | 'completed' | 'failed'
  stage: string
  progress: number
  message: string
  result: ExperimentResult | null
  error: string | null
}

export type Selection = { kind: 'node'; value: GraphNode } | { kind: 'edge'; value: GraphEdge }

export interface ComputeStatus {
  torch_version: string
  cuda_available: boolean
  cuda_runtime: string | null
  gpus: Array<{ name: string; memory_GB: number }>
  reason: string
  engines: Array<{ id: string; name: string; installed: boolean; gpu: boolean; domain: string }>
}

export interface EfficiencyMetrics {
  scenario: string; conversion: number; target_yield: number | null
  target_formed_mol: number; target_product_share: number | null; final_moles: Record<string, number>
}
export interface EfficiencyResult {
  models: Array<{ model_id: string; basis: string; nominal: EfficiencyMetrics; scenarios: EfficiencyMetrics[];
    stoichiometric_ceiling_mol: number; yield_sensitivity_range: [number, number] | null;
    samples: Array<{ time_s: number; moles: Record<string, number> }> }>
  warnings: string[]
}
export interface RefinementResult {
  records: Array<{ engine: string; edge_id: string; status: string; error?: string; note?: string;
    electronic_barrier_ev?: number; delta_energy_ev?: number; metadata?: { device: string };
    eligible_for_exploratory_kinetics: boolean }>
  warnings: string[]
  unrefined_edges: number
}
export interface RefinementJob extends Omit<ExperimentJob, 'result'> { result: RefinementResult | null }
export interface PathPopulationResult {
  status: string; warning: string
  predictions: Array<{ engine: string; state_conversion: number; target_state_population: number }>
  unavailable: Array<{ engine: string; reason: string }>
}

export interface ProtocolStageResult {
  id: string; name: string; operation: string; start_time_s: number; end_time_s: number
  temperature_K: number; volume_before_L: number; volume_after_L: number
  active_reactions: string[]; additions_mol: Record<string, number>; removed_mol: Record<string, number>
  inventory_final_mol: Record<string, number>; reaction_extents_mol: Record<string, number>
}
export interface ProtocolResult {
  status: string; protocol_id: string; total_duration_s: number
  species_identity: Array<{ id: string; name: string; formula: string; isomeric_smiles: string | null;
    stereo_status: string; stereochemistry: string[] }>
  stages: ProtocolStageResult[]
  final_amounts_mol: Record<string, number>
  material_ledger: { closed: boolean; residual: Record<string, number>; added_mol: Record<string, number>;
    removed_mol: Record<string, number> }
  sbml: { validation: { status: string; errors: number; libsbml_version: string }; sha256: string;
    network_only: boolean; document: string | null }
  roadrunner_cross_check: { status: string; engine?: string; version?: string; max_abs_difference_mol?: number;
    tolerance_mol?: number; error?: string }
  warnings: string[]
}

export interface PrebioticCycleResult {
  cycle: number; start_time_s: number; end_time_s: number
  intervention_notes: string[]
  additions_mol: Record<string, number>; sampled_mol: Record<string, number>
  evaporated_mol: Record<string, number>; condensed_mol: Record<string, number>
  returned_mol: Record<string, number>; wall_deposition_mol: Record<string, number>
  reaction_extents_mol: Record<string, number>; deposited_electrical_energy_J: number
  deposited_energy_used_J: number; deposited_energy_basis: string
  plasma_snapshot?: ApparatusPlasmaSnapshot | null
  heavy_species_snapshot?: HeavySpeciesSnapshot | null
  phase_partition?: PhasePartitionCycle | null
  compartments_mol: Record<string, Record<string, number>>
}

export interface PhasePartitionResult {
  status: string; temperature_K?: number; pressure_Pa?: number; dataset_id?: string
  transferred_to_water_mol: Record<string, number>
  entries_used: Array<{ species_id: string; kH_mol_kg_bar: number; water_mass_kg: number; source: string; selection: string }>
  limitations?: string[]
}

export interface PhasePartitionCycle {
  cycle?: number; boiler: PhasePartitionResult | null; condenser: PhasePartitionResult | null
  water_vapor_capacity?: {
    status: string; temperature_K?: number; gas_volume_L?: number; saturation_pressure_Pa?: number
    maximum_water_vapor_mol?: number; water_partial_pressure_Pa?: number; ideal_total_pressure_Pa?: number; recondensed_water_mol: number
    method?: string; limitations?: string[]
  } | null
  condenser_water_vapor_capacity?: {
    status: string; temperature_K?: number; gas_volume_L?: number; saturation_pressure_Pa?: number
    maximum_water_vapor_mol?: number; water_partial_pressure_Pa?: number; ideal_total_pressure_Pa?: number; recondensed_water_mol: number
    method?: string; limitations?: string[]
  } | null
}

export interface PlasmaSourceTerm {
  id: string; target: string; kind: string; product_label: string; role: string; reactive: boolean
  atoms_and_product_charge_balanced: boolean; unresolved_dissociation: boolean
  product_state_unresolved?: boolean
  conditional_source_rate_mol_s: number; uncoupled_active_window_opportunity_mol: number
  products: Array<{ label: string; formula: string; composition: Record<string, number>; charge: number; excited: boolean }>
}

export interface ApparatusPlasmaSnapshot {
  status: string; mole_fractions: Record<string, number>; covered_gas_fraction: number
  electron_number_density_m3?: number; declared_electron_number_density_m3?: number
  drift_power_diagnostic?: {
    drift_to_deposited_energy_ratio: number
    estimated_electron_drift_energy_J: number
    deposited_energy_J: number
    electron_density_if_all_deposit_went_to_drift_m3: number | null
    effective_electron_density_upper_bound_m3: number
    used_to_scale_source_terms: boolean
    electron_energy_fraction_known: boolean
  }
  active_discharge_duration_s?: number
  discharge_zone_volume_L?: number; discharge_zone_fraction?: number
  whole_loop_mean_residence_time_s?: number; residence_time_method?: string
  mean_electron_energy_eV?: number; electron_temperature_K?: number
  reactive_channel_count?: number; balanced_reactive_channel_count?: number
  o1d_fate_diagnostic?: {
    status: string; reason?: string; unreviewed_bath_species?: string[]
    unreviewed_bath_mole_fraction?: number
    physical_ground_O_fraction?: number; reactive_OH_H_fraction?: number
    known_physical_ground_O_fraction?: number; known_reactive_OH_H_fraction?: number
    unresolved_CO_fate_fraction?: number
    unresolved_HCN_reactive_fraction?: number
    hcn_ground_O_quench_fraction_of_all_collisions?: number
    lifetime_range_from_CO_rate_only_s?: [number, number]
    mean_first_collision_lifetime_s?: number; products_injected_into_solver: boolean
  } | null
  source_terms: PlasmaSourceTerm[]; warnings: string[]
}

export interface HeavySpeciesSnapshot {
  status: string; engine?: string; mechanism: string; mechanism_species?: number; mechanism_reactions?: number
  source_energy_used_J?: number; formed_mol: Record<string, number>; consumed_mol?: Record<string, number>
  collision_energy_allocation?: {
    all_collision_raw_threshold_demand_J: number
    overlapping_outcome_raw_threshold_J?: number
    resolved_fragment_raw_threshold_demand_J: number
    global_collision_energy_scale: number
    allocated_collision_energy_J: number
    resolved_fragment_energy_used_J: number
    nonfragment_or_inventory_limited_energy_J: number
    excluded_fragment_channels?: Array<{ channel_id: string; product_label: string; reason: string; raw_opportunity_mol: number }>
  }
  integrated_pathway_diagnostics?: {
    method: string; adaptive_interval_count: number; not_a_validated_yield_or_mechanism: boolean
    species: Record<string, {
      integrated_chemical_net_mol: number; observed_chemical_net_mol: number
      quadrature_closure_residual_mol: number
      formation: Array<{ reaction_index: number; equation: string; species_net_mol: number }>
      consumption: Array<{ reaction_index: number; equation: string; species_net_mol: number }>
    }>
  } | null
  mechanism_comparison?: {
    status: string; primary_mechanism: string
    shadow: {
      status: string; mechanism: string; mechanism_sha256: string
      mechanism_species?: number; retained_heavy_reactions?: number; excluded_plasma_rate_reactions?: number
      mapped_excitation_energy_used_J?: number
      mapped_excitation_source_extents?: Array<{ source_product_label: string; product: string; accepted_extent_mol: number }>
      formed_mol: Record<string, number>; consumed_mol?: Record<string, number>; net_change_mol?: Record<string, number>; warnings: string[]
    }
    charged_shadow: {
      status: string; mechanism: string; mechanism_sha256?: string; error?: string
      retained_reactions?: number; retained_two_temperature_reactions?: number
      disabled_embedded_electron_collision_reactions?: number
      source_charge_residual_mol_e?: number; final_charge_residual_mol_e?: number; max_atom_residual_mol?: number
      propagation_cycle?: number; carried_from_previous_cycle?: boolean; propagated_across_apparatus_cycles?: boolean
      branch_plasma_snapshot?: { status?: string; mole_fractions: Record<string, number>; covered_gas_fraction: number; mean_electron_energy_eV?: number; reactive_channel_count: number; balanced_reactive_channel_count: number }
      final_instantaneous_pathway_diagnostics?: Record<string, { formation: Array<{ reaction_index: number; equation: string; rate_type: string; species_net_rate_mol_s: number }>; consumption: Array<{ reaction_index: number; equation: string; rate_type: string; species_net_rate_mol_s: number }> }>
      source_extents?: Array<{ source_product_label: string; mapped_products: string[]; accepted_extent_mol: number; time_averaged_source_rate_mol_s: number }>
      unsupported_source_channels?: Array<{ channel_id: string; product_label: string; reason: string }>
      formed_mol: Record<string, number>; consumed_mol?: Record<string, number>; net_change_mol?: Record<string, number>; warnings: string[]
    }
    literature_hcn_shadow: {
      status: string; mechanism: string; mechanism_sha256: string
      published_species?: number; retained_heavy_reactions?: number; inert_bath_species_added?: number
      fixed_condition_electron_reactions_excluded?: number; error?: string
      requested_duration_s?: number; integrated_duration_s?: number; published_residence_time_s?: number; residence_time_domain_limited?: boolean
      formed_mol: Record<string, number>; consumed_mol?: Record<string, number>; net_change_mol?: Record<string, number>; warnings: string[]
    }
    species_comparison?: Record<string, {
      gri_net_change_mol: number; pavan_shadow_net_change_mol: number; pavan_to_gri_ratio: number | null
      pavan_charged_shadow_net_change_mol: number; pavan_charged_to_gri_ratio: number | null
      literature_hcn_shadow_net_change_mol: number; literature_hcn_to_gri_ratio: number | null
    }>
    scientific_firewall: {
      primary_inventory_mutated_by_shadow: boolean; comparison_is_independent_confidence: boolean
      quantitative_probability_claim_allowed?: boolean
    }
  }
  warnings: string[]
}

export interface AqueousCandidate {
  smiles: string; formula: string; motifs: string[]; heavy_atoms: number; ring_count: number
  stereochemistry: Array<{ atom_index: number; descriptor: string }>
}

export interface AqueousDiscovery {
  status: string
  aqueous_source_species: Array<{ id: string; amount_mol: number }>
  formula_count: number; candidate_structure_count: number
  post_discovery_motif_counts: Record<string, number>
  formula_records: Array<{ formula: string; composition: Record<string, number>;
    source_combinations: string[][]; candidate_count: number; enumeration_truncated: boolean;
    candidates: AqueousCandidate[] }>
  limitations: string[]
  scientific_firewall: { known_products_accessible_to_generator: boolean;
    reaction_templates_accessible_to_generator: boolean; historical_yields_accessible_to_generator: boolean;
    motifs_applied_after_generation: boolean; quantitative_yield_claim_allowed: boolean }
}

export interface PostDiscoveryValidation {
  status: string; discovery_artifact_sha256: string; frozen_artifact_bytes: number
  candidate_product: { name: string | null; smiles: string; present_in_frozen_discovery: boolean }
  reference_equilibrium: { intermediate_name: string; intermediate_smiles: string;
    estimated_equilibrium_mol_at_298_15_K: number; equilibrium_constant_M_minus_2: number;
    aqueous_volume_L: number; aqueous_volume_basis: string; source: string; secondary_value_source: string }
  hydrolysis_sensitivity: { duration_s: number;
    projection_method: string; cycles_replayed: number;
    scenarios: Record<string, { rate_s_inv: number; converted_fraction: number; candidate_product_mol: number;
      aminonitrile_remaining_mol: number; cumulative_aminonitrile_formed_mol: number;
      cycle_trace: Array<{ cycle: number; end_time_s: number; aminonitrile_remaining_mol: number;
        cumulative_aminonitrile_formed_mol: number; candidate_product_mol: number }> }>;
    evidence_temperature_K: number; evidence_ammonia_weight_percent: number;
    apparatus_liquid_temperature_range_K: number[]; apparatus_final_ammonia_weight_percent: number;
    source: string; applicability: string }
  limitations: string[]
  scientific_firewall: { discovery_frozen_before_evidence_load: boolean; evaluation_feedback_to_discovery: boolean;
    historical_benchmark_yields_used: boolean; known_product_identity_available_only_post_discovery: boolean;
    quantitative_yield_claim_allowed: boolean }
}

export interface HistoricalBenchmark {
  status: string; benchmark_id: string; discovery_artifact_sha256: string
  historical_conditions_matched: boolean; enumeration_coverage: number
  enumeration_coverage_is_prediction_score: boolean; validated_reproduction: boolean
  condition_checks: Record<string, { actual: number; target: number; unit: string; matched: boolean }>
  observation_comparison: Array<{ name: string; candidate_recovered: boolean; candidate_smiles: string[];
    captured_precursors_mol: Record<string, number>; precursor_supported: boolean; pathway_status: string;
    source: string }>
  remaining_blockers: string[]
  sources: { apparatus: string; protocol: string; mechanism: string }
  scientific_firewall: { benchmark_loaded_after_discovery_frozen: boolean;
    benchmark_observations_available_to_discovery: boolean; benchmark_yields_used_for_fitting: boolean;
    evaluation_feedback_to_simulation: boolean }
}

export interface PrecursorEnergyBenchmark {
  status: string; benchmark_id: string; comparable: boolean; composition_matched: boolean
  deposited_energy_J: number; energy_basis: string
  calorimetric_evidence: Record<string, unknown> | null
  pressure_self_consistent?: boolean | null
  declared_reference_pressure_Pa?: number | null
  simulated_gas_pressure_range_Pa?: [number, number] | null
  zero_time_extrapolation: Record<string, { available: boolean; covers_published_sampling_window: boolean
    sample_count: number; intercept_nmol_per_J: number | null;
    intercept_standard_error_nmol_per_J: number | null; slope_nmol_per_J_s: number | null;
    r_squared: number | null; method: string }>
  illustrative_sparse_sampling_sensitivity?: { status: string; requested_sample_times_s: number[];
    zero_time_extrapolation: PrecursorEnergyBenchmark['zero_time_extrapolation']; limitations: string }
  feed_checks: { required_species_present: Record<string, boolean>; excluded_species_absent: Record<string, boolean>
    H2_to_CH4_molar_ratio: number | null; H2_to_CH4_ratio_matched: boolean
    zero_time_yield_basis_matched: boolean; calorimetric_energy_basis_matched: boolean }
  species_comparison: Record<string, { simulated_nmol_per_J: number | null; reference_nmol_per_J: number
    observable_basis: string; final_inventory_diagnostic_nmol_per_J: number | null
    ratio_to_reference: number | null; scored: boolean }>
  remaining_blockers: string[]
  sources: { conference: string; metadata: string; journal: string }
}

export interface CandidateEnergyRanking {
  status: string; discovery_artifact_sha256: string; candidate_filter: string; candidate_count: number
  engines_requested: string[]; engines_succeeded: string[]
  engine_failures: Array<{ engine: string; error: string }>
  formula_local_consensus: Array<{ formula: string; smiles: string;
    stereochemistry: Array<{ atom_index: number; descriptor: string }>; engine_count: number;
    mean_formula_local_rank: number; rank_range: number[]; relative_energy_range_kJ_mol: number[] }>
  limitations: string[]
  scientific_firewall: { discovery_hash_verified_before_ranking: boolean; candidate_set_changed_by_ranking: boolean;
    historical_product_labels_available_to_ranker: boolean; energy_scales_mixed_between_engines: boolean;
    achiral_enantiomer_degeneracy_enforced: boolean; solution_free_energy_claimed: boolean;
    pathway_or_yield_claimed: boolean }
}

export interface CandidateHeldOutEvaluation {
  status: string; benchmark_id: string; scope: string; prediction_artifact_sha256: string
  prediction_artifact_bytes: number; discovery_artifact_sha256: string; top_k: number
  class_filter_coverage: number; top_k_formula_local_recall: number; validated_reproduction: boolean
  observation_comparison: Array<{ name: string; candidate_present_after_class_filter: boolean
    candidate_smiles: string[]; formula: string | null; best_mean_formula_local_rank: number | null
    rank_ranges: number[][]; top_k_recovered: boolean; source: string }>
  limitations: string[]
  scientific_firewall: { prediction_frozen_before_benchmark_read: boolean
    historical_product_labels_available_to_ranker: boolean; benchmark_feedback_to_prediction: boolean
    class_filter_selected_post_discovery: boolean; global_blind_product_prediction_claimed: boolean
    pathway_probability_or_yield_claimed: boolean }
}

export interface CandidateRankingJob {
  id: string; status: 'queued' | 'running' | 'completed' | 'failed'; stage: string; progress: number;
  message: string; result: CandidateEnergyRanking | null; evaluation: CandidateHeldOutEvaluation | null; error: string | null
}

export interface PrebioticResult {
  status: string; apparatus_id: string; cycles_completed: number; total_duration_s: number
  cycles: PrebioticCycleResult[]
  final_compartments_mol: Record<string, Record<string, number>>
  final_total_mol: Record<string, number>; reaction_extents_mol: Record<string, number>
  material_ledger: { closed: boolean; initial_mol: Record<string, number>; added_mol: Record<string, number>;
    sampled_mol: Record<string, number>; final_mol: Record<string, number>; residual: Record<string, number> }
  energy_ledger: { spark_mode: string; duty_fraction: number; deposited_electrical_energy_J: number;
    calorimetric_deposited_energy_J: number | null; deposited_energy_used_J: number;
    deposited_energy_basis: string;
    declared_channel_energy_J: number; heavy_species_source_energy_J: number;
    total_accounted_channel_energy_J: number; plasma_source_threshold_energy_lower_bound_J: number;
    status: string; reaction_was_energy_limited: boolean; plasma_source_was_energy_limited: boolean }
  evidence: { hypothetical_sources: string[]; all_rates_sourced: boolean }
  scientific_firewall: { known_products_accessible_to_discovery: boolean; reaction_recipe_used_for_discovery: boolean;
    blind_discovery_eligible: boolean; quantitative_yield_claim_allowed: boolean }
  stereochemistry_diagnostics: Array<{ species: string; enantiomer: string; difference_mol: number;
    racemic_within_numerical_tolerance: boolean; note: string }>
  aqueous_discovery: AqueousDiscovery
  post_discovery_validation: PostDiscoveryValidation
  historical_benchmark: HistoricalBenchmark
  precursor_energy_benchmark: PrecursorEnergyBenchmark
  hcn_formation: {
    status: string; numerical_reporting_floor_mol: number
    stoichiometric_atom_ceiling?: { status: string; available_atom_mol: Record<string, number>
      cyanide_family_total_ceiling_mol: number; molecular_hcn_total_ceiling_mol: number
      maximum_net_new_cyanide_family_mol: number; predicted_yield: boolean; interpretation: string }
    initial_hcn_mol: number; added_hcn_mol: number; sampled_hcn_mol: number
    final_hcn_by_compartment_mol: Record<string, number>
    final_hcn_total_mol: number; final_hcn_aqueous_mol: number
    final_hcn_plus_cyanide_aqueous_mol: number
    net_hcn_generated_mol: number; feed_carbon_mol_atoms: number
    net_hcn_per_feed_carbon: number | null
    final_hcn_total_nmol_per_declared_J: number | null
    final_hcn_aqueous_nmol_per_declared_J: number | null
    deposited_energy_basis: string; spark_mode: string
    material_ledger_closed: boolean
    experimental_detection_or_yield_validated: boolean
    interpretation: string
    aqueous_hydrolysis_screen?: {
      status: string; dataset_id: string; source: string
      combined_screen_mol: {
        nominal_mol: number; parameter_low_mol: number
        parameter_high_mol: number; abstract_acid_energy_mol: number
      }
      unavailable_exposure_count: number
      unavailable_exposures: Array<{cycle: number; substep?: number; compartment: string; reason: string}>
      prediction_inventory_mutated: boolean; kinetic_products_injected: boolean
      thermal_ledger_closed_by_this_screen: boolean; interpretation: string
    }
    charged_branch_aqueous_hydrolysis_screen?: {
      status: string; combined_screen_mol: {
        nominal_mol: number; parameter_low_mol: number; parameter_high_mol: number
      }; unavailable_exposure_count: number
      prediction_inventory_mutated: boolean
    }
  }
  operating_conditions: {
    operation_mode: 'staged' | 'continuous_overlap'; integration_substeps_per_cycle: number
    primary_source_forcing: 'impulse' | 'time_averaged'
    condensate_water_return_fraction: number; condensate_nonwater_return_fraction: number
    condenser_phase_transfer_model: 'step_equilibrium' | 'flow_limited_capture' | 'coupled_flow_transport'
    condenser_contact_water_holdup_L: number | null
    operator_splitting_convergence_assessed_in_this_run: boolean
    pressure_Pa: number; declared_reference_pressure_Pa: number
    simulated_gas_pressure_range_Pa: [number, number] | null; pressure_basis: string
    boiler_temperature_K: number; discharge_temperature_K: number; condenser_temperature_K: number
    boiler_liquid_volume_L: number; gas_loop_volume_L: number; condensate_trap_volume_L: number
    discharge_zone_volume_L: number; discharge_zone_fraction: number
    effective_heavy_chemistry_residence_time_s_per_cycle: number
    effective_heavy_chemistry_residence_time_s_per_integration_step: number
    discharge_zone_evidence: { value: number; unit: string; source: string; method: string; basis: string } | null
    effective_aqueous_volume_L: number; effective_aqueous_volume_basis: string; gas_circulation_rate_L_s: number
  }
  plasma_model: { mode: string; self_consistent_eedf_solved: boolean; apparatus_mixture_coupled: boolean;
    heavy_species_inventory_coupled: boolean; cycle_summaries: Array<ApparatusPlasmaSnapshot & { cycle: number }>;
    heavy_species_engine?: string | null; heavy_species_cycle_summaries: Array<HeavySpeciesSnapshot & { cycle: number }>;
    phase_partition_cycle_summaries: PhasePartitionCycle[];
    propagated_charged_branch_coupled?: boolean;
    propagated_charged_branch?: {
      status: string; cycles_completed: number
      final_compartments_mol: Record<string, Record<string, number>>
      internal_gas_state_mol: Record<string, number>
      material_ledger: { expected_conserved_mol: Record<string, number>; final_conserved_mol: Record<string, number>; residual: Record<string, number>; tolerance_mol: number; closed: boolean }
      precursor_energy_benchmark: PrecursorEnergyBenchmark
      simulated_gas_pressure_range_Pa?: [number, number] | null
      warnings: string[]
      scientific_firewall: { primary_inventory_mutated_by_shadow: boolean; comparison_is_independent_confidence: boolean; quantitative_yield_claim_allowed: boolean }
    } | null }
  warnings: string[]
  numerical_controls: { effective_aqueous_volume_L: number; effective_aqueous_volume_basis: string;
    [key: string]: number | string | boolean }
}

export interface RecordedPrebioticValidation {
  status: string
  runs: Array<{
    artifact_id: string; application_version: string; model_revision?: string; primary_source_forcing?: string
    condensate_nonwater_return_fraction?: number
    condenser_phase_transfer_model?: string; condenser_contact_water_holdup_L?: number | null
    cycles_completed: number; integration_steps: number; integration_substeps_per_cycle: number
    primary_material_ledger_closed: boolean; propagated_charged_material_ledger_closed: boolean
    quantitative_yield_claim_allowed: boolean
    primary: { hcn_zero_time_intercept_nmol_per_J: number; hcn_fit_r_squared: number;
      formaldehyde_zero_time_intercept_nmol_per_J: number; formaldehyde_fit_r_squared: number }
    propagated_charged: { hcn_zero_time_intercept_nmol_per_J: number; hcn_fit_r_squared: number;
      formaldehyde_zero_time_intercept_nmol_per_J: number; formaldehyde_fit_r_squared: number }
  }>
  convergence_studies: Array<{
    artifact_id: string; model_revision?: string; primary_source_forcing?: string; condensate_nonwater_return_fraction?: number
    condenser_phase_transfer_model?: string; condenser_contact_water_holdup_L?: number | null
    convergence_status: string
    results: Array<{ substeps: number; primary_water_hcn_mol: number; primary_water_formaldehyde_mol: number
      charged_water_hcn_mol: number; charged_water_formaldehyde_mol: number }>
  }>
  sensitivity_studies: Array<{
    artifact_id: string; status: string; scenario_parameter: string
    results: Array<{ condenser_contact_water_holdup_L: number;
      primary_water_hcn_mol: number; primary_water_formaldehyde_mol: number;
      charged_water_hcn_mol: number; charged_water_formaldehyde_mol: number }>
  }>
  scientific_firewall: { reference_yields_used_for_discovery: boolean;
    reference_yields_used_for_solver_fitting: boolean; quantitative_yield_claim_allowed: boolean }
}

export interface PrebioticCapabilities {
  apparatus: { available: boolean; compartments: string[]; features: string[] }
  reduced_order_discharge: { available: boolean; quantitative_yield_prediction: boolean; boundary: string }
  cantera_cross_section_plasma: { available: boolean; version: string | null; electron_collision_rate: boolean;
    reference_example_loaded: boolean; connected_to_apparatus: boolean; rate_coefficient_calculation: boolean;
    pavan_methane_plasma_loaded?: boolean; pavan_methane_plasma_species?: number | null;
    pavan_methane_plasma_reactions?: number | null; pavan_correlated_shadow_connected?: boolean;
    pavan_charge_balanced_shadow_connected?: boolean; pavan_retained_two_temperature_reactions?: number;
    pavan_disabled_embedded_electron_collision_reactions?: number;
    kamarinopoulou_hcn_shadow_connected?: boolean; kamarinopoulou_hcn_species?: number | null;
    kamarinopoulou_hcn_neutral_reactions?: number | null; kamarinopoulou_hcn_fixed_electron_reactions_excluded?: number;
    kamarinopoulou_hcn_mechanism_sha256?: string | null;
    eedf_modes: string[]; self_consistent_boltzmann_eedf: boolean; detail: string; required_inputs: string[] }
  boltzmann_eedf: { available: boolean; engine: string; version: string | null; targets: string[]
    opt_in_co2_o2_targets?: string[]
    mixture_composition_coupled: boolean; reduced_electric_field_coupled: boolean;
    quantitative_yield_prediction: boolean; collision_set: { available: boolean; dataset_id: string;
      sha256: string; source: string; detail: string; coverage?: CollisionCoverage };
    collision_sets?: Array<{ available: boolean; dataset_id: string; sha256: string; source: string;
      detail: string; coverage?: CollisionCoverage }> }
  blind_validation: { known_products_available_to_discovery: boolean; evaluation_data_connected: boolean }
}

export interface CollisionTargetCoverage {
  complete: boolean; completeness_scope: string; momentum_transfer_present: boolean; excitation_present: boolean
  neutral_dissociation_present?: boolean; ionization_present: boolean; attachment_present: boolean; process_count: number
  explicit_neutral_fragment_labels: string[]; unresolved_dissociation_labels: string[]
  explicit_neutral_fragment_channel_count: number; has_explicit_neutral_fragment_source: boolean
  resolved_ground_state_fragment_labels?: string[]; has_resolved_ground_state_fragment_source?: boolean
  unresolved_dissociation_channel_count: number
  channel_counts: Record<string, number>
}

export interface CollisionCoverage {
  all_required_targets_complete: boolean; completeness_scope: string
  any_unresolved_dissociation_labels: boolean; process_count: number
  all_required_targets_have_explicit_neutral_fragment_source: boolean
  all_required_targets_have_resolved_ground_state_fragment_source?: boolean
  targets: Record<string, CollisionTargetCoverage>; criteria: string[]
}

export interface CrossSectionValidationResult {
  status: string; dataset_id: string; sha256: string; cantera_version: string; cantera_rate_type: string
  points: number; energy_range_eV: [number, number]; first_nonzero_energy_eV: number | null
  maximum_cross_section_m2: number; integrated_cross_section_m2_eV: number
  scientific_firewall: { known_products_accessible_to_discovery: boolean; blind_discovery_eligible: boolean;
    usable_for_rate_without_eedf: boolean; quantitative_yield_claim_allowed: boolean }
  limitations: string[]
}

export interface CollisionRateResult {
  status: string
  eedf: { mode: string; mean_energy_eV: number; source: string; method: string; basis: string;
    cross_section_support_probability: number }
  electron_number_density_m3: number
  rate_coefficient_m3_s: number
  rate_coefficient_L_mol_s: number
  pseudo_first_order_rate_s_inv: number
  target_lifetime_s: number | null
  scientific_firewall: { known_products_accessible_to_discovery: boolean; rate_coefficient_calculated: boolean;
    cross_section_support_complete: boolean; complete_target_collision_set_declared: boolean;
    quantitative_yield_claim_allowed: boolean }
  limitations: string[]
}

export interface BoltzmannEEDFResult {
  status: string
  solver: { engine: string; version: string; method: string; iterations: number;
    final_relative_error: number; converged: boolean }
  conditions: { mole_fractions: Record<string, number>; gas_temperature_K: number;
    reduced_electric_field_Td: number; energy_max_eV: number; energy_cells: number; note: string }
  collision_set: { dataset_id: string; sha256: string; source: string; format: string;
    selected_process_count: number; eedf_process_count?: number; post_solve_outcome_process_count?: number;
    coverage: CollisionCoverage }
  eedf: { energy_levels_eV: number[]; distribution_eV_minus_3_over_2: number[];
    probability_density_eV_inv: number[]; normalization: number; mean_energy_eV: number;
    electron_temperature_K: number }
  rate_coefficients: Array<{ id: string; target: string; kind: string; product_label: string;
    threshold_eV: number; first_nonzero_energy_eV: number | null; rate_coefficient_m3_s: number;
    source?: string; source_sha256?: string; coupling?: string;
    outcome_only?: boolean; product_state_unresolved?: boolean }>
  independent_reference: { available_for_this_mixture: boolean; engine: string; reason?: string;
    unsupported_targets?: string[]; mechanism?: string; mean_energy_eV?: number;
    relative_mean_energy_difference?: number; comparison_scope?: string }
  same_collision_set_reference?: { available_for_this_mixture: boolean; engine: string; reason?: string
    collision_set_dataset_id?: string; collision_set_sha256?: string
    selected_eedf_process_count?: number; mean_energy_eV?: number
    relative_mean_energy_difference?: number; comparison_scope?: string
    rate_comparison?: { matched_process_count: number; material_rate_floor_m3_s: number
      material_process_count: number; material_within_factor_ten_count: number
      material_within_factor_ten_fraction: number | null
      largest_reactive_channel_by_target: Array<{ target: string; kind: string; product_label: string
        bolos_rate_coefficient_m3_s: number; cantera_eedf_rate_coefficient_m3_s: number
        factor_difference: number | null; outcome_only: boolean }>
      comparison_scope: string }
    private_collision_data_persisted?: boolean }
  scientific_firewall: { known_products_accessible_to_solver: boolean; reaction_recipe_used: boolean;
    self_consistent_eedf_solved: boolean; complete_starting_gas_collision_coverage: boolean;
    heavy_species_chemistry_complete: boolean; quantitative_yield_claim_allowed: boolean }
  numerical_warnings: string[]; limitations: string[]
}
