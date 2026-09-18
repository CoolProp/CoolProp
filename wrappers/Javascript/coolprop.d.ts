// Type definitions for the official CoolProp WebAssembly build.
//
// The package is an emscripten MODULARIZE + EXPORT_ES6 module: the default
// export is the module factory, and awaiting it yields the CoolProp API.
// Mirrors the embind bindings in src/emscripten_interface.cxx.

/** An embind enum entry (e.g. coolprop.parameters.iT). */
export interface EnumValue {
  readonly value: number;
}

/** An embind enum table: named entries, each carrying a numeric `.value`. */
export type EnumTable = Record<string, EnumValue>;

/** Phase-envelope vectors returned by AbstractState.get_phase_envelope_data(). */
export interface PhaseEnvelopeData {
  T: number[];
  p: number[];
  lnT: number[];
  lnp: number[];
  rhomolar_liq: number[];
  rhomolar_vap: number[];
  lnrhomolar_liq: number[];
  lnrhomolar_vap: number[];
  hmolar_liq: number[];
  hmolar_vap: number[];
  smolar_liq: number[];
  smolar_vap: number[];
  Q: number[];
  cpmolar_liq: number[];
  cpmolar_vap: number[];
  cvmolar_liq: number[];
  cvmolar_vap: number[];
  viscosity_liq: number[];
  viscosity_vap: number[];
  conductivity_liq: number[];
  conductivity_vap: number[];
  speed_sound_vap: number[];
}

/** Low-level state-based interface (see the CoolProp docs). */
export interface AbstractState {
  backend_name(): string;
  using_mole_fractions(): boolean;
  using_mass_fractions(): boolean;
  using_volu_fractions(): boolean;

  set_mole_fractions(fractions: number[]): void;
  set_mass_fractions(fractions: number[]): void;
  set_volu_fractions(fractions: number[]): void;
  get_mole_fractions(): number[];
  get_mass_fractions(): number[];
  mole_fractions_liquid(): number[];
  mole_fractions_vapor(): number[];

  update(pair: EnumValue, value1: number, value2: number): void;

  T(): number;
  rhomolar(): number;
  rhomass(): number;
  p(): number;
  Q(): number;
  tau(): number;
  delta(): number;
  molar_mass(): number;
  acentric_factor(): number;
  gas_constant(): number;
  Bvirial(): number;
  Cvirial(): number;
  compressibility_factor(): number;
  hmolar(): number;
  hmass(): number;
  smolar(): number;
  smass(): number;
  umolar(): number;
  umass(): number;
  cpmolar(): number;
  cpmass(): number;
  cvmolar(): number;
  cvmass(): number;
  gibbsmolar(): number;
  gibbsmass(): number;
  helmholtzmolar(): number;
  helmholtzmass(): number;
  speed_sound(): number;
  isothermal_compressibility(): number;
  isobaric_expansion_coefficient(): number;
  isentropic_expansion_coefficient(): number;
  viscosity(): number;
  conductivity(): number;
  surface_tension(): number;
  Prandtl(): number;

  keyed_output(key: EnumValue): number;
  trivial_keyed_output(key: EnumValue): number;
  saturated_liquid_keyed_output(key: EnumValue): number;
  saturated_vapor_keyed_output(key: EnumValue): number;

  first_partial_deriv(Of: EnumValue, Wrt: EnumValue, Constant: EnumValue): number;
  second_partial_deriv(Of: EnumValue, Wrt1: EnumValue, Constant1: EnumValue, Wrt2: EnumValue, Constant2: EnumValue): number;
  first_saturation_deriv(Of: EnumValue, Wrt: EnumValue): number;
  second_saturation_deriv(Of: EnumValue, Wrt1: EnumValue, Wrt2: EnumValue): number;
  first_two_phase_deriv(Of: EnumValue, Wrt: EnumValue, Constant: EnumValue): number;
  second_two_phase_deriv(Of: EnumValue, Wrt1: EnumValue, Constant1: EnumValue, Wrt2: EnumValue, Constant2: EnumValue): number;
  first_two_phase_deriv_splined(Of: EnumValue, Wrt: EnumValue, Constant: EnumValue, x_end: number): number;

  build_phase_envelope(level: string): void;
  get_phase_envelope_data(): PhaseEnvelopeData;

  melting_line(param: number, given: number, value: number): number;
  saturation_ancillary(param: EnumValue, Q: number, given: EnumValue, value: number): number;

  T_critical(): number;
  p_critical(): number;
  rhomolar_critical(): number;
  rhomass_critical(): number;
  T_reducing(): number;
  rhomolar_reducing(): number;
  rhomass_reducing(): number;
  p_triple(): number;
  Ttriple(): number;
  Tmin(): number;
  Tmax(): number;
  pmax(): number;
  dipole_moment(): number;

  /** Frees the underlying C++ object (embind). */
  delete(): void;
}

/** Options accepted by the emscripten module factory. */
export interface CoolPropModuleOptions {
  /** Custom resolution of the `.wasm` asset, e.g. `(f) => '/wasm/' + f`. */
  locateFile?: (path: string, prefix: string) => string;
}

/** The CoolProp API surface exposed by the initialized module. */
export interface CoolPropModule {
  /** Convert degrees Fahrenheit to Kelvin. */
  F2K(F: number): number;
  /** Simple fluid properties that do not depend on the state. */
  Props1SI(fluid: string, output: string): number;
  /** High-level state-based property evaluation. */
  PropsSI(output: string, name1: string, value1: number, name2: string, value2: number, fluid: string): number;
  /** Humid-air properties. */
  HAPropsSI(output: string, name1: string, value1: number, name2: string, value2: number, name3: string, value3: number): number;
  get_global_param_string(name: string): string;
  get_fluid_param_string(fluid: string, param: string): string;
  apply_simple_mixing_rule(identifier1: string, identifier2: string, rule: string): void;
  get_mixture_binary_pair_data(CAS1: string, CAS2: string, param: string): string;
  add_fluids_as_JSON(backend: string, json: string): boolean;

  /** Construct an AbstractState for `backend` and `&`-separated fluid names. */
  factory(backend: string, fluidNames: string): AbstractState;

  parameters: EnumTable;
  input_pairs: EnumTable;
  phases: EnumTable;
  backend_families: EnumTable;
}

declare function createCoolPropModule(options?: CoolPropModuleOptions): Promise<CoolPropModule>;

export default createCoolPropModule;
