from collections.abc import Sequence
import enum
from typing import Annotated, TypeAlias, overload

import numpy
from numpy.typing import NDArray


import numpy
from numpy.typing import NDArray
from collections.abc import Sequence
_Scalar = float
_Num = float | Sequence[float] | NDArray[numpy.float64]

class SimpleState:
    def __init__(self) -> None: ...

    @property
    def T(self) -> float: ...

    @T.setter
    def T(self, arg: float, /) -> None: ...

    @property
    def p(self) -> float: ...

    @p.setter
    def p(self, arg: float, /) -> None: ...

    @property
    def rhomolar(self) -> float: ...

    @rhomolar.setter
    def rhomolar(self, arg: float, /) -> None: ...

class GuessesStructure:
    def __init__(self) -> None: ...

    @property
    def T(self) -> float: ...

    @T.setter
    def T(self, arg: float, /) -> None: ...

    @property
    def p(self) -> float: ...

    @p.setter
    def p(self, arg: float, /) -> None: ...

    @property
    def rhomolar(self) -> float: ...

    @rhomolar.setter
    def rhomolar(self, arg: float, /) -> None: ...

    @property
    def hmolar(self) -> float: ...

    @hmolar.setter
    def hmolar(self, arg: float, /) -> None: ...

    @property
    def smolar(self) -> float: ...

    @smolar.setter
    def smolar(self, arg: float, /) -> None: ...

    @property
    def rhomolar_liq(self) -> float: ...

    @rhomolar_liq.setter
    def rhomolar_liq(self, arg: float, /) -> None: ...

    @property
    def rhomolar_vap(self) -> float: ...

    @rhomolar_vap.setter
    def rhomolar_vap(self, arg: float, /) -> None: ...

    @property
    def x(self) -> list[float]: ...

    @x.setter
    def x(self, arg: Sequence[float], /) -> None: ...

    @property
    def y(self) -> list[float]: ...

    @y.setter
    def y(self, arg: Sequence[float], /) -> None: ...

    def clear(self) -> None: ...

class CriticalState(SimpleState):
    def __init__(self) -> None: ...

    @property
    def stable(self) -> bool: ...

    @stable.setter
    def stable(self, arg: bool, /) -> None: ...

class PhaseEnvelopeData:
    def __init__(self) -> None: ...

    @property
    def K(self) -> list[list[float]]: ...

    @K.setter
    def K(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def lnK(self) -> list[list[float]]: ...

    @lnK.setter
    def lnK(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def x(self) -> list[list[float]]: ...

    @x.setter
    def x(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def y(self) -> list[list[float]]: ...

    @y.setter
    def y(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def T(self) -> list[float]: ...

    @T.setter
    def T(self, arg: Sequence[float], /) -> None: ...

    @property
    def p(self) -> list[float]: ...

    @p.setter
    def p(self, arg: Sequence[float], /) -> None: ...

    @property
    def lnT(self) -> list[float]: ...

    @lnT.setter
    def lnT(self, arg: Sequence[float], /) -> None: ...

    @property
    def lnp(self) -> list[float]: ...

    @lnp.setter
    def lnp(self, arg: Sequence[float], /) -> None: ...

    @property
    def rhomolar_liq(self) -> list[float]: ...

    @rhomolar_liq.setter
    def rhomolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def rhomolar_vap(self) -> list[float]: ...

    @rhomolar_vap.setter
    def rhomolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def lnrhomolar_liq(self) -> list[float]: ...

    @lnrhomolar_liq.setter
    def lnrhomolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def lnrhomolar_vap(self) -> list[float]: ...

    @lnrhomolar_vap.setter
    def lnrhomolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def hmolar_liq(self) -> list[float]: ...

    @hmolar_liq.setter
    def hmolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def hmolar_vap(self) -> list[float]: ...

    @hmolar_vap.setter
    def hmolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def smolar_liq(self) -> list[float]: ...

    @smolar_liq.setter
    def smolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def smolar_vap(self) -> list[float]: ...

    @smolar_vap.setter
    def smolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def Q(self) -> list[float]: ...

    @Q.setter
    def Q(self, arg: Sequence[float], /) -> None: ...

    @property
    def cpmolar_liq(self) -> list[float]: ...

    @cpmolar_liq.setter
    def cpmolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def cpmolar_vap(self) -> list[float]: ...

    @cpmolar_vap.setter
    def cpmolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def cvmolar_liq(self) -> list[float]: ...

    @cvmolar_liq.setter
    def cvmolar_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def cvmolar_vap(self) -> list[float]: ...

    @cvmolar_vap.setter
    def cvmolar_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def viscosity_liq(self) -> list[float]: ...

    @viscosity_liq.setter
    def viscosity_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def viscosity_vap(self) -> list[float]: ...

    @viscosity_vap.setter
    def viscosity_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def conductivity_liq(self) -> list[float]: ...

    @conductivity_liq.setter
    def conductivity_liq(self, arg: Sequence[float], /) -> None: ...

    @property
    def conductivity_vap(self) -> list[float]: ...

    @conductivity_vap.setter
    def conductivity_vap(self, arg: Sequence[float], /) -> None: ...

    @property
    def speed_sound_vap(self) -> list[float]: ...

    @speed_sound_vap.setter
    def speed_sound_vap(self, arg: Sequence[float], /) -> None: ...

class configuration_keys(enum.IntEnum):
    NORMALIZE_GAS_CONSTANTS = 0

    CRITICAL_WITHIN_1UK = 1

    CRITICAL_SPLINES_ENABLED = 2

    SAVE_RAW_TABLES = 3

    ALTERNATIVE_TABLES_DIRECTORY = 4

    ALTERNATIVE_SVDTABLES_DIRECTORY = 5

    ALTERNATIVE_REFPROP_PATH = 6

    ALTERNATIVE_REFPROP_HMX_BNC_PATH = 7

    ALTERNATIVE_REFPROP_LIBRARY_PATH = 8

    REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS = 9

    REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS = 10

    REFPROP_USE_GERG = 11

    REFPROP_ERROR_THRESHOLD = 12

    REFPROP_USE_PENGROBINSON = 13

    REFPROP_RESOLVE_COOLPROP_ALIASES = 14

    MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB = 15

    DONT_CHECK_PROPERTY_LIMITS = 16

    HENRYS_LAW_TO_GENERATE_VLE_GUESSES = 17

    PHASE_ENVELOPE_STARTING_PRESSURE_PA = 18

    R_U_CODATA = 19

    VTPR_UNIFAC_PATH = 20

    SPINODAL_MINIMUM_DELTA = 21

    OVERWRITE_FLUIDS = 22

    OVERWRITE_DEPARTURE_FUNCTION = 23

    OVERWRITE_BINARY_INTERACTION = 24

    USE_GUESSES_IN_PROPSSI = 25

    ASSUME_CRITICAL_POINT_STABLE = 26

    VTPR_ALWAYS_RELOAD_LIBRARY = 27

    FLOAT_PUNCTUATION = 28

    ENABLE_SUPERANCILLARIES = 29

    ENABLE_MELTING_CALORIC_HS = 30

    HSU_D_TWOPHASE_EOS_POLISH = 31

    LIST_STRING_DELIMITER = 32

    ALLOW_SVDSBTL_IN_PROPSSI = 33

    SVDSBTL_SAMPLING_THREADS = 34

    TABULAR_NX = 35

    MIXTURE_STABILITY_ALGORITHM = 36

    TABULAR_NY = 37

NORMALIZE_GAS_CONSTANTS: configuration_keys = configuration_keys.NORMALIZE_GAS_CONSTANTS

CRITICAL_WITHIN_1UK: configuration_keys = configuration_keys.CRITICAL_WITHIN_1UK

CRITICAL_SPLINES_ENABLED: configuration_keys = configuration_keys.CRITICAL_SPLINES_ENABLED

SAVE_RAW_TABLES: configuration_keys = configuration_keys.SAVE_RAW_TABLES

ALTERNATIVE_TABLES_DIRECTORY: configuration_keys = configuration_keys.ALTERNATIVE_TABLES_DIRECTORY

ALTERNATIVE_SVDTABLES_DIRECTORY: configuration_keys = ...

ALTERNATIVE_REFPROP_PATH: configuration_keys = configuration_keys.ALTERNATIVE_REFPROP_PATH

ALTERNATIVE_REFPROP_HMX_BNC_PATH: configuration_keys = ...

ALTERNATIVE_REFPROP_LIBRARY_PATH: configuration_keys = ...

REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS: configuration_keys = ...

REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS: configuration_keys = ...

REFPROP_USE_GERG: configuration_keys = configuration_keys.REFPROP_USE_GERG

REFPROP_ERROR_THRESHOLD: configuration_keys = configuration_keys.REFPROP_ERROR_THRESHOLD

REFPROP_USE_PENGROBINSON: configuration_keys = configuration_keys.REFPROP_USE_PENGROBINSON

REFPROP_RESOLVE_COOLPROP_ALIASES: configuration_keys = ...

MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB: configuration_keys = ...

DONT_CHECK_PROPERTY_LIMITS: configuration_keys = configuration_keys.DONT_CHECK_PROPERTY_LIMITS

HENRYS_LAW_TO_GENERATE_VLE_GUESSES: configuration_keys = ...

PHASE_ENVELOPE_STARTING_PRESSURE_PA: configuration_keys = ...

R_U_CODATA: configuration_keys = configuration_keys.R_U_CODATA

VTPR_UNIFAC_PATH: configuration_keys = configuration_keys.VTPR_UNIFAC_PATH

SPINODAL_MINIMUM_DELTA: configuration_keys = configuration_keys.SPINODAL_MINIMUM_DELTA

OVERWRITE_FLUIDS: configuration_keys = configuration_keys.OVERWRITE_FLUIDS

OVERWRITE_DEPARTURE_FUNCTION: configuration_keys = configuration_keys.OVERWRITE_DEPARTURE_FUNCTION

OVERWRITE_BINARY_INTERACTION: configuration_keys = configuration_keys.OVERWRITE_BINARY_INTERACTION

USE_GUESSES_IN_PROPSSI: configuration_keys = configuration_keys.USE_GUESSES_IN_PROPSSI

ASSUME_CRITICAL_POINT_STABLE: configuration_keys = configuration_keys.ASSUME_CRITICAL_POINT_STABLE

VTPR_ALWAYS_RELOAD_LIBRARY: configuration_keys = configuration_keys.VTPR_ALWAYS_RELOAD_LIBRARY

FLOAT_PUNCTUATION: configuration_keys = configuration_keys.FLOAT_PUNCTUATION

ENABLE_SUPERANCILLARIES: configuration_keys = configuration_keys.ENABLE_SUPERANCILLARIES

ENABLE_MELTING_CALORIC_HS: configuration_keys = configuration_keys.ENABLE_MELTING_CALORIC_HS

HSU_D_TWOPHASE_EOS_POLISH: configuration_keys = configuration_keys.HSU_D_TWOPHASE_EOS_POLISH

LIST_STRING_DELIMITER: configuration_keys = configuration_keys.LIST_STRING_DELIMITER

ALLOW_SVDSBTL_IN_PROPSSI: configuration_keys = configuration_keys.ALLOW_SVDSBTL_IN_PROPSSI

SVDSBTL_SAMPLING_THREADS: configuration_keys = configuration_keys.SVDSBTL_SAMPLING_THREADS

TABULAR_NX: configuration_keys = configuration_keys.TABULAR_NX

MIXTURE_STABILITY_ALGORITHM: configuration_keys = configuration_keys.MIXTURE_STABILITY_ALGORITHM

TABULAR_NY: configuration_keys = configuration_keys.TABULAR_NY

class parameters(enum.IntEnum):
    igas_constant = 1

    imolar_mass = 2

    iacentric_factor = 3

    irhomolar_reducing = 4

    irhomolar_critical = 5

    iT_reducing = 6

    iT_critical = 7

    irhomass_reducing = 8

    irhomass_critical = 9

    iP_critical = 10

    iP_reducing = 11

    iT_triple = 12

    iP_triple = 13

    iT_min = 14

    iT_max = 15

    iP_max = 16

    iP_min = 17

    idipole_moment = 18

    iT = 19

    iP = 20

    iQ = 21

    iTau = 23

    iDelta = 24

    iDmolar = 25

    iHmolar = 26

    iSmolar = 27

    iCpmolar = 28

    iCp0molar = 29

    iCvmolar = 30

    iUmolar = 31

    iGmolar = 32

    iHelmholtzmolar = 33

    iSmolar_residual = 35

    iHmolar_residual = 34

    iGmolar_residual = 36

    iDmass = 40

    iHmass = 41

    iSmass = 42

    iCpmass = 43

    iCp0mass = 44

    iCvmass = 45

    iUmass = 46

    iGmass = 47

    iHelmholtzmass = 48

    iviscosity = 52

    iconductivity = 53

    isurface_tension = 54

    iPrandtl = 55

    ispeed_sound = 56

    iisothermal_compressibility = 57

    iisobaric_expansion_coefficient = 58

    ifundamental_derivative_of_gas_dynamics = 60

    ialphar = 61

    idalphar_ddelta_consttau = 63

    idalpha0_dtau_constdelta = 65

    iBvirial = 69

    iCvirial = 70

    idBvirial_dT = 71

    idCvirial_dT = 72

    iZ = 73

    iPIP = 74

    ifraction_min = 75

    ifraction_max = 76

    iT_freeze = 77

    iGWP20 = 78

    iGWP100 = 79

    iGWP500 = 80

    iFH = 81

    iHH = 82

    iPH = 83

    iODP = 84

    iPhase = 85

    INVALID_PARAMETER = 0

    iQmass = 22

    iHmolar_idealgas = 37

    iSmolar_idealgas = 38

    iUmolar_idealgas = 39

    iHmass_idealgas = 49

    iSmass_idealgas = 50

    iUmass_idealgas = 51

    iisentropic_expansion_coefficient = 59

    idalphar_dtau_constdelta = 62

    ialpha0 = 64

    idalpha0_ddelta_consttau = 66

    id2alpha0_ddelta2_consttau = 67

    id3alpha0_ddelta3_consttau = 68

    iundefined_parameter = 86

igas_constant: parameters = parameters.igas_constant

imolar_mass: parameters = parameters.imolar_mass

iacentric_factor: parameters = parameters.iacentric_factor

irhomolar_reducing: parameters = parameters.irhomolar_reducing

irhomolar_critical: parameters = parameters.irhomolar_critical

iT_reducing: parameters = parameters.iT_reducing

iT_critical: parameters = parameters.iT_critical

irhomass_reducing: parameters = parameters.irhomass_reducing

irhomass_critical: parameters = parameters.irhomass_critical

iP_critical: parameters = parameters.iP_critical

iP_reducing: parameters = parameters.iP_reducing

iT_triple: parameters = parameters.iT_triple

iP_triple: parameters = parameters.iP_triple

iT_min: parameters = parameters.iT_min

iT_max: parameters = parameters.iT_max

iP_max: parameters = parameters.iP_max

iP_min: parameters = parameters.iP_min

idipole_moment: parameters = parameters.idipole_moment

iT: parameters = parameters.iT

iP: parameters = parameters.iP

iQ: parameters = parameters.iQ

iTau: parameters = parameters.iTau

iDelta: parameters = parameters.iDelta

iDmolar: parameters = parameters.iDmolar

iHmolar: parameters = parameters.iHmolar

iSmolar: parameters = parameters.iSmolar

iCpmolar: parameters = parameters.iCpmolar

iCp0molar: parameters = parameters.iCp0molar

iCvmolar: parameters = parameters.iCvmolar

iUmolar: parameters = parameters.iUmolar

iGmolar: parameters = parameters.iGmolar

iHelmholtzmolar: parameters = parameters.iHelmholtzmolar

iSmolar_residual: parameters = parameters.iSmolar_residual

iHmolar_residual: parameters = parameters.iHmolar_residual

iGmolar_residual: parameters = parameters.iGmolar_residual

iDmass: parameters = parameters.iDmass

iHmass: parameters = parameters.iHmass

iSmass: parameters = parameters.iSmass

iCpmass: parameters = parameters.iCpmass

iCp0mass: parameters = parameters.iCp0mass

iCvmass: parameters = parameters.iCvmass

iUmass: parameters = parameters.iUmass

iGmass: parameters = parameters.iGmass

iHelmholtzmass: parameters = parameters.iHelmholtzmass

iviscosity: parameters = parameters.iviscosity

iconductivity: parameters = parameters.iconductivity

isurface_tension: parameters = parameters.isurface_tension

iPrandtl: parameters = parameters.iPrandtl

ispeed_sound: parameters = parameters.ispeed_sound

iisothermal_compressibility: parameters = parameters.iisothermal_compressibility

iisobaric_expansion_coefficient: parameters = parameters.iisobaric_expansion_coefficient

ifundamental_derivative_of_gas_dynamics: parameters = ...

ialphar: parameters = parameters.ialphar

idalphar_ddelta_consttau: parameters = parameters.idalphar_ddelta_consttau

idalpha0_dtau_constdelta: parameters = parameters.idalpha0_dtau_constdelta

iBvirial: parameters = parameters.iBvirial

iCvirial: parameters = parameters.iCvirial

idBvirial_dT: parameters = parameters.idBvirial_dT

idCvirial_dT: parameters = parameters.idCvirial_dT

iZ: parameters = parameters.iZ

iPIP: parameters = parameters.iPIP

ifraction_min: parameters = parameters.ifraction_min

ifraction_max: parameters = parameters.ifraction_max

iT_freeze: parameters = parameters.iT_freeze

iGWP20: parameters = parameters.iGWP20

iGWP100: parameters = parameters.iGWP100

iGWP500: parameters = parameters.iGWP500

iFH: parameters = parameters.iFH

iHH: parameters = parameters.iHH

iPH: parameters = parameters.iPH

iODP: parameters = parameters.iODP

iPhase: parameters = parameters.iPhase

INVALID_PARAMETER: parameters = parameters.INVALID_PARAMETER

iQmass: parameters = parameters.iQmass

iHmolar_idealgas: parameters = parameters.iHmolar_idealgas

iSmolar_idealgas: parameters = parameters.iSmolar_idealgas

iUmolar_idealgas: parameters = parameters.iUmolar_idealgas

iHmass_idealgas: parameters = parameters.iHmass_idealgas

iSmass_idealgas: parameters = parameters.iSmass_idealgas

iUmass_idealgas: parameters = parameters.iUmass_idealgas

iisentropic_expansion_coefficient: parameters = parameters.iisentropic_expansion_coefficient

idalphar_dtau_constdelta: parameters = parameters.idalphar_dtau_constdelta

ialpha0: parameters = parameters.ialpha0

idalpha0_ddelta_consttau: parameters = parameters.idalpha0_ddelta_consttau

id2alpha0_ddelta2_consttau: parameters = parameters.id2alpha0_ddelta2_consttau

id3alpha0_ddelta3_consttau: parameters = parameters.id3alpha0_ddelta3_consttau

iundefined_parameter: parameters = parameters.iundefined_parameter

class input_pairs(enum.IntEnum):
    QT_INPUTS = 1

    PQ_INPUTS = 3

    QSmolar_INPUTS = 5

    QSmass_INPUTS = 7

    HmolarQ_INPUTS = 9

    HmassQ_INPUTS = 11

    DmolarQ_INPUTS = 13

    DmassQ_INPUTS = 15

    PT_INPUTS = 17

    DmassT_INPUTS = 18

    DmolarT_INPUTS = 19

    HmolarT_INPUTS = 20

    HmassT_INPUTS = 21

    SmolarT_INPUTS = 22

    SmassT_INPUTS = 23

    TUmolar_INPUTS = 24

    TUmass_INPUTS = 25

    DmassP_INPUTS = 26

    DmolarP_INPUTS = 27

    HmassP_INPUTS = 28

    HmolarP_INPUTS = 29

    PSmass_INPUTS = 30

    PSmolar_INPUTS = 31

    PUmass_INPUTS = 32

    PUmolar_INPUTS = 33

    HmassSmass_INPUTS = 34

    HmolarSmolar_INPUTS = 35

    SmassUmass_INPUTS = 36

    SmolarUmolar_INPUTS = 37

    DmassHmass_INPUTS = 38

    DmolarHmolar_INPUTS = 39

    DmassSmass_INPUTS = 40

    DmolarSmolar_INPUTS = 41

    DmassUmass_INPUTS = 42

    DmolarUmolar_INPUTS = 43

    INPUT_PAIR_INVALID = 0

    QmassT_INPUTS = 2

    PQmass_INPUTS = 4

    QmassSmolar_INPUTS = 6

    QmassSmass_INPUTS = 8

    HmolarQmass_INPUTS = 10

    HmassQmass_INPUTS = 12

    DmolarQmass_INPUTS = 14

    DmassQmass_INPUTS = 16

QT_INPUTS: input_pairs = input_pairs.QT_INPUTS

PQ_INPUTS: input_pairs = input_pairs.PQ_INPUTS

QSmolar_INPUTS: input_pairs = input_pairs.QSmolar_INPUTS

QSmass_INPUTS: input_pairs = input_pairs.QSmass_INPUTS

HmolarQ_INPUTS: input_pairs = input_pairs.HmolarQ_INPUTS

HmassQ_INPUTS: input_pairs = input_pairs.HmassQ_INPUTS

DmolarQ_INPUTS: input_pairs = input_pairs.DmolarQ_INPUTS

DmassQ_INPUTS: input_pairs = input_pairs.DmassQ_INPUTS

PT_INPUTS: input_pairs = input_pairs.PT_INPUTS

DmassT_INPUTS: input_pairs = input_pairs.DmassT_INPUTS

DmolarT_INPUTS: input_pairs = input_pairs.DmolarT_INPUTS

HmolarT_INPUTS: input_pairs = input_pairs.HmolarT_INPUTS

HmassT_INPUTS: input_pairs = input_pairs.HmassT_INPUTS

SmolarT_INPUTS: input_pairs = input_pairs.SmolarT_INPUTS

SmassT_INPUTS: input_pairs = input_pairs.SmassT_INPUTS

TUmolar_INPUTS: input_pairs = input_pairs.TUmolar_INPUTS

TUmass_INPUTS: input_pairs = input_pairs.TUmass_INPUTS

DmassP_INPUTS: input_pairs = input_pairs.DmassP_INPUTS

DmolarP_INPUTS: input_pairs = input_pairs.DmolarP_INPUTS

HmassP_INPUTS: input_pairs = input_pairs.HmassP_INPUTS

HmolarP_INPUTS: input_pairs = input_pairs.HmolarP_INPUTS

PSmass_INPUTS: input_pairs = input_pairs.PSmass_INPUTS

PSmolar_INPUTS: input_pairs = input_pairs.PSmolar_INPUTS

PUmass_INPUTS: input_pairs = input_pairs.PUmass_INPUTS

PUmolar_INPUTS: input_pairs = input_pairs.PUmolar_INPUTS

HmassSmass_INPUTS: input_pairs = input_pairs.HmassSmass_INPUTS

HmolarSmolar_INPUTS: input_pairs = input_pairs.HmolarSmolar_INPUTS

SmassUmass_INPUTS: input_pairs = input_pairs.SmassUmass_INPUTS

SmolarUmolar_INPUTS: input_pairs = input_pairs.SmolarUmolar_INPUTS

DmassHmass_INPUTS: input_pairs = input_pairs.DmassHmass_INPUTS

DmolarHmolar_INPUTS: input_pairs = input_pairs.DmolarHmolar_INPUTS

DmassSmass_INPUTS: input_pairs = input_pairs.DmassSmass_INPUTS

DmolarSmolar_INPUTS: input_pairs = input_pairs.DmolarSmolar_INPUTS

DmassUmass_INPUTS: input_pairs = input_pairs.DmassUmass_INPUTS

DmolarUmolar_INPUTS: input_pairs = input_pairs.DmolarUmolar_INPUTS

INPUT_PAIR_INVALID: input_pairs = input_pairs.INPUT_PAIR_INVALID

QmassT_INPUTS: input_pairs = input_pairs.QmassT_INPUTS

PQmass_INPUTS: input_pairs = input_pairs.PQmass_INPUTS

QmassSmolar_INPUTS: input_pairs = input_pairs.QmassSmolar_INPUTS

QmassSmass_INPUTS: input_pairs = input_pairs.QmassSmass_INPUTS

HmolarQmass_INPUTS: input_pairs = input_pairs.HmolarQmass_INPUTS

HmassQmass_INPUTS: input_pairs = input_pairs.HmassQmass_INPUTS

DmolarQmass_INPUTS: input_pairs = input_pairs.DmolarQmass_INPUTS

DmassQmass_INPUTS: input_pairs = input_pairs.DmassQmass_INPUTS

class phases(enum.IntEnum):
    iphase_liquid = 0

    iphase_supercritical = 1

    iphase_supercritical_gas = 2

    iphase_supercritical_liquid = 3

    iphase_critical_point = 4

    iphase_gas = 5

    iphase_twophase = 6

    iphase_unknown = 7

    iphase_not_imposed = 8

iphase_liquid: phases = phases.iphase_liquid

iphase_supercritical: phases = phases.iphase_supercritical

iphase_supercritical_gas: phases = phases.iphase_supercritical_gas

iphase_supercritical_liquid: phases = phases.iphase_supercritical_liquid

iphase_critical_point: phases = phases.iphase_critical_point

iphase_gas: phases = phases.iphase_gas

iphase_twophase: phases = phases.iphase_twophase

iphase_unknown: phases = phases.iphase_unknown

iphase_not_imposed: phases = phases.iphase_not_imposed

class fluid_types(enum.IntEnum):
    FLUID_TYPE_PURE = 0

    FLUID_TYPE_PSEUDOPURE = 1

    FLUID_TYPE_REFPROP = 2

    FLUID_TYPE_INCOMPRESSIBLE_LIQUID = 3

    FLUID_TYPE_INCOMPRESSIBLE_SOLUTION = 4

    FLUID_TYPE_UNDEFINED = 5

FLUID_TYPE_PURE: fluid_types = fluid_types.FLUID_TYPE_PURE

FLUID_TYPE_PSEUDOPURE: fluid_types = fluid_types.FLUID_TYPE_PSEUDOPURE

FLUID_TYPE_REFPROP: fluid_types = fluid_types.FLUID_TYPE_REFPROP

FLUID_TYPE_INCOMPRESSIBLE_LIQUID: fluid_types = fluid_types.FLUID_TYPE_INCOMPRESSIBLE_LIQUID

FLUID_TYPE_INCOMPRESSIBLE_SOLUTION: fluid_types = fluid_types.FLUID_TYPE_INCOMPRESSIBLE_SOLUTION

FLUID_TYPE_UNDEFINED: fluid_types = fluid_types.FLUID_TYPE_UNDEFINED

class fast_evaluate_status(enum.IntEnum):
    fast_evaluate_ok = 0

    fast_evaluate_out_of_range = 1

    fast_evaluate_two_phase_disallowed = 2

    fast_evaluate_unsupported_input = 3

    fast_evaluate_unsupported_output = 4

    fast_evaluate_internal_error = 5

fast_evaluate_ok: fast_evaluate_status = fast_evaluate_status.fast_evaluate_ok

fast_evaluate_out_of_range: fast_evaluate_status = fast_evaluate_status.fast_evaluate_out_of_range

fast_evaluate_two_phase_disallowed: fast_evaluate_status = ...

fast_evaluate_unsupported_input: fast_evaluate_status = ...

fast_evaluate_unsupported_output: fast_evaluate_status = ...

fast_evaluate_internal_error: fast_evaluate_status = fast_evaluate_status.fast_evaluate_internal_error

class SpinodalData:
    def __init__(self) -> None: ...

    @property
    def tau(self) -> list[float]: ...

    @property
    def delta(self) -> list[float]: ...

    @property
    def M1(self) -> list[float]: ...

PyCriticalState: TypeAlias = CriticalState

PyGuessesStructure: TypeAlias = GuessesStructure

PyPhaseEnvelopeData: TypeAlias = PhaseEnvelopeData

PySpinodalData: TypeAlias = SpinodalData

class AbstractState:
    def __init__(self, arg0: str, arg1: str, /) -> None: ...

    def set_T(self, arg: float, /) -> None: ...

    def backend_name(self) -> str: ...

    def using_mole_fractions(self) -> bool: ...

    def using_mass_fractions(self) -> bool: ...

    def using_volu_fractions(self) -> bool: ...

    def set_mole_fractions(self, arg: Sequence[float], /) -> None: ...

    def set_mass_fractions(self, arg: Sequence[float], /) -> None: ...

    def set_volu_fractions(self, arg: Sequence[float], /) -> None: ...

    def mole_fractions_liquid(self) -> list[float]: ...

    def mole_fractions_liquid_double(self) -> list[float]: ...

    def mole_fractions_vapor(self) -> list[float]: ...

    def mole_fractions_vapor_double(self) -> list[float]: ...

    def get_mole_fractions(self) -> list[float]: ...

    def get_mass_fractions(self) -> list[float]: ...

    def update(self, arg0: input_pairs, arg1: float, arg2: float, /) -> None: ...

    def update_with_guesses(self, arg0: input_pairs, arg1: float, arg2: float, arg3: GuessesStructure, /) -> None: ...

    def available_in_high_level(self) -> bool: ...

    def build_options_json(self) -> str: ...

    def fluid_param_string(self, arg: str, /) -> str: ...

    def fluid_names(self) -> list[str]: ...

    @overload
    def set_binary_interaction_double(self, arg0: str, arg1: str, arg2: str, arg3: float, /) -> None: ...

    @overload
    def set_binary_interaction_double(self, arg0: int, arg1: int, arg2: str, arg3: float, /) -> None: ...

    @overload
    def set_binary_interaction_string(self, arg0: str, arg1: str, arg2: str, arg3: str, /) -> None: ...

    @overload
    def set_binary_interaction_string(self, arg0: int, arg1: int, arg2: str, arg3: str, /) -> None: ...

    @overload
    def get_binary_interaction_double(self, arg0: str, arg1: str, arg2: str, /) -> float: ...

    @overload
    def get_binary_interaction_double(self, arg0: int, arg1: int, arg2: str, /) -> float: ...

    def get_binary_interaction_string(self, arg0: str, arg1: str, arg2: str, /) -> str: ...

    def apply_simple_mixing_rule(self, arg0: int, arg1: int, arg2: str, /) -> None: ...

    def set_fluid_parameter_double(self, arg0: int, arg1: str, arg2: float, /) -> None: ...

    def clear(self) -> bool: ...

    def get_reducing_state(self) -> SimpleState: ...

    def get_state(self, arg: str, /) -> SimpleState: ...

    def Tmin(self) -> float: ...

    def Tmax(self) -> float: ...

    def pmax(self) -> float: ...

    def Ttriple(self) -> float: ...

    def phase(self) -> phases: ...

    def specify_phase(self, arg: phases, /) -> None: ...

    def unspecify_phase(self) -> None: ...

    def T_critical(self) -> float: ...

    def p_critical(self) -> float: ...

    def rhomolar_critical(self) -> float: ...

    def rhomass_critical(self) -> float: ...

    def all_critical_points(self) -> list[CriticalState]: ...

    def build_spinodal(self) -> None: ...

    def get_spinodal_data(self) -> object: ...

    def criticality_contour_values(self) -> tuple: ...

    def tangent_plane_distance(self, arg0: float, arg1: float, arg2: Sequence[float], arg3: float, /) -> float: ...

    def T_reducing(self) -> float: ...

    def rhomolar_reducing(self) -> float: ...

    def rhomass_reducing(self) -> float: ...

    def p_triple(self) -> float: ...

    def name(self) -> str: ...

    def dipole_moment(self) -> float: ...

    def keyed_output(self, arg: parameters, /) -> float: ...

    def trivial_keyed_output(self, arg: parameters, /) -> float: ...

    def fast_evaluate(self, input_pair: input_pairs, val1: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], val2: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], outputs: Annotated[NDArray[numpy.int32], dict(shape=(None,), order='C', writable=False)], out: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C')], status: Annotated[NDArray[numpy.int32], dict(shape=(None,), order='C')], imposed_phase: phases = phases.iphase_not_imposed) -> None: ...

    def saturated_liquid_keyed_output(self, arg: parameters, /) -> float: ...

    def saturated_vapor_keyed_output(self, arg: parameters, /) -> float: ...

    def T(self) -> float: ...

    def rhomolar(self) -> float: ...

    def rhomass(self) -> float: ...

    def p(self) -> float: ...

    def Q(self) -> float: ...

    def Qmass(self) -> float: ...

    def tau(self) -> float: ...

    def delta(self) -> float: ...

    def molar_mass(self) -> float: ...

    def acentric_factor(self) -> float: ...

    def gas_constant(self) -> float: ...

    def Bvirial(self) -> float: ...

    def dBvirial_dT(self) -> float: ...

    def Cvirial(self) -> float: ...

    def dCvirial_dT(self) -> float: ...

    def compressibility_factor(self) -> float: ...

    def hmolar(self) -> float: ...

    def hmass(self) -> float: ...

    def hmolar_excess(self) -> float: ...

    def hmass_excess(self) -> float: ...

    def smolar(self) -> float: ...

    def smass(self) -> float: ...

    def smolar_excess(self) -> float: ...

    def smass_excess(self) -> float: ...

    def umolar(self) -> float: ...

    def umass(self) -> float: ...

    def umolar_excess(self) -> float: ...

    def umass_excess(self) -> float: ...

    def cpmolar(self) -> float: ...

    def cpmass(self) -> float: ...

    def cp0molar(self) -> float: ...

    def cp0mass(self) -> float: ...

    def cvmolar(self) -> float: ...

    def cvmass(self) -> float: ...

    def gibbsmolar(self) -> float: ...

    def gibbsmass(self) -> float: ...

    def gibbsmolar_excess(self) -> float: ...

    def gibbsmass_excess(self) -> float: ...

    def helmholtzmolar(self) -> float: ...

    def helmholtzmass(self) -> float: ...

    def helmholtzmolar_excess(self) -> float: ...

    def helmholtzmass_excess(self) -> float: ...

    def volumemolar_excess(self) -> float: ...

    def volumemass_excess(self) -> float: ...

    def hmolar_idealgas(self) -> float: ...

    def hmass_idealgas(self) -> float: ...

    def hmolar_residual(self) -> float: ...

    def smolar_idealgas(self) -> float: ...

    def smass_idealgas(self) -> float: ...

    def smolar_residual(self) -> float: ...

    def umolar_idealgas(self) -> float: ...

    def umass_idealgas(self) -> float: ...

    def gibbsmolar_residual(self) -> float: ...

    def neff(self) -> float: ...

    def get_fluid_constant(self, arg0: int, arg1: parameters, /) -> float: ...

    def get_fluid_parameter_double(self, arg0: int, arg1: str, /) -> float: ...

    def set_cubic_alpha_C(self, arg0: int, arg1: str, arg2: float, arg3: float, arg4: float, /) -> None: ...

    def update_QT_pure_superanc(self, arg0: float, arg1: float, /) -> None: ...

    def speed_sound(self) -> float: ...

    def isothermal_compressibility(self) -> float: ...

    def isobaric_expansion_coefficient(self) -> float: ...

    def fugacity_coefficient(self, arg: int, /) -> float: ...

    def fugacity(self, arg: int, /) -> float: ...

    def chemical_potential(self, arg: int, /) -> float: ...

    def fundamental_derivative_of_gas_dynamics(self) -> float: ...

    def PIP(self) -> float: ...

    def true_critical_point(self) -> tuple: ...

    def ideal_curve(self, arg: str, /) -> tuple: ...

    def first_partial_deriv(self, arg0: parameters, arg1: parameters, arg2: parameters, /) -> float: ...

    def second_partial_deriv(self, arg0: parameters, arg1: parameters, arg2: parameters, arg3: parameters, arg4: parameters, /) -> float: ...

    def first_saturation_deriv(self, arg0: parameters, arg1: parameters, /) -> float: ...

    def second_saturation_deriv(self, arg0: parameters, arg1: parameters, arg2: parameters, /) -> float: ...

    def first_two_phase_deriv(self, arg0: parameters, arg1: parameters, arg2: parameters, /) -> float: ...

    def second_two_phase_deriv(self, arg0: parameters, arg1: parameters, arg2: parameters, arg3: parameters, arg4: parameters, /) -> float: ...

    def first_two_phase_deriv_splined(self, arg0: parameters, arg1: parameters, arg2: parameters, arg3: float, /) -> float: ...

    def build_phase_envelope(self, arg: str, /) -> None: ...

    def get_phase_envelope_data(self) -> PhaseEnvelopeData: ...

    def has_melting_line(self) -> bool: ...

    def melting_line(self, arg0: int, arg1: int, arg2: float, /) -> float: ...

    def saturation_ancillary(self, arg0: parameters, arg1: int, arg2: parameters, arg3: float, /) -> float: ...

    def viscosity(self) -> float: ...

    def viscosity_contributions(self) -> dict: ...

    def conductivity(self) -> float: ...

    def conductivity_contributions(self) -> dict: ...

    def surface_tension(self) -> float: ...

    def Prandtl(self) -> float: ...

    def conformal_state(self, arg0: str, arg1: float, arg2: float, /) -> dict: ...

    def change_EOS(self, arg0: int, arg1: str, /) -> None: ...

    def alpha0(self) -> float: ...

    def dalpha0_dDelta(self) -> float: ...

    def dalpha0_dTau(self) -> float: ...

    def d2alpha0_dDelta2(self) -> float: ...

    def d2alpha0_dDelta_dTau(self) -> float: ...

    def d2alpha0_dTau2(self) -> float: ...

    def d3alpha0_dTau3(self) -> float: ...

    def d3alpha0_dDelta_dTau2(self) -> float: ...

    def d3alpha0_dDelta2_dTau(self) -> float: ...

    def d3alpha0_dDelta3(self) -> float: ...

    def alphar(self) -> float: ...

    def dalphar_dDelta(self) -> float: ...

    def dalphar_dTau(self) -> float: ...

    def d2alphar_dDelta2(self) -> float: ...

    def d2alphar_dDelta_dTau(self) -> float: ...

    def d2alphar_dTau2(self) -> float: ...

    def d3alphar_dDelta3(self) -> float: ...

    def d3alphar_dDelta2_dTau(self) -> float: ...

    def d3alphar_dDelta_dTau2(self) -> float: ...

    def d3alphar_dTau3(self) -> float: ...

    def d4alphar_dDelta4(self) -> float: ...

    def d4alphar_dDelta3_dTau(self) -> float: ...

    def d4alphar_dDelta2_dTau2(self) -> float: ...

    def d4alphar_dDelta_dTau3(self) -> float: ...

    def d4alphar_dTau4(self) -> float: ...

def get_config_as_json_string() -> str: ...

def set_config_as_json_string(arg: str, /) -> None: ...

@overload
def config_key_description(arg: configuration_keys, /) -> str: ...

@overload
def config_key_description(arg: str, /) -> str: ...

def set_config_string(arg0: configuration_keys, arg1: str, /) -> None: ...

def set_config_double(arg0: configuration_keys, arg1: float, /) -> None: ...

def set_departure_functions(arg: str, /) -> None: ...

def set_config_bool(arg0: configuration_keys, arg1: bool, /) -> None: ...

def get_config_string(arg: configuration_keys, /) -> str: ...

def get_config_double(arg: configuration_keys, /) -> float: ...

def get_config_bool(arg: configuration_keys, /) -> bool: ...

def get_config_int(arg: configuration_keys, /) -> int: ...

def set_config_int(arg0: configuration_keys, arg1: int, /) -> None: ...

def get_parameter_information(arg0: int, arg1: str, /) -> str: ...

def get_parameter_index(arg: str, /) -> parameters: ...

def get_phase_index(arg: str, /) -> phases: ...

def is_trivial_parameter(arg: int, /) -> bool: ...

def generate_update_pair(arg0: parameters, arg1: float, arg2: parameters, arg3: float, /) -> tuple: ...

def Props1SI(arg0: str, arg1: str, /) -> float: ...

@overload
def PropsSI(Output: str, FluidName: str, /) -> float: ...
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str, /) -> float: ...  # type: ignore[overload-overlap]
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Num, Name2: str, Prop2: _Num, FluidName: str, /) -> NDArray[numpy.float64]: ...
@overload
def PropsSI(Output: list[str] | tuple[str, ...] | NDArray[numpy.str_], Name1: str, Prop1: _Num, Name2: str, Prop2: _Num, FluidName: str, /) -> NDArray[numpy.float64]: ...

def PhaseSI(arg0: str, arg1: float, arg2: str, arg3: float, arg4: str, /) -> str: ...

def PropsSImulti(arg0: Sequence[str], arg1: str, arg2: Sequence[float], arg3: str, arg4: Sequence[float], arg5: str, arg6: Sequence[str], arg7: Sequence[float], /) -> list[list[float]]: ...

def get_global_param_string(arg: str, /) -> str: ...

def get_debug_level() -> int: ...

def set_debug_level(arg: int, /) -> None: ...

def get_fluid_param_string(arg0: str, arg1: str, /) -> str: ...

def extract_backend(arg: str, /) -> tuple: ...

def extract_fractions(arg: str, /) -> tuple: ...

def set_reference_stateS(arg0: str, arg1: str, /) -> None: ...

def set_reference_stateD(arg0: str, arg1: float, arg2: float, arg3: float, arg4: float, /) -> None: ...

def saturation_ancillary(arg0: str, arg1: str, arg2: int, arg3: str, arg4: float, /) -> float: ...

def add_fluids_as_JSON(arg0: str, arg1: str, /) -> bool: ...

@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Scalar, Name2: str, Value2: _Scalar, Name3: str, Value3: _Scalar, /) -> float: ...  # type: ignore[overload-overlap]
@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Num, Name2: str, Value2: _Num, Name3: str, Value3: _Num, /) -> NDArray[numpy.float64]: ...

def HAProps_Aux(arg0: str, arg1: float, arg2: float, arg3: float, /) -> tuple: ...

def cair_sat(arg: float, /) -> float: ...

def get_mixture_binary_pair_data(arg0: str, arg1: str, arg2: str, /) -> str: ...

def set_mixture_binary_pair_data(arg0: str, arg1: str, arg2: str, arg3: float, /) -> None: ...

def apply_simple_mixing_rule(arg0: str, arg1: str, arg2: str, /) -> None: ...

def set_interaction_parameters(arg: str, /) -> None: ...

def set_predefined_mixtures(arg: str, /) -> None: ...

def get_mixture_binary_pair_pcsaft(arg0: str, arg1: str, arg2: str, /) -> str: ...

def set_mixture_binary_pair_pcsaft(arg0: str, arg1: str, arg2: str, arg3: float, /) -> None: ...

def FluidsList() -> list[str]: ...

def get_aliases(arg: str, /) -> list[str]: ...

def get_REFPROPname(arg: str, /) -> str: ...

def get_BibTeXKey(arg0: str, arg1: str, /) -> str: ...

def get_errstr() -> str: ...

@overload
def set_reference_state(FluidName: str, reference_state: str, /) -> None: ...
@overload
def set_reference_state(FluidName: str, T0: float, rhomolar: float, hmolar0: float, smolar0: float, /) -> None: ...

class MonotonicExpansionMatch:
    @property
    def idx(self) -> int: ...

    @property
    def ymin(self) -> float: ...

    @property
    def ymax(self) -> float: ...

    @property
    def xmin(self) -> float: ...

    @property
    def xmax(self) -> float: ...

class IntervalMatch:
    @property
    def expansioninfo(self) -> list[MonotonicExpansionMatch]: ...

    @property
    def xmin(self) -> float: ...

    @property
    def xmax(self) -> float: ...

    @property
    def ymin(self) -> float: ...

    @property
    def ymax(self) -> float: ...

class ChebyshevExpansion:
    def __init__(self, xmin: float, xmax: float, coef: Sequence[float]) -> None: ...

    def xmin(self) -> float: ...

    def xmax(self) -> float: ...

    def coeff(self) -> list[float]: ...

    def eval_many(self, x: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], y: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]) -> None: ...

    def solve_for_x(self, y: float, a: float, b: float, bits: int, max_iter: int, boundstytol: float) -> float: ...

    def solve_for_x_many(self, y: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], a: float, b: float, bits: int, max_iter: int, boundstytol: float, x: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], counts: Annotated[NDArray[numpy.uint64], dict(shape=(None,), order='C')]) -> None: ...

class ChebyshevApproximation1D:
    def __init__(self, expansions: Sequence[ChebyshevExpansion]) -> None: ...

    def xmin(self) -> float: ...

    def xmax(self) -> float: ...

    def is_monotonic(self) -> bool: ...

    def eval_many(self, x: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], y: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]) -> None: ...

    def get_x_for_y(self, y: float, bits: int, max_iter: int, boundstytol: float) -> list[tuple[float, int]]: ...

    def count_x_for_y_many(self, y: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], bits: int, max_iter: int, boundstytol: float, counts: Annotated[NDArray[numpy.uint64], dict(shape=(None,), order='C')]) -> None: ...

    def monotonic_intervals(self) -> list[IntervalMatch]: ...

class SuperAncillary:
    def __init__(self, json_as_string: str) -> None: ...

    def eval_sat(self, T: float, prop: str, Q: int) -> float: ...

    def eval_sat_many(self, T: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], prop: str, Q: int, y: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]) -> None: ...

_capi: object = ...
