import timeit
import CoolProp.CoolProp as CP


def time_check(N, h, p, TTSE=False, mode='TTSE'):

    if TTSE:
        if mode == 'TTSE':
            setup = "import CoolProp; import CoolProp.CoolProp as CP; CP.enable_TTSE_LUT('Water'); CP.set_TTSE_mode('Water','TTSE'); CP.Props('T','H',500,'P',10000,'Water'); IWater = CP.get_Fluid_index('Water'); from CoolProp.param_constants import iT,iH,iP,iD"
        elif mode == 'BICUBIC':
            setup = "import CoolProp; import CoolProp.CoolProp as CP; CP.enable_TTSE_LUT('Water'); CP.set_TTSE_mode('Water','BICUBIC'); CP.Props('T','H',500,'P',10000,'Water'); IWater = CP.get_Fluid_index('Water'); from CoolProp.param_constants import  iT,iH,iP,iD"
        else:
            raise ValueError()
    else:
        setup = "import CoolProp.CoolProp as CP; IWater = CP.get_Fluid_index('Water'); CP.disable_TTSE_LUT('Water'); from CoolProp.param_constants import  iT,iH,iP,iD"

    time = timeit.Timer("CP.IProps(iD,iH," + str(h) + ",iP," + str(p) + ",IWater)", setup).timeit(N) / N * 1e6
    value = CP.Props('D', 'H', h, 'P', p, 'Water')

    return time, value


values = dict(subcooled=(500, 10000), twophase=(2000, 10000), superheated=(3000, 10000), supercritical=(2000, 30000))

N = 10000
for k in ['subcooled', 'twophase', 'superheated', 'supercritical']:
    h, p = values[k]

    time_EOS, value_EOS = time_check(N, h, p, TTSE=False)
    time_TTSE, value_TTSE = time_check(N, h, p, TTSE=True)
    time_BICUBIC, value_BICUBIC = time_check(N, h, p, TTSE=True, mode='BICUBIC')

    print("%s %s %s %s %s %s %s" % (k, h, p, (value_TTSE / value_EOS - 1.0) * 100, (value_BICUBIC / value_EOS - 1.0) * 100, time_EOS / time_TTSE, time_EOS / time_BICUBIC))
