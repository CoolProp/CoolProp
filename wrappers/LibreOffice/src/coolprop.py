import uno
import unohelper

from org.coolprop.wrappers.libreoffice.CalcAddIn import XCoolProp

try:
    from CoolProp import CoolProp
except:
    pass


class CoolPropCalcAddin(unohelper.Base, XCoolProp):
    """CoolProp AddIn for LibreOffice Calc.

    This class wraps the interface declared in XCoolProp.idl to the appropriate
    python function calls.
    """

    def __init__(self, ctx):
        self.ctx = ctx

    def PropsSI(self, output, name1, prop1, name2, prop2, fluid_name):
        """Calculate fluid properties for a given state from SI inputs."""
        try:
            return CoolProp.PropsSI(output, name1, prop1, name2, prop2, fluid_name)
        except Exception as e:
            return str(e)

    def Props1SI(self, fluid_name, output):
        """Get trivial fluid properties (e.g. critical temperature)."""
        try:
            return CoolProp.PropsSI(fluid_name, output)
        except Exception as e:
            return str(e)

    def PhaseSI(self, name1, prop1, name2, prop2, fluid_name):
        """Return the phase of a given state from SI inputs"""
        try:
            return CoolProp.PhaseSI(name1, prop1, name2, prop2, fluid_name)
        except Exception as e:
            return str(e)

    def HAPropsSI(self, output, name1, prop1, name2, prop2, name3, prop3):
        """Calculate properties of humid air from SI inputs."""
        try:
            return CoolProp.HAPropsSI(output, name1, prop1, name2, prop2, name3, prop3)
        except Exception as e:
            return str(e)

    def get_fluid_param_string(self, fluid_name, param_name):
        """Get fluid parameter string from CoolProp."""
        try:
            return CoolProp.get_fluid_param_string(fluid_name, param_name)
        except Exception as e:
            return str(e)

    def get_global_param_string(self, param_name, split=False):
        """Get global parameter string from CoolProp.

        The split option was added for convenience and can be used to output
        comma separated values as a column vector in Calc. Thus, the function
        can be directly used as input for dropdown fields.
        """
        try:
            param_str = CoolProp.get_global_param_string(param_name)
            if split:
                return zip(param_str.split(','))
            else:
                return [[param_str]]
        except Exception as e:
            return [[str(e)]]


g_ImplementationHelper = unohelper.ImplementationHelper()
g_ImplementationHelper.addImplementation(CoolPropCalcAddin,
                                         "org.coolprop.wrappers.libreoffice.CalcAddIn",
                                         ("com.sun.star.sheet.AddIn",),
                                        )
