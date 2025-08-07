import CoolProp
from CoolProp.Plots import PropertyPlot
import CoolProp.CoolProp as CP
import math
import sys

def format_number(value):
    """Format a number for C++ output, handling NaN and scientific notation"""
    if math.isnan(value):
        return 'std::nan("")'
    elif abs(value) >= 1e7 or (abs(value) <= 1e-4 and value != 0):
        # Use scientific notation with high precision for very large/small numbers
        return f'{value:.15e}'
    else:
        # Use regular format with high precision 
        return f'{value:.15g}'

def generate_value_at_test():
    print('TEST_CASE("Check value_at for p-h plots", "[Plot]") {')
    print('    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);')
    print()
    print('    CHECK_THAT(plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));')
    print('    CHECK_THAT(plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));')
    print(f'    CHECK_THAT(plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs({CP.PropsSI("T", "P",200000, "H",300000, "R134a")}, 1e-10));')
    print(f'    CHECK_THAT(plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs({CP.PropsSI("Q", "P",200000, "H",300000, "R134a")}, 1e-10));')
    print('}')

def generate_ph_test():
    print('TEST_CASE("Check that the isolines are the same as from Python", "[Plot]") {')
    print('    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);')
    print('    const int isoline_count = 5;')
    print('    const int points_per_isoline = 5;')
    print()
    print('    // CHECK(plot.xkey_ == CoolProp::iHmass);')
    print('    // CHECK(plot.ykey_ == CoolProp::iP);')
    print()
    
    # Create PH plot
    plot = PropertyPlot('HEOS::R134a', 'ph', unit_system='SI', tp_limits='ACHP')
    plot.calc_isolines(num=5, points=5)
    
    # Axis properties
    axis_values = plot.axis.axis()
    print('    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);')
    print('    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Log);')
    print(f'    CHECK_THAT(plot.xaxis.min, WithinAbs({axis_values[0]}, 1));')
    print(f'    CHECK_THAT(plot.xaxis.max, WithinAbs({axis_values[1]}, 1));')
    print(f'    CHECK_THAT(plot.yaxis.min, WithinAbs({axis_values[2]}, 1));')
    print(f'    CHECK_THAT(plot.yaxis.max, WithinAbs({axis_values[3]}, 1));')
    print()
    
    # Supported isoline types
    print('    std::vector<CoolProp::parameters> iso_types = plot.supported_isoline_keys();')
    print('    REQUIRE(iso_types.size() == 4);')
    print('    CHECK(iso_types[0] == CoolProp::iT);')
    print('    CHECK(iso_types[1] == CoolProp::iQ);')
    print('    CHECK(iso_types[2] == CoolProp::iDmass);')
    print('    CHECK(iso_types[3] == CoolProp::iSmass);')
    print()
    
    # Generate isolines for each type
    isolines_map = [
        ('Q', 'iQ', CoolProp.iQ),
        ('T', 'iT', CoolProp.iT), 
        ('S', 'iSmass', CoolProp.iSmass),
        ('D', 'iDmass', CoolProp.iDmass)
    ]
    
    for cpp_name, cpp_param, py_param in isolines_map:
        print('    {')
        print(f'        // {cpp_name} isolines')
        print(f'        CoolProp::Plot::Range {cpp_name.lower()}_range = plot.isoline_range(CoolProp::{cpp_param});')
        print(f'        std::vector<double> {cpp_name.lower()}_values = CoolProp::Plot::generate_values_in_range(CoolProp::{cpp_param}, {cpp_name.lower()}_range, isoline_count);')
        print(f'        CoolProp::Plot::Isolines {cpp_name.lower()}_isolines = plot.calc_isolines(CoolProp::{cpp_param}, {cpp_name.lower()}_values, points_per_isoline);')
        print(f'        REQUIRE({cpp_name.lower()}_isolines.size() == isoline_count);')
        
        # Generate value checks
        for i in range(5):
            value = plot._isolines[py_param][i].value
            print(f'        CHECK_THAT({cpp_name.lower()}_isolines[{i}].value, WithinAbs({value}, 1e-10));')
        
        # Generate expected x values
        print('        const double expected_x[isoline_count][points_per_isoline] = {')
        for i in range(5):
            x_values = [format_number(x) for x in plot._isolines[py_param][i].x]
            print('            {' + ', '.join(x_values) + '},')
        print('        };')
        
        # Generate expected y values
        print('        const double expected_y[isoline_count][points_per_isoline] = {')
        for i in range(5):
            y_values = [format_number(y) for y in plot._isolines[py_param][i].y]
            print('            {' + ', '.join(y_values) + '},')
        print('        };')
        
        # Generate validation loop
        print(f'        for (int i = 0; i < {cpp_name.lower()}_isolines.size(); ++i) {{')
        print(f'            REQUIRE({cpp_name.lower()}_isolines[i].size() == points_per_isoline);')
        print(f'            for (int j = 0; j < {cpp_name.lower()}_isolines[i].size(); ++j) {{')
        
        # Check for NaN values in this isoline type
        has_nan_x = any(math.isnan(x) for isoline in plot._isolines[py_param] for x in isoline.x)
        
        if has_nan_x:
            print(f'                if (std::isnan({cpp_name.lower()}_isolines[i].x[j]))')
            print('                    CHECK(std::isnan(expected_x[i][j]));')
            print('                else')
            print(f'                    CHECK_THAT({cpp_name.lower()}_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));')
            print(f'                CHECK_THAT({cpp_name.lower()}_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));')
        else:
            print(f'                CHECK_THAT({cpp_name.lower()}_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));')
            print(f'                CHECK_THAT({cpp_name.lower()}_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));')
            
        print('            }')
        print('        }')
        print('    }')
    
    print('}')
    print()

def generate_ts_test():
    print('TEST_CASE("Basic TS Plot has same output as Python", "[Plot]") {')
    print('    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);')
    print('    const int isoline_count = 5;')
    print('    const int points_per_isoline = 5;')
    print()
    print('    // CHECK(plot.xkey_ == CoolProp::iSmass);')
    print('    // CHECK(plot.ykey_ == CoolProp::iT);')
    print()
    
    # Create TS plot
    plot = PropertyPlot('HEOS::R134a', 'Ts', unit_system='SI', tp_limits='ACHP')
    plot.calc_isolines(num=5, points=5)
    
    # Axis properties
    axis_values = plot.axis.axis()
    print('    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);')
    print('    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Lin);')
    print(f'    CHECK_THAT(plot.xaxis.min, WithinAbs({axis_values[0]}, 1));')
    print(f'    CHECK_THAT(plot.xaxis.max, WithinAbs({axis_values[1]}, 1));')
    print(f'    CHECK_THAT(plot.yaxis.min, WithinAbs({axis_values[2]}, 1));')
    print(f'    CHECK_THAT(plot.yaxis.max, WithinAbs({axis_values[3]}, 1));')
    print()
    
    # Supported isoline types
    print('    std::vector<CoolProp::parameters> iso_types = plot.supported_isoline_keys();')
    print('    REQUIRE(iso_types.size() == 4);')
    print('    CHECK(iso_types[0] == CoolProp::iP);')
    print('    CHECK(iso_types[1] == CoolProp::iQ);')
    print('    CHECK(iso_types[2] == CoolProp::iDmass);')
    print('    CHECK(iso_types[3] == CoolProp::iHmass);')
    print()
    
    # Generate isolines for each type
    isolines_map = [
        ('Q', 'iQ', CoolProp.iQ),
        ('P', 'iP', CoolProp.iP),
        ('H', 'iHmass', CoolProp.iHmass),
        ('D', 'iDmass', CoolProp.iDmass)
    ]
    
    for cpp_name, cpp_param, py_param in isolines_map:
        print('    {')
        print(f'        // {cpp_name} isolines')
        print(f'        CoolProp::Plot::Range {cpp_name.lower()}_range = plot.isoline_range(CoolProp::{cpp_param});')
        print(f'        std::vector<double> {cpp_name.lower()}_values = CoolProp::Plot::generate_values_in_range(CoolProp::{cpp_param}, {cpp_name.lower()}_range, isoline_count);')
        print(f'        CoolProp::Plot::Isolines {cpp_name.lower()}_isolines = plot.calc_isolines(CoolProp::{cpp_param}, {cpp_name.lower()}_values, points_per_isoline);')
        print(f'        REQUIRE({cpp_name.lower()}_isolines.size() == isoline_count);')
        
        # Generate value checks
        for i in range(5):
            value = plot._isolines[py_param][i].value
            if cpp_name == 'P':
                print(f'        CHECK_THAT({cpp_name.lower()}_isolines[{i}].value, WithinAbs({value}, 1e-7));')
            else:
                print(f'        CHECK_THAT({cpp_name.lower()}_isolines[{i}].value, WithinAbs({value}, 1e-10));')
        
        # Add blank line only for the first Q isolines section
        if cpp_name == 'Q':
            print()
        
        # Generate expected x values
        print('        const double expected_x[isoline_count][points_per_isoline] = {')
        for i in range(5):
            x_values = [format_number(x) for x in plot._isolines[py_param][i].x]
            print('            {' + ', '.join(x_values) + '},')
        print('        };')
        
        # Generate expected y values
        print('        const double expected_y[isoline_count][points_per_isoline] = {')
        for i in range(5):
            y_values = [format_number(y) for y in plot._isolines[py_param][i].y]
            print('            {' + ', '.join(y_values) + '},')
        print('        };')
        
        # Add blank line only for the first Q isolines section
        if cpp_name == 'Q':
            print()
        
        # Generate validation loop
        print(f'        for (int i = 0; i < {cpp_name.lower()}_isolines.size(); ++i) {{')
        print(f'            REQUIRE({cpp_name.lower()}_isolines[i].size() == points_per_isoline);')
        print(f'            for (int j = 0; j < {cpp_name.lower()}_isolines[i].size(); ++j) {{')
        
        # Check for NaN values in this isoline type
        has_nan_y = any(math.isnan(y) for isoline in plot._isolines[py_param] for y in isoline.y)
        
        if has_nan_y:
            print(f'                if (std::isnan(expected_y[i][j])) {{')
            print(f'                    CHECK(std::isnan({cpp_name.lower()}_isolines[i].y[j]));')
            print('                } else {')
            print(f'                    CHECK_THAT({cpp_name.lower()}_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));')
            print(f'                    CHECK_THAT({cpp_name.lower()}_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));')
            print('                }')
        else:
            print(f'                CHECK_THAT({cpp_name.lower()}_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));')
            print(f'                CHECK_THAT({cpp_name.lower()}_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));')
            
        print('            }')
        print('        }')
        print('    }')
    
    print('}')

if __name__ == '__main__':

    from contextlib import redirect_stdout
    import io
    from pathlib import Path
    import re

    f = io.StringIO()
    with redirect_stdout(f):
        generate_value_at_test()
        generate_ph_test()
        generate_ts_test()
    captured_output = f.getvalue()
    if len(captured_output) == 0:
        raise ValueError('no captured output in stdout')
    file_path = (Path(__file__).parent.parent.parent / 'src' / "CoolPropPlot.cpp")
    with file_path.open('r') as fp:
        contents = fp.read()
        pattern = re.compile(r"\/\/ <autogenerated>(.*)\/\/ <\/autogenerated>", flags=re.MULTILINE | re.DOTALL)
        # for match in re.findall(pattern, contents):
        #     print(match)
        contents_new = re.sub(pattern, f'// <autogenerated>\n{captured_output}\n// </autogenerated>', contents)
        if contents_new == contents:
            raise ValueError("No changes were made")
    with file_path.open('w') as fp:
        fp.write(contents_new)