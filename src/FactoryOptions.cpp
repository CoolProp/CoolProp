#include "CoolProp/FactoryOptions.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>

#include "CoolProp/detail/filepaths.h"
#include "CoolProp/Exceptions.h"

namespace CoolProp {

namespace {

bool is_all_whitespace(const std::string& s) {
    return std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c) != 0; });
}

}  // namespace

FactoryOptionsParseResult parse_factory_options(const std::string& factory_string) {
    FactoryOptionsParseResult result;

    const std::size_t pos = factory_string.find_first_of('?');
    if (pos == std::string::npos) {
        result.clean_string = factory_string;
        return result;
    }

    result.clean_string = factory_string.substr(0, pos);
    std::string tail = factory_string.substr(pos + 1);

    if (tail.empty() || is_all_whitespace(tail)) {
        return result;
    }

    if (tail.front() == '@') {
        const std::string path = tail.substr(1);
        if (path.empty()) {
            throw ValueError("factory-string options: '@' indirection requires a non-empty path");
        }
        // get_file_contents() reads the file verbatim and throws on
        // read failure.  The path may contain '?' / '&' / spaces /
        // etc.  — no further parsing of the path string.
        result.options_json = get_file_contents(path.c_str());
        return result;
    }

    result.options_json = std::move(tail);
    return result;
}

}  // namespace CoolProp
