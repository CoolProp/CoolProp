#include "CPstrings.h"

std::string strjoin(const std::vector<std::string> &strings, const std::string &delim)
{
    // Empty input vector
    if (strings.empty()){return "";}

    std::string output = strings[0];
    for (unsigned int i = 1; i < strings.size(); i++)
    {
        output += format("%s%s",delim.c_str(),strings[i].c_str());
    }
    return output;
}

std::vector<std::string> strsplit(const std::string &s, char del)
{
    std::vector<std::string> v;
    std::string::const_iterator i1 = s.begin(), i2;
    while (true){
        i2 = std::find(i1, s.end(), del);
        v.push_back(std::string(i1, i2));
        if (i2 == s.end())
            break;
        i1 = i2+1;
    }
    return v;
}