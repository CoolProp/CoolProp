#include "PlatformDetermination.h"
#include "Exceptions.h"
#include "CPfilepaths.h"
#include "CPstrings.h"
#include "CoolPropTools.h"

#include <fstream>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <cerrno>

// This will kill the horrible min and max macros
#ifndef NOMINMAX
#    define NOMINMAX
#endif

#if defined(__ISWINDOWS__)
#    define UNICODE
#    define _UNICODE
#    include "Windows.h"
#    include <windows.h>  // for the CreateDirectory function
#else
#    include <unistd.h>
#    if !defined(__powerpc__)
#        include <pwd.h>
#    endif
#endif

#if defined(__ISWINDOWS__)
/// From http://stackoverflow.com/a/17827724/1360263
bool IsBrowsePath(const std::wstring& path) {
    return (path == L"." || path == L"..");
}
unsigned long long CalculateDirSize(const std::wstring& path, std::vector<std::wstring>* errVect) {
    unsigned long long size = 0;
    WIN32_FIND_DATAW data;
    HANDLE sh = NULL;
    sh = FindFirstFileW((path + L"\\*").c_str(), &data);

    if (sh == INVALID_HANDLE_VALUE) {
        //if we want, store all happened error
        if (errVect != NULL) errVect->push_back(path);
        return size;
    }

    do {
        // skip current and parent
        if (!IsBrowsePath(data.cFileName)) {
            // if found object is ...
            if (data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
                // directory, then search it recursievly
                size += CalculateDirSize(path + L"\\" + data.cFileName, NULL);
            else
                // otherwise get object size and add it to directory size
                size += data.nFileSizeHigh * (unsigned long long)(MAXDWORD) + data.nFileSizeLow;
        }

    } while (FindNextFileW(sh, &data));  // do

    FindClose(sh);

    return size;
}
#elif defined(__ANDROID__) || defined(__powerpc__)
// Android doesn't have ftw.h, also doesn't accept not having this file
unsigned long long CalculateDirSize(const std::string& path) {
    return 0;
}
#else
#    include <ftw.h>
#    include <stdint.h>
#    include <iostream>

static thread_local unsigned long long ftw_summer;  // An evil global variable for the ftw function
int ftw_function(const char* fpath, const struct stat* sb, int tflag, struct FTW* ftwbuf) {
    ftw_summer += sb->st_size;
    return 0; /* To tell nftw() to continue */
}
unsigned long long CalculateDirSize(const std::string& path) {
    ftw_summer = 0;
    int flags = 0 | FTW_DEPTH | FTW_PHYS;
    nftw(path.c_str(), ftw_function, 20, flags);
    double temp = ftw_summer;
    ftw_summer = 0;
    return temp;
}
#endif

std::vector<char> get_binary_file_contents(const char* filename) {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in) {
        std::vector<char> contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned int)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return (contents);
    }
    throw(errno);
}

void make_dirs(std::string file_path) {
    std::replace(file_path.begin(), file_path.end(), '\\', '/');  // replace all '\' with '/'

#if defined(__ISWINDOWS__)
    const char sep = '\\';  // well, Windows (and DOS) allows forward slash "/", too :)
#else
    const char sep = '/';
#endif

    std::vector<std::string> pathsplit = strsplit(file_path, '/');
    std::string path = pathsplit[0];  // will throw if pathsplit.size() == 0
    for (std::size_t i = 0, sz = pathsplit.size(); i < sz; i++) {
        if (!path_exists(path)) {
#if defined(__ISWINDOWS__)  // Defined for 32-bit and 64-bit windows
            int errcode = CreateDirectoryA((LPCSTR)path.c_str(), NULL);
            if (errcode == 0) {
                switch (GetLastError()) {
                    case ERROR_ALREADY_EXISTS:
                        break;
                    case ERROR_PATH_NOT_FOUND:
                        throw CoolProp::ValueError(format("Unable to make the directory %s", path.c_str()));
                    default:
                        break;
                }
            }
#else
#    if defined(__powerpc__)
#    else
            mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#    endif
#endif
        }
        if (i < (sz - 1)) path += format("%c%s", sep, pathsplit[i + 1].c_str());
    }
};

std::string get_separator(void) {
#if defined(__ISLINUX__)
    return std::string("/");
#elif defined(__ISAPPLE__)
    return std::string("/");
#elif defined(__ISWINDOWS__)
    return std::string("\\");
#else
    throw CoolProp::NotImplementedError("This function is not defined for your platform.");
#endif
}

std::string get_home_dir(void) {
// See http://stackoverflow.com/questions/2552416/how-can-i-find-the-users-home-dir-in-a-cross-platform-manner-using-c
#if defined(__ISLINUX__)
    char* home = NULL;
    home = getenv("HOME");
    return std::string(home);
#elif defined(__ISAPPLE__)
    char* home = NULL;
    home = getenv("HOME");
    if (home == NULL) {
        struct passwd* pwd = getpwuid(getuid());
        if (pwd) {
            home = pwd->pw_dir;
        }
    }
    if (home == NULL) {
        throw CoolProp::NotImplementedError("Could not detect home directory.");
    }
    return std::string(home);
#elif defined(__ISWINDOWS__)
#    if defined(_MSC_VER)
#        pragma warning(push)
#        pragma warning(disable : 4996)
#    endif
    char* pUSERPROFILE = getenv("USERPROFILE");
    if (pUSERPROFILE != NULL) {
        return std::string(pUSERPROFILE);
    } else {
        char* pHOMEDRIVE = getenv("HOMEDRIVE");
        char* pHOMEPATH = getenv("HOMEPATH");
        if (pHOMEDRIVE != NULL && pHOMEPATH != NULL) {
            return std::string(pHOMEDRIVE) + std::string(pHOMEPATH);
        } else {
            return std::string("");
        }
    }
#    if defined(_MSC_VER)
#        pragma warning(pop)
#    endif
#else
    throw CoolProp::NotImplementedError("This function is not defined for your platform.");
#endif
};

bool path_exists(const std::string& path) {
    std::string path_cpy;
    if (endswith(path, get_separator())) {
        path_cpy = path.substr(0, path.size() - 1);
    } else {
        path_cpy = path;
    }

#if defined(__ISWINDOWS__)  // Defined for 32-bit and 64-bit windows
    struct _stat buf;
    // Get data associated with path using the windows libraries,
    // and if you can (result == 0), the path exists
    if (_stat(path_cpy.c_str(), &buf) == 0) {
        return true;
    } else {
        return false;
    }
#elif defined(__ISLINUX__) || defined(__ISAPPLE__)
    struct stat st;
    if (lstat(path_cpy.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) return true;
        if (S_ISREG(st.st_mode)) return true;
        return false;
    } else {
        return false;
    }
#else
    throw CoolProp::NotImplementedError("This function is not defined for your platform.");
#endif
};

std::string join_path(const std::string& one, const std::string& two) {
    std::string result;
    std::string separator = get_separator();
    if (!endswith(one, separator) && !one.empty()) {
        result = one + separator;
    } else {
        result = one;
    }
    result.append(two);
    return result;
}

std::string get_file_contents(const char* filename) {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned int)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return (contents);
    }
    throw(errno);
}
