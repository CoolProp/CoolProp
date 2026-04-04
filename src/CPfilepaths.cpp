#include "PlatformDetermination.h"
#include "Exceptions.h"
#include "CPfilepaths.h"
#include "CPstrings.h"
#include "CoolPropTools.h"

#include <fstream>
#include <algorithm>
#include <cerrno>
#include <filesystem>

// This will kill the horrible min and max macros
#ifndef NOMINMAX
#    define NOMINMAX
#endif

namespace fs = std::filesystem;

// Platform-specific includes needed by get_home_dir
#if defined(__ISWINDOWS__)
#    define UNICODE
#    define _UNICODE
#    include "Windows.h"
#    include <windows.h>
#else
#    include <sys/types.h>
#    include <sys/stat.h>
#    include <unistd.h>
#    if !defined(__powerpc__)
#        include <pwd.h>
#    endif
#endif

// ---- CalculateDirSize -------------------------------------------------------

unsigned long long CalculateDirSize(const std::string& path) {
    unsigned long long size = 0;
    std::error_code ec;
    for (const auto& entry : fs::recursive_directory_iterator(path, fs::directory_options::skip_permission_denied, ec)) {
        if (!ec && entry.is_regular_file(ec) && !ec) {
            size += entry.file_size(ec);
            ec.clear();
        }
    }
    return size;
}

// Legacy implementations kept for reference (powerpc / pre-NDK-r23 Android):
//
// #if defined(__ISWINDOWS__)
// /// From http://stackoverflow.com/a/17827724/1360263
// static bool IsBrowsePath(const std::wstring& path) {
//     return (path == L"." || path == L"..");
// }
// unsigned long long CalculateDirSize(const std::wstring& path, std::vector<std::wstring>* errVect) {
//     unsigned long long size = 0;
//     WIN32_FIND_DATAW data;
//     HANDLE sh = NULL;
//     sh = FindFirstFileW((path + L"\\*").c_str(), &data);
//     if (sh == INVALID_HANDLE_VALUE) {
//         if (errVect != NULL) errVect->push_back(path);
//         return size;
//     }
//     do {
//         if (!IsBrowsePath(data.cFileName)) {
//             if (data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
//                 size += CalculateDirSize(path + L"\\" + data.cFileName, NULL);
//             else
//                 size += data.nFileSizeHigh * (unsigned long long)(MAXDWORD) + data.nFileSizeLow;
//         }
//     } while (FindNextFileW(sh, &data));
//     FindClose(sh);
//     return size;
// }
// #else
// // Android / powerpc fallback — stub returning 0
// unsigned long long CalculateDirSize(const std::string& path) { return 0; }
// #endif

// ---- get_binary_file_contents -----------------------------------------------

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

// ---- make_dirs --------------------------------------------------------------

void make_dirs(std::string file_path) {
    std::error_code ec;
    fs::create_directories(fs::path(file_path), ec);
    // Ignore errors (same behaviour as legacy code which silently skipped failures)
}

// Legacy implementation kept for reference (powerpc / pre-NDK-r23 Android):
//
// void make_dirs(std::string file_path) {
//     std::replace(file_path.begin(), file_path.end(), '\\', '/');
// #    if defined(__ISWINDOWS__)
//     const char sep = '\\';
// #    else
//     const char sep = '/';
// #    endif
//     std::vector<std::string> pathsplit = strsplit(file_path, '/');
//     std::string path = pathsplit[0];
//     for (std::size_t i = 0, sz = pathsplit.size(); i < sz; i++) {
//         if (!path_exists(path)) {
// #        if defined(__ISWINDOWS__)
//             int errcode = CreateDirectoryA((LPCSTR)path.c_str(), NULL);
//             if (errcode == 0) {
//                 switch (GetLastError()) {
//                     case ERROR_ALREADY_EXISTS: break;
//                     case ERROR_PATH_NOT_FOUND:
//                         throw CoolProp::ValueError(format("Unable to make the directory %s", path.c_str()));
//                     default: break;
//                 }
//             }
// #        elif !defined(__powerpc__)
//             mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
// #        endif
//         }
//         if (i < (sz - 1)) path += format("%c%s", sep, pathsplit[i + 1].c_str());
//     }
// }

// ---- get_separator ----------------------------------------------------------

std::string get_separator(void) {
    return std::string(1, fs::path::preferred_separator);
}

// Legacy implementation kept for reference (powerpc / pre-NDK-r23 Android):
//
// std::string get_separator(void) {
// #if defined(__ISLINUX__)
//     return std::string("/");
// #elif defined(__ISAPPLE__)
//     return std::string("/");
// #elif defined(__ISWINDOWS__)
//     return std::string("\\");
// #else
//     throw CoolProp::NotImplementedError("This function is not defined for your platform.");
// #endif
// }

// ---- get_home_dir -----------------------------------------------------------

std::string get_home_dir(void) {
// See http://stackoverflow.com/questions/2552416/how-can-i-find-the-users-home-dir-in-a-cross-platform-manner-using-c
#if defined(__ISLINUX__) || defined(__ISAPPLE__)
    char* home = NULL;
    home = getenv("HOME");
#    if defined(__ISAPPLE__)
    if (home == NULL) {
        struct passwd* pwd = getpwuid(getuid());
        if (pwd) {
            home = pwd->pw_dir;
        }
    }
#    endif
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

// ---- path_exists ------------------------------------------------------------

bool path_exists(const std::string& path) {
    std::error_code ec;
    return fs::exists(fs::path(path), ec);
}

// Legacy implementation kept for reference (powerpc / pre-NDK-r23 Android):
//
// bool path_exists(const std::string& path) {
//     std::string path_cpy;
//     if (endswith(path, get_separator())) {
//         path_cpy = path.substr(0, path.size() - 1);
//     } else {
//         path_cpy = path;
//     }
// #    if defined(__ISWINDOWS__)
//     struct _stat buf;
//     return (_stat(path_cpy.c_str(), &buf) == 0);
// #    elif defined(__ISLINUX__) || defined(__ISAPPLE__)
//     struct stat st;
//     if (lstat(path_cpy.c_str(), &st) == 0)
//         return S_ISDIR(st.st_mode) || S_ISREG(st.st_mode);
//     return false;
// #    else
//     throw CoolProp::NotImplementedError("This function is not defined for your platform.");
// #    endif
// }

// ---- join_path --------------------------------------------------------------

std::string join_path(const std::string& one, const std::string& two) {
    return (fs::path(one) / two).string();
}

// Legacy implementation kept for reference (powerpc / pre-NDK-r23 Android):
//
// std::string join_path(const std::string& one, const std::string& two) {
//     std::string result;
//     std::string separator = get_separator();
//     if (!endswith(one, separator) && !one.empty()) {
//         result = one + separator;
//     } else {
//         result = one;
//     }
//     result.append(two);
//     return result;
// }

// ---- get_file_contents ------------------------------------------------------

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
