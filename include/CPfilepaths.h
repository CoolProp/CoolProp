#ifndef COOLPROP_FILE_PATH_H
#define COOLPROP_FILE_PATH_H

#include <string>
#include <vector>

/// Get directory separator
std::string get_separator(void);

/// Get the user's home directory;  It is believed that is is always a place that files can be written
std::string get_home_dir(void);

/// Return true if path exists
bool path_exists(const std::string& path);

/// Return merged path, append separator if string two is empty
std::string join_path(const std::string& one, const std::string& two);

/// Make directory and all required intermediate directories
void make_dirs(std::string file_path);

/// Get the size of a directory in bytes
#if defined(__ISWINDOWS__)
unsigned long long CalculateDirSize(const std::wstring& path, std::vector<std::wstring>* errVect = NULL);
#else
unsigned long long CalculateDirSize(const std::string& path);
#endif

// Get all the contents of a file and dump into a STL string
// Thanks to http://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
std::string get_file_contents(const char* filename);

/// Get all the contents of a binary file
std::vector<char> get_binary_file_contents(const char* filename);

#endif