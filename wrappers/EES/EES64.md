# Adding a 64-bit library for the EES wrapper

Status: investigation / implementation plan, September 2026.

This note collects what F-Chart documents about 64-bit EES, what the current
CoolProp EES wrapper does, and what has to change to ship a 64-bit external
function next to the existing 32-bit one.

## 1. What F-Chart requires for 64-bit EES

Sources (F-Chart EES help):

- EES 64-bit Professional License: <https://fchartsoftware.com/ees/eeshelp/hs713.htm>
  and <https://fchartsoftware.com/ees/64-bit.php>
- File types and filename extensions: <https://fchartsoftware.com/ees/eeshelp/ees_file_types_and_filename_extensions.htm>
- External functions: <https://fchartsoftware.com/ees/eeshelp/external_functions.htm>
- DLL file skeleton in Visual C++: <https://fchartsoftware.com/ees/eeshelp/dll_file_skeleton_in_visual_c__.htm>
- Example .DLF external function in C++: <https://fchartsoftware.com/ees/eeshelp/example_.dlf_external_function_in_c__.htm>
- Library files: <https://fchart.com/ees/eeshelp/library_files.htm>
- EES_REFPROP installation (shows the folder layout for both bitnesses):
  <https://fchartsoftware.com/ees/ees_refprop/installation.htm>

The relevant facts:

1. **Separate binaries, separate extensions.** The 64-bit license runs
   `EES64.exe` and refuses 32-bit external libraries. External routines have to
   be recompiled with a 64-bit compiler and saved with `.DLF64` (function),
   `.DLP64` (procedure), `.FDL64` (FORTRAN) or `.DLL64` (several routines in one
   DLL). Our wrapper is a single external function, so the target is
   `COOLPROP_EES.DLF64`.
2. **Separate library folder.** Files placed in `USERLIB64` are loaded
   automatically at startup by the 64-bit program, `USERLIB` is the 32-bit
   equivalent. The EES_REFPROP installation page spells the layout out as
   `EES32\USERLIB\EES_REFPROP` and `EES64\USERLIB64\EES_REFPROP`, that is, one
   subfolder per add-on inside the user library folder. The default install
   folder of the 32-bit program is `C:\EES32` for historical reasons.
3. **EES library files are versioned too.** The auto-loaded set in `USERLIB64` is
   documented as `.LIB64`, `.FDL64`, `.DLF64`, `.DLP64` and `.DLL64`, so our
   `CoolProp.LIB` (plain EES source, it only wraps the external call and does the
   unit checks) should also be shipped as `CoolProp.LIB64`. The 64-bit program
   can read 32-bit files, so the content does not have to change, only the name.
4. **The function name is the file name.** EES takes the name of a `.DLF`/`.DLF64`
   from the file name, so the file has to stay `COOLPROP_EES.DLF64` for the
   existing `CoolProp.LIB` calls to keep resolving.
5. **Calling convention.** x64 Windows has a single calling convention, so
   `__cdecl` and `__stdcall` collapse into one and no name decoration happens.
   Nothing special is needed, but the `-DCONVENTION=__cdecl` that the EES target
   passes today should become empty for a 64-bit build, in line with the
   `CONVENTION` logic in `CMakeLists.txt:701-719`.
6. **No ABI change in the call itself.** The documented signature stays
   `double FUNC(char s[256], int& mode, struct EesParamRec* input_rec)` with
   `struct EesParamRec { double value; struct EesParamRec* next; }`. On x64 the
   record is 16 bytes (8 byte double, 8 byte pointer), which is what both EES
   (Delphi) and our C++ produce, so the linked list walk in `main.cpp` is fine.
7. **Precision.** The 64-bit program uses 64-bit doubles instead of the 80-bit
   extended type of the 32-bit program. We exchange plain doubles, so this only
   means slightly different round-off in EES itself, not in CoolProp.

## 2. Current state in this repository

| Item | Location | 32-bit assumption |
|---|---|---|
| Wrapper source | `wrappers/EES/main.cpp` | none, except the `mode` signature (see below) |
| EES library file | `wrappers/EES/CoolProp.LIB` | name only |
| Help file | `wrappers/EES/CoolProp.htm` | none |
| Sample | `wrappers/EES/CoolProp_EES_Sample.EES` | none |
| Build target | `CMakeLists.txt:1149-1212` | hard `FATAL_ERROR` for 64-bit, `.dlf` suffix, `-m32` |
| Windows package | `CMakeLists.txt:1386-1404` | sub-build forced to `-AWin32` |
| Stand-alone installer | `wrappers/EES/BuildInnoInstaller.iss.in:21,34` | `c:\ees32\Userlib\COOLPROP_EES`, `COOLPROP_EES.dlf` |
| Combined installer | `CoolProp/ExcelAddinInstaller`, `addin-installer.iss:91-94` and `cmake-templates/config.iss:19` | single `EESINSDIR` pointing at `C:\EES32\Userlib\COOLPROP_EES` |

Notable details found while reading the build:

- `CMakeLists.txt:1150-1153` aborts the configure step with
  "You cannot build the EES wrapper as a 64-bit library." That is the only hard
  blocker, the rest is packaging.
- `CMakeLists.txt:1169-1172` sets `COMPILE_FLAGS` to `-m32` for non-MSVC
  compilers, which **overwrites** the `-DCOOLPROP_LIB -DCONVENTION=__cdecl` set a
  few lines earlier instead of appending to it. Worth fixing while touching the
  block.

## 3. ABI finding: the `mode` argument

F-Chart documents the second argument as passed **by reference**, both in the
Delphi skeleton (`var Mode: integer`) and in the C++ skeleton and example
(`int& mode`). Our wrapper declares it by value:

```cpp
// wrappers/EES/main.cpp:81
__declspec(dllexport) double COOLPROP_EES(char fluid[256], int mode, struct EesParamRec* input_rec)
```

Consequences:

- The stack/register slot is the same size on both bitnesses, so nothing is
  corrupted, but `mode` holds the low bits of a pointer instead of the mode
  value. The `mode == -1` branch (return an example call string for the Function
  Information dialog) can therefore never be taken, and `mode == -2` / `-3` (unit
  strings for inputs and outputs) are not handled at all.
- This is not new, the same signature is in the museum version of the wrapper,
  so the 32-bit library has always behaved this way.

Recommendation: fix it to `int& mode` as part of, or just before, the 64-bit
work, but **verify it on a machine with EES installed first**. If the
documentation were wrong and EES really passed the mode by value, the
dereference would fault, so this needs one manual run in EES 32-bit and EES
64-bit before it ships. If `mode` is honoured, the `-2` and `-3` cases must also
be answered (or at least answered with an empty string and `mode` set to 0),
otherwise EES will take whatever we leave in the buffer as a unit string.

This is independent of the 64-bit port and can be split into its own change.

## 4. Proposed changes

### 4.1 `CMakeLists.txt`, EES module block

1. Drop the `FATAL_ERROR` and derive the artefact name from `BITNESS`:
   `.dlf` for 32-bit, `.dlf64` for 64-bit.
2. Only apply `-m32` for a 32-bit non-MSVC build, and append to the existing
   `COMPILE_FLAGS` instead of replacing them.
3. Use `-DCONVENTION=` (empty) for the 64-bit build, `-DCONVENTION=__cdecl`
   stays for 32-bit.
4. Copy `CoolProp.LIB` to the build directory as `CoolProp.LIB64` for a 64-bit
   build (same content, different name), and install into
   `${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}/64bit` so the two bitnesses
   do not overwrite each other.

### 4.2 Windows package

Add a `COOLPROP_WINDOWS_PACKAGE_EES64` target next to
`COOLPROP_WINDOWS_PACKAGE_EES` (`CMakeLists.txt:1386-1404`) that runs the same
sub-build with `-A x64` into a separate binary directory and copies the result to
`InnoScript/source/EES64`. `COOLPROP_WINDOWS_PACKAGE_INSTALLER` then depends on
both.

### 4.3 Installers

- `wrappers/EES/BuildInnoInstaller.iss.in`: parameterise the default directory
  and the file list, or add a second script for the 64-bit case. The 64-bit
  target directory is `C:\EES64\Userlib64\COOLPROP_EES`.
- `CoolProp/ExcelAddinInstaller` (separate repository, read-only from here):
  `cmake-templates/config.iss` needs a second define, for example
  `#define EESINSDIR64 "C:\EES64\Userlib64\COOLPROP_EES"`, and
  `addin-installer.iss` needs four more `Source:` lines under a second task
  (`EesUserLib64`) for `CoolProp.htm`, `CoolProp.LIB64`, `COOLPROP_EES.dlf64` and
  the sample. This has to go in as a separate pull request there, and the
  `GIT_TAG` in `cmake/dependencies.cmake:120` has to be bumped afterwards.
- Both tasks should ideally only be offered when the matching EES folder exists.
  Inno Setup can check that with `DirExists()` in a `Check:` parameter, which
  avoids installing a 64-bit library for a user who only has the 32-bit program.

### 4.4 Documentation

`Web/coolprop/wrappers/EES/index.rst` needs the 64-bit build command, the
`USERLIB64` paths and the debugging instructions for `EES64.exe`.
`wrappers/EES/README.rst` mentions `c:\EES32\Userlib` only.

## 5. Verification

Nothing here can be verified on Linux, the wrapper is Windows-only. The
checklist for a Windows machine with both EES licences:

1. `cmake -G "Visual Studio 17 2022" -A x64 .. -DCOOLPROP_EES_MODULE=ON` and
   build the `COOLPROP_EES` target, confirm `COOLPROP_EES.dlf64` is produced and
   that `dumpbin /exports` shows an undecorated `COOLPROP_EES`.
2. Copy `COOLPROP_EES.dlf64`, `CoolProp.LIB64` and `CoolProp.htm` into
   `C:\EES64\Userlib64\COOLPROP_EES`, start `EES64.exe`, check that the function
   shows up in the Function Information dialog.
3. Run `CoolProp_EES_Sample.EES` in EES64 and compare the numbers against the
   32-bit run.
4. Repeat the 32-bit build to make sure the existing artefact is unchanged.
5. Check the `$DEBUG` path, it writes `log.txt` and `log_stdout.txt` into the
   working directory of the EES process.

## 6. Open questions

- Does EES64 accept a `.LIB` file in `USERLIB64`, or is `.LIB64` mandatory? The
  help says the `.LIB64` set is auto-loaded, and that the 64-bit program can read
  32-bit files, but this should be confirmed by experiment. Shipping the copy
  under both names costs 7 kB and removes the doubt.
- Should the sample be re-saved as `CoolProp_EES_Sample.EES64`? Not required, the
  64-bit program reads `.EES` files, but the conversion is reported to be slow
  for large files.
- Is the `mode` fix (section 3) wanted in the same change, or separately? It
  changes 32-bit behaviour as well.
