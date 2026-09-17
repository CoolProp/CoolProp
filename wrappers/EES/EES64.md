# Adding a 64-bit library for the EES wrapper

Status: implemented in the build, packaging and release workflow, September
2026. The parts that need a Windows machine with both EES licences are listed
in section 5 and are still open.

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

## 2. State before this change

| Item | Location | 32-bit assumption |
|---|---|---|
| Wrapper source | `wrappers/EES/main.cpp` | none, except the `mode` signature (see below) |
| EES library file | `wrappers/EES/CoolProp.LIB` | name only |
| Help file | `wrappers/EES/CoolProp.htm` | none |
| Sample | `wrappers/EES/CoolProp_EES_Sample.EES` | none |
| Build target | `CMakeLists.txt`, EES module block | hard `FATAL_ERROR` for 64-bit, `.dlf` suffix, `-m32` |
| Windows package | `CMakeLists.txt`, `COOLPROP_WINDOWS_PACKAGE_EES` | sub-build forced to `-AWin32` |
| Stand-alone installer | `wrappers/EES/BuildInnoInstaller.iss.in` | `c:\ees32\Userlib\COOLPROP_EES`, `COOLPROP_EES.dlf` |
| Combined installer | `CoolProp/ExcelAddinInstaller`, `addin-installer.iss` and `cmake-templates/config.iss` | single `EESINSDIR` pointing at `C:\EES32\Userlib\COOLPROP_EES` |

Two defects found while reading the build:

- The EES module aborted the configure step with "You cannot build the EES
  wrapper as a 64-bit library." That was the only hard blocker, the rest was
  packaging.
- For non-MSVC compilers the module **assigned** `COMPILE_FLAGS` as `-m32`,
  which dropped the `-DCOOLPROP_LIB -DCONVENTION=__cdecl` set a few lines
  earlier instead of appending to it. Both are fixed.

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

The mode is also the channel for the result status. The F-Chart help says it
plainly:

> Under normal operation, S is returned from the function as the null string and
> the function should set Mode to 0. [On error] S should be set to an
> appropriate error message and Mode should be set to a positive integer. EES
> will then terminate calculations and display this error message. [For a
> non-fatal warning] Mode should be set to a negative integer and the warning
> should appear in string S.

Because the wrapper never wrote to `mode`, none of its error messages could
reach the user either. It now declares `int& mode` and follows the contract:

- `-1` returns the example call, `-2` and `-3` return an empty string because
  the units of the two inputs depend on the property keys encoded in the fluid
  string,
- a normal call returns the null string with `mode` set to 0,
- every error path sets `mode` to 1 and leaves its message in the string,
- a CoolProp warning sets `mode` to -1 and leaves the warning in the string.

The error case is a real behaviour change: a failing call used to return 0 and
let the EES solve continue with that number, it now stops the calculation and
shows the message. That is what the documented contract asks for, and a silent
zero is the worse of the two.

**This still needs one manual run in EES before it ships, in 32-bit as much as
in 64-bit.** If the F-Chart documentation were wrong and EES really passed the
mode by value, the dereference would fault on every call. The 32-bit library is
therefore not bit-identical to the one before this change. The change is kept in
its own commit so it can be reverted on its own.

## 4. What was changed

### 4.1 `CMakeLists.txt`, EES module block

1. The `FATAL_ERROR` is gone, the artefact name follows `BITNESS`: `.dlf` for
   32-bit, `.dlf64` for 64-bit.
2. `-m32` is only applied for a 32-bit non-MSVC build, and appended instead of
   replacing `COMPILE_FLAGS`.
3. `-DCONVENTION=__cdecl` stays for 32-bit. The 64-bit build defines nothing and
   lets `CoolPropLib.h` pick the empty default, the same way the shared library
   does.
4. `CoolProp.LIB` is copied to the build directory as `CoolProp.LIB64` for a
   64-bit build (same content, different name). Only one of the two names is
   shipped per folder, EES would otherwise load the same functions twice.
5. The install destination gained a bitness folder,
   `${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}/<32|64>bit`, so the two
   builds no longer overwrite each other.
6. A non-Windows configure is refused up front. The old 64-bit `FATAL_ERROR`
   happened to block that too, now it is checked on purpose.

For MSVC the 32-bit compile flags come out exactly as before. For a 32-bit
MinGW build they do not: the old code dropped `-DCOOLPROP_LIB` (see the second
defect above), so that library exported `COOLPROP_EES` alone, and it now also
exports the CoolProp C API, the same way the MSVC build always did.

### 4.2 Windows package and release

`COOLPROP_WINDOWS_PACKAGE_EES64` runs the same sub-build with `-A x64` into
`EES64/` and copies the result to `InnoScript/source/EES64`.
`COOLPROP_WINDOWS_PACKAGE_INSTALLER` depends on both EES targets, so the
Windows installer job builds both bitnesses. That job runs on every push to
`master` and on `v*` tags, and `windows_installer.yml` is already in the builder
list of `release_all_files.yml`, which is what the nightly and the tagged file
drops collect.

`dev/ci/check_ees_artifacts.py` verifies after the build that both libraries
exist, that their PE headers really report i386 and amd64, and that each one
exports an undecorated `COOLPROP_EES`, which is the name EES derives from the
file name. It then stages the eight files for the new `EES` artifact, and both
uploads use `if-no-files-found: error`. A library of the wrong bitness or with a
decorated export is ignored by EES without a message, so this gate fails the job
rather than shipping a broken wrapper.

What the gate does not cover: it inspects the staged files, not the compiled
installer, and it does not notice a stale `InnoScript/source` tree from an
earlier local build. In CI the tree is always fresh.

### 4.3 Installers

- `wrappers/EES/BuildInnoInstaller.iss.in` installs both sets, one task per
  bitness, into `c:\EES32\Userlib\COOLPROP_EES` and
  `c:\EES64\Userlib64\COOLPROP_EES`. Note that this script is not wired into the
  build, the shipped installer comes from the repository below.
- `CoolProp/ExcelAddinInstaller` (separate repository): `config.iss` gained
  `EESINSDIR64`, `addin-installer.iss` the four `Source:` lines and the
  `EesUserLib64` task, `messages.iss` the task descriptions in all three
  languages. That change lives on the branch `chp/ees-64bit-userlib64`, commit
  `0e3974c`, and `cmake/dependencies.cmake` pins exactly that commit. Re-pin to
  the merge commit once the branch lands on master there. Pinning a commit that
  is not on the default branch works because the package is fetched without
  `GIT_SHALLOW`, so the full history is cloned before the checkout. **Do not
  delete that branch before the re-pin has landed**: a squash merge followed by
  the usual branch deletion makes the commit unreachable, and every Windows
  package build then fails at the configure step, tagged releases included.

### 4.4 Documentation

`Web/coolprop/wrappers/EES/index.rst` and `wrappers/EES/README.rst` describe the
two flavours, the `Userlib64` paths and the 64-bit build command.

## 5. Verification

Nothing here can be verified on Linux, the wrapper is Windows-only. The
checklist for a Windows machine with both EES licences:

1. `cmake -G "Visual Studio 17 2022" -A x64 .. -DCOOLPROP_EES_MODULE=ON` and
   build the `COOLPROP_EES` target, confirm `COOLPROP_EES.dlf64` is produced and
   that `dumpbin /exports` shows an undecorated `COOLPROP_EES`. The CI gate
   checks both, this is the manual equivalent.
2. Copy `COOLPROP_EES.dlf64`, `CoolProp.LIB64` and `CoolProp.htm` into
   `C:\EES64\Userlib64\COOLPROP_EES`, start `EES64.exe`, check that the function
   shows up in the Function Information dialog.
3. Run `CoolProp_EES_Sample.EES` in EES64 and compare the numbers against the
   32-bit run.
4. Run the same sample with the 32-bit build. The mode fix changes that library
   too, so it needs the same pass, not just a rebuild.
5. Exercise the mode contract in both flavours: the Function Information dialog
   must show the example call (mode -1), a bad fluid name must stop the
   calculation with the CoolProp message rather than return 0 (positive mode),
   and the units shown for the arguments must be acceptable, since the empty
   answer to mode -2 and -3 is our reading of the documentation, not something
   F-Chart spells out. Watch in particular whether any of the three requests
   comes back as a warning: they answer in the string and leave the mode at the
   value EES passed in, following the F-Chart example, and a negative mode is
   also what a warning looks like on a normal call.
6. Check the `$DEBUG` path, it writes `log.txt` and `log_stdout.txt` into the
   working directory of the EES process.

## 6. Open questions

- Does EES64 accept a `.LIB` file in `USERLIB64`, or is `.LIB64` mandatory? The
  help says the `.LIB64` set is auto-loaded, and that the 64-bit program can read
  32-bit files. We ship `.LIB64` only, because shipping both names in the same
  folder would define every function twice.
- Should the sample be re-saved as `CoolProp_EES_Sample.EES64`? Not required, the
  64-bit program reads `.EES` files, but the conversion is reported to be slow
  for large files.
- The `mode` fix of section 3 changes 32-bit behaviour as well and is the one
  part of this change that cannot be checked without a running EES.
