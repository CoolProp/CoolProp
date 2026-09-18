# The CoolProp wrapper for EES

Developer reference for `wrappers/EES`. EES comes in a 32-bit and a 64-bit
flavour and the wrapper is built for both. Everything below applies to both
unless a paragraph says otherwise.

## 1. What EES expects from an external function

Sources (F-Chart EES help):

- External functions: <https://fchartsoftware.com/ees/eeshelp/external_functions.htm>
- File types and filename extensions: <https://fchartsoftware.com/ees/eeshelp/ees_file_types_and_filename_extensions.htm>
- EES 64-bit Professional License: <https://fchartsoftware.com/ees/eeshelp/hs713.htm>
  and <https://fchartsoftware.com/ees/64-bit.php>
- DLL file skeleton in Visual C++: <https://fchartsoftware.com/ees/eeshelp/dll_file_skeleton_in_visual_c__.htm>
- Example .DLF external function in C++: <https://fchartsoftware.com/ees/eeshelp/example_.dlf_external_function_in_c__.htm>
- Library files: <https://fchart.com/ees/eeshelp/library_files.htm>
- EES_REFPROP installation, which shows the folder layout for both bitnesses:
  <https://fchartsoftware.com/ees/ees_refprop/installation.htm>

**One library per bitness.** The 32-bit program loads `.dlf` and `.LIB` files
from `USERLIB`, the 64-bit program (`EES64.exe`) loads `.dlf64` and `.LIB64`
files from `USERLIB64`. Neither can read the library of the other, so the
wrapper is compiled twice. Per add-on the files live in a subfolder of the user
library folder, which is `C:\EES32\Userlib\COOLPROP_EES` and
`C:\EES64\Userlib64\COOLPROP_EES` in a default installation. The `EES32` folder
name is historical and is used by the 32-bit program whatever the operating
system.

**The function name is the file name.** EES derives the name of the external
function from the file name, so the libraries have to be called
`COOLPROP_EES.dlf` and `COOLPROP_EES.dlf64`, and the exported symbol has to be
an undecorated `COOLPROP_EES`.

**The call.** Both bitnesses use the same signature:

```cpp
__declspec(dllexport) double COOLPROP_EES(char fluid[256], int& mode, struct EesParamRec* input_rec)
```

with `struct EesParamRec { double value; struct EesParamRec* next; }`. `fluid`
is a 255-character C string plus its terminator, used in both directions. The
linked list carries the inputs and ends on a null `next`. On x64 the record is
16 bytes, which is what both EES (Delphi) and this C++ produce.

**The mode says what EES wants and what it got.** EES passes it by reference:

| Mode in | EES is asking for | The wrapper answers |
|---|---|---|
| `-1` | an example of the call format | `T = PropsSI('T','P',101325,'Q',0,'Water')` |
| `-2` | the units of the inputs | empty string |
| `-3` | the units of the output | empty string |
| other | the calculation | see below |

On the way back from a calculation the mode carries the status: 0 with an empty
string for a normal result, a positive value with a message for an error, which
makes EES stop the calculation and show the message, and a negative value with a
message for a warning, which does not stop it. The three requests answer in the
string and leave the mode as EES set it, the way the F-Chart example does.

Taking `mode` by value reads the low bits of the pointer instead, which serves
none of the requests and lets no error message reach the user.

**Calling convention.** 32-bit Windows has several, and EES uses the C one, so
the 32-bit build defines `CONVENTION=__cdecl`. x64 has a single convention and
no name decoration, so the 64-bit build defines nothing.

**Precision.** The 64-bit program uses 64-bit doubles where the 32-bit program
uses the 80-bit extended type. The wrapper exchanges plain doubles, so this
changes round-off inside EES only, not inside CoolProp.

## 2. What is shipped

Four files per bitness, installed into the folder named above:

| 32-bit | 64-bit | What it is |
|---|---|---|
| `COOLPROP_EES.dlf` | `COOLPROP_EES.dlf64` | the external function, CoolProp linked in statically |
| `CoolProp.LIB` | `CoolProp.LIB64` | EES source: the `PropsSI` and `PropsSIZ` wrappers, the unit-system check, and the string protocol that packs the fluid and the property keys into the one string argument |
| `CoolProp.htm` | `CoolProp.htm` | the help file EES shows for the library |
| `CoolProp_EES_Sample.EES` | `CoolProp_EES_Sample.EES` | the example file |

Only one library name is shipped per folder. Both names in one folder would make
EES define every function twice. The `.LIB64` file is a copy of the `.LIB` file,
the content does not depend on the bitness.

**The `.LIB` file is not plain text.** EES wraps the source in a 35-byte header
and a 14-byte trailer. Four of those header bytes, at offset 31, are the length
of the text body as a little-endian unsigned integer, so a hand edit that changes
the length has to write that field back. `dev/ci/check_ees_artifacts.py` checks
it. The remaining header and trailer bytes are not decoded and are left alone.
If EES ever refuses the file, open it in EES and save it again, which rewrites
the whole wrapper.

## 3. Building it

The wrapper is Windows-only; a configure on any other system is refused. The
bitness comes from the generator, and `BITNESS` then drives the artefact suffix,
the EES library name, the calling convention and the install folder:

```
git clone https://github.com/CoolProp/CoolProp
mkdir CoolProp/build && cd CoolProp/build
cmake .. -G "Visual Studio 17 2022" -A Win32 -DCOOLPROP_EES_MODULE=ON
cmake --build . --target COOLPROP_EES --config Release
```

```
cd .. && mkdir build64 && cd build64
cmake .. -G "Visual Studio 17 2022" -A x64 -DCOOLPROP_EES_MODULE=ON
cmake --build . --target COOLPROP_EES --config Release
```

Each build leaves the four files of its bitness in the build directory.
`cmake --install` puts them under `EES/<system>/32bit` or `EES/<system>/64bit`,
relative to the install prefix, so the two do not overwrite each other.

To debug, point the output at the user library folder and start EES from the
debugger: set the output file to
`C:\EES32\Userlib\COOLPROP_EES\COOLPROP_EES.dlf` (or
`C:\EES64\Userlib64\COOLPROP_EES\COOLPROP_EES.dlf64`), set the debug command to
`c:\EES32\ees` (or `c:\EES64\ees64`), set a breakpoint in `COOLPROP_EES` and run
the project.

## 4. Packaging and release

`COOLPROP_WINDOWS_PACKAGE_EES` and `COOLPROP_WINDOWS_PACKAGE_EES64` run the two
sub-builds with `-A Win32` and `-A x64` and stage their output in
`InnoScript/source/EES` and `InnoScript/source/EES64`.
`COOLPROP_WINDOWS_PACKAGE_INSTALLER` depends on both, so the Windows installer
job builds both libraries. `windows_installer.yml` runs it on a push to the
mainline branches, on `v*` tags and on a pull request against those branches,
and the file is in the builder list of `release_all_files.yml`, which is what
the nightly and the tagged file drops collect. The libraries are also published
on their own, as the `EES` artifact with one subfolder per bitness.

The installer offers one task per flavour, so a machine with only one EES
installed gets only the files it can use. The shipped installer is built from
`CoolProp/ExcelAddinInstaller`, which `cmake/dependencies.cmake` pins by commit;
`wrappers/EES/BuildInnoInstaller.iss.in` is a stand-alone script for a local
package of the two libraries and is not wired into the build.

`dev/ci/check_ees_artifacts.py` runs after the build and fails the job unless
both libraries exist, their PE headers report i386 and amd64 respectively, each
one exports an undecorated `COOLPROP_EES`, and each `CoolProp.LIB` agrees with
its own header and still defines `propssi`, `propssiz` and
`coolprop_assert_si_units` and none of the removed functions. EES ignores a
library of the wrong bitness or with a decorated export without saying anything,
so this is checked before release rather than by the user. It runs as part of
the Windows installer job, so it covers the triggers listed above and not every
branch. What it inspects is the `InnoScript/source` tree that the packaging
targets fill, before it copies the files on to the artifact folder, and never
the compiled installer; a missing or misnamed file in the installer script fails
the ISCC step of the build instead.

## 5. Behaviour worth knowing

- **`PropsSI` and `PropsSIZ` are the only functions.** The deprecated
  `coolprop()` and `coolpropsi()` were removed, see below.
- **A failed call stops the EES calculation**, with a positive mode and the
  message in the string, the way F-Chart documents it. The alternative is
  returning 0 and letting the solve continue with that number.
- **EES does not unit-check a `COOLPROP_EES` call.** The units of the two
  arguments depend on the property keys encoded in the fluid string, which EES
  does not pass when it asks for units, so the wrapper answers modes -2 and -3
  with an empty string. F-Chart: "EES will provide unit checking if the external
  procedure provides the strings of inputs and outputs requested with Modes=-2
  and -3. Otherwise it will skip unit checking for this equation." The unit
  system itself is still checked, by `CoolProp.LIB` (see
  `coolprop_assert_si_units`).
- **One unit system, checked before the call.** The string protocol carries no
  units, so `CoolProp.LIB` asserts the EES unit system before it calls in, and
  the last field of the string names the system the numbers are in:

  | Function | Tag in the string | C++ path | Unit system asserted |
  |---|---|---|---|
  | `PropsSI`, `PropsSIZ` | `SI` | `PropsSI` / `PropsSImulti` | K, Pa, J, mass |

  The wrapper refuses any other tag with a message about mismatched releases.
  Nothing in the shipped library file sends one, so a tag other than `SI` means
  an old `CoolProp.LIB` is sitting next to a new `COOLPROP_EES`.

  **Why `coolprop()` and `coolpropsi()` are gone.** `coolprop()` tagged its
  string `kSI` and reached the deprecated v4 `Props` API, which reads kPa and kJ,
  while the library file asserted the *SI* system. A model that did what the
  error message asked handed pressures in Pa to a function reading kPa and got a
  state a factor of 1000 away without any complaint. `coolpropsi()` could not run
  at all: it called `COOLPROP_EES(f6$, ...)` and `f6$` was never assigned, and
  behind that sat a key-translation block it never used and a set of 1000x
  conversions that contradicted its own `SI` tag. The v4 `Props` call went with
  them; `PropsSI` covers what they did.
- **`$DEBUG` writes to the working directory of the EES process.** Append
  `'$DEBUG'` to the fluid name and the wrapper writes `log.txt` and
  `log_stdout.txt` by relative path, so they land where EES runs, not in the
  user library folder.
- **The sample file stays a `.EES` file.** The 64-bit program reads it; only
  large files are slow to convert.

## 6. Checking a build

CI covers the bitness and the exported name; everything below needs a licence of
EES. Run it for a bitness after a change to `main.cpp` or `CoolProp.LIB`:

1. Build the target as in section 3 and confirm that `dumpbin /exports` lists an
   undecorated `COOLPROP_EES`. This is the manual equivalent of the CI gate.
2. Copy the four files into the matching user library folder and start EES. The
   function appears in the Function Information dialog, with the example call
   from mode -1.
3. Run `CoolProp_EES_Sample.EES` and compare the numbers against the other
   bitness.
4. Give a fluid name that CoolProp rejects: the calculation stops and shows the
   CoolProp message instead of returning 0.
5. Press F8. `Check Units` reports the `COOLPROP_EES` equation as skipped rather
   than flagging it.
6. Exercise the unit assertion, which is the one part of `CoolProp.LIB` with no
   automated cover. On K, Pa, J and mass a `PropsSI` call goes through; set the
   unit system to kPa or kJ and it stops with the message from
   `coolprop_assert_si_units`.
7. Confirm that EES accepted the library file: the Function Information dialog
   lists `PropsSI` and `PropsSIZ` and nothing else. This is the check that
   matters after any edit of the wrapped `.LIB` format.
8. Append `'$DEBUG'` to a fluid name and look for `log.txt` and `log_stdout.txt`
   in the working directory of the EES process.
