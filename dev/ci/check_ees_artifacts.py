#!/usr/bin/env python
"""Verify the EES wrapper artifacts of the Windows package.

EES loads external functions by bitness: the 32-bit program reads
COOLPROP_EES.dlf from its Userlib folder, the 64-bit program (EES64.exe)
reads COOLPROP_EES.dlf64 from Userlib64.  EES ignores a library of the wrong
bitness without a message, so a build that produced two files of the same
machine type would ship a broken 64-bit wrapper unnoticed.

For each bitness this checks that the four files exist, that the PE header
reports the expected machine type, that the library exports an undecorated
COOLPROP_EES, and that CoolProp.LIB agrees with its own header and defines
what it should.  With --stage the checked files are also copied into a clean
folder, which is what the release workflow uploads.

Usage:
    python dev/ci/check_ees_artifacts.py <source-dir> [--stage <dir>]

where <source-dir> is the "InnoScript/source" folder of the Windows package
build, the one holding the EES and EES64 subfolders.
"""

import argparse
import os
import shutil
import struct
import sys

# Machine types as defined by the PE format, see the Microsoft PE documentation
IMAGE_FILE_MACHINE_I386 = 0x014C
IMAGE_FILE_MACHINE_AMD64 = 0x8664

MACHINE_NAMES = {
    IMAGE_FILE_MACHINE_I386: "i386 (32-bit)",
    IMAGE_FILE_MACHINE_AMD64: "amd64 (64-bit)",
}

# EES derives the name of the external function from the file name
EXPORTED_FUNCTION = "COOLPROP_EES"

# The EES library file is not plain text: the source sits between a header and
# a trailer written by EES, and four header bytes hold the length of that text.
# An edit that changes the length without writing the field back leaves EES
# reading the wrong number of bytes.  These offsets hold for the file as EES
# writes it; a different header fails the check, which is the safe direction.
EES_LIB_HEADER_LEN = 35
EES_LIB_TRAILER_LEN = 14
EES_LIB_LENGTH_OFFSET = 31

# A correct length only proves the header agrees with the body, not that the
# body is the library, so the body is checked for what it must and must not
# define.  The version string and the undecoded header fields are deliberately
# not pinned, because a newer EES may write a different header.
EES_LIB_REQUIRED = ["Function propssi(", "Function propssiz(", "Function coolprop_assert_si_units("]
EES_LIB_REMOVED = ["Function coolprop(", "Function coolpropsi(", "coolprop_assert_ksi_units", "coolprop_assert_cpsi_units"]

# What a malformed image can raise while it is inspected.  struct.error is not
# a ValueError, so it has to be listed separately.
PE_ERRORS = (ValueError, OSError, struct.error)

# One entry per bitness: folder, the library file and the expected machine
LAYOUT = [
    {
        "folder": "EES",
        "library": "COOLPROP_EES.dlf",
        "ees_library": "CoolProp.LIB",
        "machine": IMAGE_FILE_MACHINE_I386,
        "files": ["COOLPROP_EES.dlf", "CoolProp.LIB", "CoolProp.htm", "CoolProp_EES_Sample.EES"],
    },
    {
        "folder": "EES64",
        "library": "COOLPROP_EES.dlf64",
        "ees_library": "CoolProp.LIB64",
        "machine": IMAGE_FILE_MACHINE_AMD64,
        "files": ["COOLPROP_EES.dlf64", "CoolProp.LIB64", "CoolProp.htm", "CoolProp_EES_Sample.EES"],
    },
]


# check_entry only reports a missing file for the names in "files", so both
# libraries have to stay in that list; otherwise renaming one would turn its
# check into a silent pass.
for _entry in LAYOUT:
    assert _entry["ees_library"] in _entry["files"], "{0} is not in the file list of {1}".format(_entry["ees_library"], _entry["folder"])
    assert _entry["library"] in _entry["files"], "{0} is not in the file list of {1}".format(_entry["library"], _entry["folder"])


def _pe_offset(image):
    """Return the offset of the PE header, raising ValueError if there is none."""
    if len(image) < 0x40 or image[0:2] != b"MZ":
        raise ValueError("no MZ signature, this is not a Windows library")
    # The offset of the PE header is stored at 0x3C in the DOS header
    (offset,) = struct.unpack_from("<I", image, 0x3C)
    if offset + 24 > len(image) or image[offset : offset + 4] != b"PE\0\0":
        raise ValueError("no PE signature at offset {0}".format(offset))
    return offset


def read_pe_machine(path):
    """Return the machine type stored in the PE header of path.

    Raises one of PE_ERRORS if the file is not a PE image.
    """
    with open(path, "rb") as handle:
        image = handle.read()
    (machine,) = struct.unpack_from("<H", image, _pe_offset(image) + 4)
    return machine


def read_pe_exports(path):
    """Return the names exported by the PE image at path.

    Walks the export directory by hand so that the check needs no third party
    module on the build agent.  Raises one of PE_ERRORS when the image cannot
    be read, which the caller reports as a failure: a library whose exports
    cannot be inspected is never silently accepted.
    """
    with open(path, "rb") as handle:
        image = handle.read()
    pe = _pe_offset(image)
    section_count, = struct.unpack_from("<H", image, pe + 6)
    optional_size, = struct.unpack_from("<H", image, pe + 20)
    optional = pe + 24
    magic, = struct.unpack_from("<H", image, optional)
    if magic == 0x10B:  # PE32
        directories = optional + 96
    elif magic == 0x20B:  # PE32+, the 64-bit variant
        directories = optional + 112
    else:
        raise ValueError("unknown optional header magic 0x{0:04x}".format(magic))
    export_rva, export_size = struct.unpack_from("<II", image, directories)
    if export_rva == 0 or export_size == 0:
        return []

    # The section table follows the optional header and maps the addresses
    # used inside the image (RVAs) to offsets in the file.
    sections = []
    for index in range(section_count):
        entry = optional + optional_size + index * 40
        virtual_size, virtual_address, raw_size, raw_offset = struct.unpack_from("<IIII", image, entry + 8)
        sections.append((virtual_address, max(virtual_size, raw_size), raw_offset))

    def to_offset(rva):
        for virtual_address, size, raw_offset in sections:
            if virtual_address <= rva < virtual_address + size:
                return rva - virtual_address + raw_offset
        raise ValueError("address 0x{0:08x} is outside every section".format(rva))

    export_dir = to_offset(export_rva)
    name_count, = struct.unpack_from("<I", image, export_dir + 24)
    names_rva, = struct.unpack_from("<I", image, export_dir + 32)
    if name_count == 0:
        return []
    names_table = to_offset(names_rva)

    names = []
    for index in range(name_count):
        (name_rva,) = struct.unpack_from("<I", image, names_table + index * 4)
        start = to_offset(name_rva)
        end = image.find(b"\0", start)
        if end < 0:
            raise ValueError("unterminated export name at offset {0}".format(start))
        names.append(image[start:end].decode("ascii", "replace"))
    return names


def check_ees_library(path, errors):
    """Check the header of an EES library file against its own contents."""
    with open(path, "rb") as handle:
        blob = handle.read()
    wrapper = EES_LIB_HEADER_LEN + EES_LIB_TRAILER_LEN
    if len(blob) <= wrapper:
        errors.append("{0} is {1} bytes, too short to hold an EES header and a body".format(path, len(blob)))
        return
    (declared,) = struct.unpack_from("<I", blob, EES_LIB_LENGTH_OFFSET)
    actual = len(blob) - wrapper
    if declared != actual:
        errors.append("{0} declares a body of {1} bytes in its header but holds {2}".format(path, declared, actual))
    else:
        print("ok: {0} declares its body length of {1} bytes correctly".format(path, actual))

    body = blob[EES_LIB_HEADER_LEN:-EES_LIB_TRAILER_LEN].decode("latin-1").lower()
    for needle in EES_LIB_REQUIRED:
        if needle.lower() not in body:
            errors.append("{0} does not define {1}".format(path, needle.rstrip("(")))
    for needle in EES_LIB_REMOVED:
        if needle.lower() in body:
            errors.append("{0} still defines {1}, which was removed".format(path, needle.rstrip("(")))


def check_entry(source_dir, entry, errors):
    """Check one bitness and append a message to errors for every problem."""
    folder = os.path.join(source_dir, entry["folder"])
    for name in entry["files"]:
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            errors.append("missing file: {0}".format(path))
        elif os.path.getsize(path) == 0:
            errors.append("empty file: {0}".format(path))

    ees_library = os.path.join(folder, entry["ees_library"])
    if os.path.isfile(ees_library):
        try:
            check_ees_library(ees_library, errors)
        except (OSError, struct.error) as err:
            errors.append("cannot read the header of {0}: {1}".format(ees_library, err))

    library = os.path.join(folder, entry["library"])
    if not os.path.isfile(library):
        # Already reported above, nothing left to inspect
        return
    try:
        machine = read_pe_machine(library)
    except PE_ERRORS as err:
        errors.append("cannot read the PE header of {0}: {1}".format(library, err))
        return
    if machine != entry["machine"]:
        errors.append(
            "{0} reports machine 0x{1:04x} ({2}) but {3} is required".format(
                library,
                machine,
                MACHINE_NAMES.get(machine, "unknown"),
                MACHINE_NAMES[entry["machine"]],
            )
        )
    else:
        print("ok: {0} is {1}".format(library, MACHINE_NAMES[machine]))

    # EES takes the function name from the file name and looks it up in the
    # export table, so an image without an undecorated COOLPROP_EES is of no
    # use even when the bitness is right.
    try:
        exports = read_pe_exports(library)
    except PE_ERRORS as err:
        errors.append("cannot read the export table of {0}: {1}".format(library, err))
        return
    if EXPORTED_FUNCTION not in exports:
        errors.append(
            "{0} does not export {1}, it exports {2}".format(library, EXPORTED_FUNCTION, ", ".join(exports) or "nothing")
        )
    else:
        print("ok: {0} exports {1}".format(library, EXPORTED_FUNCTION))


def stage(source_dir, stage_dir):
    """Copy the checked files into stage_dir, one subfolder per bitness."""
    for entry in LAYOUT:
        target = os.path.join(stage_dir, entry["folder"])
        if not os.path.isdir(target):
            os.makedirs(target)
        for name in entry["files"]:
            shutil.copyfile(os.path.join(source_dir, entry["folder"], name), os.path.join(target, name))
        print("staged {0} in {1}".format(entry["folder"], target))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source_dir", help="the InnoScript/source folder of the Windows package build")
    parser.add_argument("--stage", dest="stage_dir", default=None, help="copy the checked files into this folder")
    args = parser.parse_args(argv)

    errors = []
    for entry in LAYOUT:
        check_entry(args.source_dir, entry, errors)

    if errors:
        for message in errors:
            print("error: {0}".format(message), file=sys.stderr)
        return 1

    if args.stage_dir is not None:
        stage(args.source_dir, args.stage_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
