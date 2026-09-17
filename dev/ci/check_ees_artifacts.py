#!/usr/bin/env python
"""Verify the EES wrapper artifacts of the Windows package.

EES loads external functions by bitness: the 32-bit program reads
COOLPROP_EES.dlf from its Userlib folder, the 64-bit program (EES64.exe)
reads COOLPROP_EES.dlf64 from Userlib64.  A library of the wrong bitness is
silently ignored by EES, so a build that produces two files with the same
machine type would ship a broken 64-bit wrapper without anybody noticing.

This script checks that both files exist and that their PE headers really
report the expected machine type.  With --stage it also copies the checked
files into a clean folder, which is what the release workflow uploads.

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

# One entry per bitness: folder, the library file and the expected machine
LAYOUT = [
    {
        "folder": "EES",
        "library": "COOLPROP_EES.dlf",
        "machine": IMAGE_FILE_MACHINE_I386,
        "files": ["COOLPROP_EES.dlf", "CoolProp.LIB", "CoolProp.htm", "CoolProp_EES_Sample.EES"],
    },
    {
        "folder": "EES64",
        "library": "COOLPROP_EES.dlf64",
        "machine": IMAGE_FILE_MACHINE_AMD64,
        "files": ["COOLPROP_EES.dlf64", "CoolProp.LIB64", "CoolProp.htm", "CoolProp_EES_Sample.EES"],
    },
]


def read_pe_machine(path):
    """Return the machine type stored in the PE header of path.

    Raises ValueError if the file is not a PE image.
    """
    with open(path, "rb") as handle:
        header = handle.read(0x1000)
    if len(header) < 0x40 or header[0:2] != b"MZ":
        raise ValueError("no MZ signature, this is not a Windows library")
    # The offset of the PE header is stored at 0x3C in the DOS header
    (pe_offset,) = struct.unpack_from("<I", header, 0x3C)
    if pe_offset + 6 > len(header) or header[pe_offset : pe_offset + 4] != b"PE\0\0":
        raise ValueError("no PE signature at offset {0}".format(pe_offset))
    (machine,) = struct.unpack_from("<H", header, pe_offset + 4)
    return machine


def check_entry(source_dir, entry, errors):
    """Check one bitness and append a message to errors for every problem."""
    folder = os.path.join(source_dir, entry["folder"])
    for name in entry["files"]:
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            errors.append("missing file: {0}".format(path))
        elif os.path.getsize(path) == 0:
            errors.append("empty file: {0}".format(path))

    library = os.path.join(folder, entry["library"])
    if not os.path.isfile(library):
        # Already reported above, nothing left to inspect
        return
    try:
        machine = read_pe_machine(library)
    except (ValueError, OSError) as err:
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
