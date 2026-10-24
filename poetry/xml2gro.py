#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from shape_model.xml_io import XML, backimage
from .read_itp_sel import read_line

def _load_atoms(itp_files):
    mol_itps = []
    for itp in itp_files:
        atoms, _, _, _, _ = read_line(itp)
        mol_itps.append(atoms)
    return mol_itps


def xml2gro(xml, itp_files, counts, out=None, title="mono_all", use_image=True, shift=True):
    """Convert GALAMOST/pygamd XML to GROMACS GRO."""
    itp_files = list(itp_files)
    counts = [int(i) for i in counts]
    if len(itp_files) != len(counts):
        raise ValueError("itp_files and counts must have the same length")

    mol_itps = _load_atoms(itp_files)
    info = XML(xml).return_info()
    pos = [[float(x) for x in p[:3]] for p in info["position"]]
    image = info.get("image", [])
    box = [float(x) for x in info["box"][:3]]

    if use_image and image:
        for i in range(len(pos)):
            backimage(pos[i], image[i], box)

    if shift:
        for p in pos:
            p[0] += box[0] / 2.0
            p[1] += box[1] / 2.0
            p[2] += box[2] / 2.0

    expected = sum(len(atoms) * n for atoms, n in zip(mol_itps, counts))
    if expected != len(pos):
        raise ValueError(f"atom number mismatch: XML={len(pos)}, ITP/counts={expected}")

    out = out or str(Path(xml).with_suffix(".gro"))
    Path(out).parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w") as f:
        f.write(str(title) + "\n")
        f.write(str(len(pos)) + "\n")
        k = 0
        mol_n = 0
        for xi, atoms in enumerate(mol_itps):
            for i in range(counts[xi]):
                for atom in atoms:
                    f.write(
                        f"{(mol_n + i + 1) % 100000:5d}"
                        f"{atom[3]:<5.5}"
                        f"{atom[4]:>5.5}"
                        f"{(k + 1) % 100000:5d}"
                        f"{pos[k][0]:8.3f}{pos[k][1]:8.3f}{pos[k][2]:8.3f}\n"
                    )
                    k += 1
            mol_n += counts[xi]
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")

    return out


def main():
    parser = argparse.ArgumentParser(description="Convert GALAMOST/pygamd XML to GROMACS GRO")
    parser.add_argument("--xml", required=True)
    parser.add_argument("--itp", nargs="+", required=True)
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--title", default="mono_all")
    parser.add_argument("--no-image", action="store_true")
    parser.add_argument("--no-shift", action="store_true")
    args = parser.parse_args()

    out = xml2gro(
        xml=args.xml,
        itp_files=args.itp,
        counts=args.counts,
        out=args.out,
        title=args.title,
        use_image=not args.no_image,
        shift=not args.no_shift,
    )
    print(out, "success!")


if __name__ == "__main__":
    main()
