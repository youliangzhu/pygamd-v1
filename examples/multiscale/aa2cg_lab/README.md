# Lightweight AA-to-CG example

This self-contained example maps 100 AC10 molecules (9200 AA atoms) to a
20-bead-per-molecule CG model.

- `aa2cg_nordf.json` runs the mapping and GALAMOST export without reading a
  trajectory.
- `aa2cg_rdf.json` also computes one `BBE-CO` RDF from
  `md_nvt_sample.xtc`.
- `md_nvt_sample.xtc` contains three representative frames extracted from
  source frame indices 0, 500, and 1000 (0, 5000, and 10000 ps).

From the repository root, install the RDF dependencies and validate both
configurations:

```console
python -m pip install ".[rdf]"
python -m multiscale aa2cg validate-config examples/multiscale/aa2cg_lab/aa2cg_nordf.json --json
python -m multiscale aa2cg validate-config examples/multiscale/aa2cg_lab/aa2cg_rdf.json --json
```

Run either workflow with:

```console
python -m multiscale aa2cg run examples/multiscale/aa2cg_lab/aa2cg_nordf.json --json
python -m multiscale aa2cg run examples/multiscale/aa2cg_lab/aa2cg_rdf.json --json
```

Generated files are written below `outputs/` in this example directory.
