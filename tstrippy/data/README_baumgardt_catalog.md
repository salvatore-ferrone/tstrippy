# Baumgardt Catalog Update Workflow

This package reads Baumgardt MWGC data directly from two raw ASCII exports:

- `mwgcs-kinematics`
- `mwgcs-structural`

## Update steps

1. Download/export the latest ASCII tables from the source website.
2. Replace the contents of `mwgcs-kinematics` with the kinematics table.
3. Replace the contents of `mwgcs-structural` with the structural table.
4. Run tests:
   - `pytest tests/test_baumgardt_ascii_catalog.py -q`

No manual FITS editing is required.
