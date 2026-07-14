# utils/ — shared helpers

| file | what it holds |
|---|---|
| `paths.py` | The **single source of truth** for every file path in the repo. All scripts import paths from here — no script hard-codes a location. Change a path once, here. |
| `variant_keys.py` | The variant key used to join tables across stages: `(gene, chrom, pos, ref, alt)`. |
