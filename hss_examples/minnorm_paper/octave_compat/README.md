# octave_compat (only for GNU Octave; MATLAB users can ignore this folder)

The suite runs in GNU Octave 8.4 or later on the repo's `@hss` code, unmodified except for three
syntax-level incompatibilities that the files here work around:

1. `dictionary.m`: a minimal stand-in for MATLAB's `dictionary` (used by `hss_matvec`),
   built on `containers.Map`.
2. `make_octave_copy.py`: writes an Octave-parsable *copy* of `@hss` and `+hssutil` into a
   folder of your choice. The copy replaces the `arguments` blocks in `hss_constructor.m` and
   `inter_decompv3.m` by explicit parsing, uses `qr(A,0)` for pivoted economy QR, and swaps
   the `sparsesign` MEX for `sparsesign_octave.m`. The solver files
   (`hss_ulvminnormsolve.m`, `hss_matvec.m`, ...) are copied byte for byte.
3. Name=value calls: Octave accepts `f(a, name=value)` syntactically but drops the names, so
   the suite always uses `'name', value`.

```bash
python3 make_octave_copy.py /path/to/structmats /tmp/hss_octave
octave --eval "addpath('octave_compat'); addpath('/tmp/hss_octave'); cd hss_examples/minnorm_paper; mp_selftest"
```
Run from the repo root. Octave timings are not representative; use these runs for
correctness only.
