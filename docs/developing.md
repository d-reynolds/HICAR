# Developing HICAR

Contributions are welcome. HICAR is structured so that common additions — a new
output variable, a new physics option — are straightforward, and the maintainers
are happy to help with larger changes. For a map of the source tree, start with
the [Code overview](code_overview.md).


## Contribution workflow

- Make changes and open pull requests against the **`develop`** branch.
- Every PR is gated by the automated CI suite (build, unit/invariant tests,
  smoke runs, decomposition/restart reproducibility, and bit-for-bit
  regression). See the [CI/CD pipeline](ci_cd_pipeline.md) for what each gate
  checks and [Testing](testing.md) for how to reproduce a lane locally before
  pushing.
- For substantial additions, please get in touch first so we can help fit the
  change into the model's structure.

## Adding a new output variable

HICAR's variable handling is **metadata-driven**: once a variable is registered
and described, allocation and output happen automatically — there is no
per-variable code to add in the domain or output objects. To add a new variable, three things are needed:

1. **Give it an index.** Add an entry for the variable in the `kVARS` structure
   in `src/constants/icar_constants.F90`. Every model variable is referred to
   throughout the code by its `kVARS%name` index.

2. **Describe it.** Add a branch for the new `kVARS` index to the `get_varmeta`
   function in `src/io/default_output_metadata.F90`. This is where you set the
   output `name`, the `dimensions` (e.g. `three_d_t_dimensions`), the `units`
   and other netCDF `attributes`, and the valid `minval`/`maxval` range. The
   dimensions you set here are what drives allocation, so this step is required —
   a `kVARS` entry with no metadata branch returns an empty name and is silently
   skipped.

3. **Request it.** Have the relevant physics module ask for the variable in its
   `*_var_request` routine (e.g. `adv_std_var_request` in
   `src/physics/advect.F90`, `lsm_var_request` in `lsm_driver.F90`). Use
   `options%alloc_vars([...])` to have it allocated and `options%restart_vars([...])`
   to have it carried through restarts and written to output:

   ```fortran
   call options%alloc_vars(   [kVARS%my_new_var, ...] )
   call options%restart_vars( [kVARS%my_new_var, ...] )
   ```

That's all. The domain object's `create_variables` loops over the requested
variables, reads each one's metadata, and allocates the correctly-shaped array;
the output path writes whatever is in the output set. (`vars_to_allocate`
controls which `kVARS` are allocated; `vars_for_output` controls which are
written.)

If the variable you want is *already* registered, simply request it (step 3) in
the module that needs it — or, for output, it is exposed through the namelist's
output-variable selection.

## Domain decomposition conventions

When writing physics or I/O code, note the index conventions:

- `ims:ime` — memory bounds (include halos)
- `its:ite` — tile/interior bounds (this rank's computed cells)
- `ids:ide` — global domain bounds

Halo (ghost-cell) data is exchanged via the domain's halo methods. Operations
that need up-to-date neighbor data must run *after* the relevant halo exchange.
Be especially careful with anything whose result depends on the local subdomain
extent (local reductions used as global parameters, search windows sized from
local bounds) — these break bit-for-bit reproducibility across different MPI rank
counts.

## Building the documentation

These docs are built with [MkDocs](https://www.mkdocs.org/). To preview them
locally:

```bash
pip3 install mkdocs
mkdocs serve          # or: python3 -m mkdocs serve
```

Then open <http://127.0.0.1:8000> in a browser. The site layout (page order and
titles) is defined in `mkdocs.yml`. Use `mkdocs build --strict` to catch broken
intra-doc links before committing.
