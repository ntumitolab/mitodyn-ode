# AGENTS.md

Julia workspace that models mitochondrial dynamics as a ModelingToolkit ODE system.
All docs are executable Literate.jl scripts; there is no test suite.

## Commands

- Instantiate deps: `julia --project=@. -e 'using Pkg; Pkg.instantiate()'`
- Run one doc directly (no notebook/cache, writes TIFFs to cwd):
  `julia --project=@. docs/fig3-1dss.jl`
- Execute docs as notebooks like CI does:
  `NBCACHE=.cache julia --project=@. -p 4 .github/ci.jl`
  -> executed `.ipynb` under `.cache/docs/`
- Limit CI to a directory (not a single file): `DOCDIR=docs NBCACHE=.cache julia --project=@. .github/ci.jl`
- Render the site (after notebooks are copied into `docs/`): `cd docs && jupyter-book build --html`
  (CI installs `jupyter-book` globally via npm.)

No lint, formatter, typecheck, or `runtests.jl`. "Verification" means the doc scripts run
without error. CI runs every `docs/*.jl` and `docs/*.ipynb` in parallel, then renders.

## Layout

- Root `Project.toml` is both the app env and a Julia 1.13 `[workspace]`
  (`projects = ["Startup"]`). `Startup/` is a local package wired via `[sources]`.
  Run everything with `julia --project=@.` from the repo root. Compat: `julia = "1.13"`.
- `Startup/` has no `Manifest.toml` (gitignored); it resolves through the workspace.
- `Startup/src/Startup.jl`: `make_model(; name, caceq=cac_atp(), simplify=true)`.
  `name` is required for a direct call; docs use `@named sys = make_model()` so MTK
  supplies it. `simplify=true` runs `mtkcompile`; `false` returns the raw `System`.
  `caceq` is the cytosolic Ca equation (default `cac_atp()`); pass a custom one
  (e.g. `cac_wave` in `docs/ca-osci.jl`) for oscillating Ca.
  `cac_atp(; ca_base=0.09μM, ca_act=0.25μM, n=4, katp=25)` is not exported.
- `Startup/src/utils.jl`: units are plain numeric scale factors, not Unitful; SI base is
  mM/mV. Includes `second`, `ms`, `minute`, `Hz`, `kHz`, `metre`, `cm`, `cm²`, `mL`,
  `Liter`, `μL`, `pL`, `mM`, `Molar`, `μM`, `nM`, `Volt`, `mV`, `T0=310`, `F`, `R`,
  `VT`, `iVT`, plus `hil` (two methods), `extract`, `indexof`, `aecdev`,
  `aec2atp`/`aec2adp` (model uses `kEqAK=0.931`). Only `make_model`, `indexof`,
  `extract`, and `exportTIF` are exported.
- `docs/*.jl` are Literate sources: `#=== ... ===#` is prose, `#---` starts a code cell.
  Preserve these markers.
- `docs/index.md` documents model equations/parameters; keep it in sync with the model.
  `docs/packages.md` lists root dependencies with versions.
- `docs/myst.yml` configures the MyST/Jupyter Book site; `*.ipynb` in `docs/` is the toc
  pattern, so executed notebooks copied there are published automatically.
- `Manifest.toml` is committed and CI reads `julia_version` from it (currently `1.13.2`)
  to pick the toolchain. Do not hand-edit it; use `Pkg.resolve()`/`Pkg.update()`.

## Gotchas

- `.github/ci.jl` reads `DOCDIR` (default `docs`) and `NBCACHE` (default `.cache`).
  There is no per-file selector; to test one doc, run it directly instead.
- CI hashes each doc and skips execution on an unchanged cache hit. Delete
  `.cache/docs/<name>.sha` (and the cached notebook) to force a rerun.
- `ci.jl` converts any `docs/*.ipynb` to a `.jl` with Literate markers and deletes the
  original `.ipynb`; prefer adding/editing `docs/*.jl`.
- Plotting uses `PythonPlot` (PythonCall/CondaPkg), not PyCall. Locally CondaPkg manages
  a matplotlib env (`.CondaPkg/`, gitignored); without a Python matplotlib the plotting
  docs fail. CI sets `JULIA_CONDAPKG_BACKEND=Null` and installs `requirements.txt`
  (matplotlib) with `uv`, plus `GKSwstype=100` for headless rendering.
- `exportTIF` writes TIFFs to the current working directory. Generated `*.tif`,
  `*.png`, and executed `.ipynb` files are gitignored; do not commit them.
- CI (`.github/workflows/ci.yml`) runs on `main` and PRs to `main`; default branch is
  `main`. Dependency bumps arrive as bot PRs from `update-manifest.yml` / Dependabot.
- `SHA` is still a root dependency because `.github/ci.jl` hashes notebooks; do not drop
  it just because no doc imports it.
- If a package appears missing locally, instantiate first; fall back to adding doc-only
  packages to `~/.julia/environments/v1.13`.

## Environment variables

- `DOCDIR`, `NBCACHE`: input directory / cache directory for `.github/ci.jl`.
- `JULIA_CONDAPKG_BACKEND`, `JULIA_CONDAPKG_OFFLINE`: Python backend selection.
- `JULIA_HEAP_SIZE_HINT`, `JULIA_NUM_THREADS`, `OPENBLAS_NUM_THREADS`: performance.
- `GKSwstype`: headless plotting (CI uses `100`).
- `JULIA_CPU_TARGET=sysimage`: set on hosted x64 runners only.

## Conventions

- `.editorconfig`: 4-space indent for code, 2 spaces for yml/toml/json/md; trim trailing
  whitespace; final newline. No formatter/linter is configured.
- Keep Literate markers (`#===`, `#---`) and update `docs/index.md` when equations or
  parameters change.
- Run at least one representative doc locally (e.g. `julia --project=@. docs/fig3-1dss.jl`)
  before claiming a model change works; real verification is all docs passing in CI.
