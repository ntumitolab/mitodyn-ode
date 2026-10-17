# AGENTS.md

Julia workspace that models mitochondrial dynamics as a ModelingToolkit ODE system.
All figures/notebooks are executable docs. There is no test suite.

## Quick start
- Instantiate deps: `julia --project=@. -e 'using Pkg; Pkg.instantiate()'`
- Run one doc script: `julia --project=@. docs/fig-dm-step.jl`
- Execute a doc as a notebook (what CI does): `NB=docs/fig-dm-step.jl julia --project=@. .github/ci.jl`
  -> executed `.ipynb` in `.cache/docs/`
- Render the site: `cd docs && jupyter-book build --html` (Node/`jupyter-book`; CI installs it globally)

## Commands
- Instantiate deps: `julia --project=@. -e 'using Pkg; Pkg.instantiate()'`
- Run one doc script: `julia --project=@. docs/fig-dm-step.jl`
- Execute a doc as a notebook: `NB=docs/fig-dm-step.jl julia --project=@. .github/ci.jl`
- Render the site: `cd docs && jupyter-book build --html`

No lint, formatter, typecheck, or `runtests.jl`. "Verification" means the doc scripts
run without error. CI runs every `docs/*.jl` and `docs/*.ipynb` in parallel, then renders
the book. Any file added under `docs/` is automatically picked up by the CI matrix.

## Layout
- Root `Project.toml` is both the app env and a Julia 1.13 `[workspace]`
  (`projects = ["Startup"]`). `Startup/` is a local package wired via `[sources]`.
  Run scripts with `julia --project=@.` from the repo root. Compat: `julia = "1.13"`.
- `Startup/` has no `Manifest.toml` (gitignored); it resolves through the workspace.
- `Startup/src/Startup.jl` defines `make_model(; name, caceq=cac_atp(), simplify=true)`,
  returning a (compiled if `simplify=true`) `ModelingToolkit.System`.
  - `name` is required for a direct call; docs use `@named sys = make_model()` so MTK
    injects `name=:sys`.
  - `caceq` is an equation for cytosolic Ca vs ATP/ADP. Default `cac_atp()`;
    pass a custom one (e.g. the doc-local `cac_wave`) for oscillating Ca.
  - `simplify=false` returns the raw system for inspection/tearing.
  - Also defines `cac_atp(; ca_base=0.09μM, ca_act=0.25μM, n=4, katp=25)` (not exported;
    import via `using Startup: cac_atp`).
- `Startup/src/utils.jl`: unit constants are plain numeric scale factors, not Unitful;
  SI base is mM/mV. Constants include `second`, `ms`, `minute`, `Hz`, `kHz`, `metre`,
  `cm`, `cm²`, `mL`, `Liter`, `μL`, `pL`, `mM`, `Molar`, `μM`, `nM`, `Volt`, `mV`,
  `T0=310`, `F`, `R`, `VT`, `iVT`. Also `hil` (two methods), `extract`, `indexof`,
  `aecdev`, `aec2atp`/`aec2adp`. Only `make_model`, `indexof`, `extract`, `exportTIF`
  are exported; the rest are importable via `using Startup: ...`.
  - `aec2atp(a, keq=1)` / `aec2adp(a, keq=1)`; the model uses `kEqAK=0.931`.
- `docs/*.jl` are Literate.jl sources: `#=== ... ===#` is prose, `#---` starts a new
  code cell. Keep these markers.
- `docs/index.md` documents the model equations/parameters; keep it in sync when the
  model changes.
- `docs/myst.yml` configures the Jupyter Book/MyST site; `requirements.txt` pins Python
  deps for CI.

## Gotchas
- `Manifest.toml` is committed and CI reads `julia_version` from it (currently `1.13.2`)
  to pick the toolchain. Do not hand-edit; change it with `Pkg.resolve()`/`Pkg.update()`.
  Local Julia may differ.
- `JSON`, `Literate`, `MarkdownTables`, and `Tables` are now declared in `Project.toml`
  (no longer doc-only). `SHA` is unused in docs and has been dropped. CI no longer
  `Pkg.add`s doc-only packages; all dependencies resolve through `Pkg.instantiate()`.
  If a doc-only package is missing locally, add it to `~/.julia/environments/v1.13`.
- Plotting uses `PythonPlot` (PythonCall/CondaPkg), not PyCall. Locally PythonCall manages
  a CondaPkg env with matplotlib; CI sets `JULIA_CONDAPKG_BACKEND=Null` and gets
  matplotlib from `uv pip install -r requirements.txt` (pinned `3.11.2`). Plotting scripts
  fail without a Python matplotlib. CI also sets `GKSwstype=100` for headless rendering.
- CI (`.github/workflows/ci.yml`) runs only on `main` and PRs to `main`; default branch
  is `main`. Julia dependency bumps arrive as bot PRs (`update-manifest.yml`).
- Generated artifacts (`*.tif`, `*.png`, executed `.ipynb`) are not tracked; CI regenerates
  notebooks and copies them back to `docs/` during render. Do not commit them.
- `exportTIF` writes to the current working directory.

## Environment variables
- `NB`, `NBCACHE`: select input file / cache dir for `.github/ci.jl`.
- `JULIA_CONDAPKG_BACKEND`, `JULIA_CONDAPKG_OFFLINE`: Python backend selection.
- `JULIA_HEAP_SIZE_HINT`, `JULIA_NUM_THREADS`, `OPENBLAS_NUM_THREADS`: performance.
- `GKSwstype`: headless plotting (CI uses `100`).
- `JULIA_CPU_TARGET=sysimage`: set on hosted x64 runners.

## Conventions
- `.editorconfig`: 4-space indent for code, 2 spaces for yml/toml/json/md; trim trailing
  whitespace; final newline. No formatter/linter is configured, so style is by convention.
- Literate markers (`#===`, `#---`) must be preserved.

## When changing the model
1. Update `docs/index.md` (equations and parameter tables) to match the code.
2. Run at least one representative doc locally
   (e.g. `julia --project=@. docs/fig3-1dss.jl`).
3. Confirm all docs still run in CI.

## Troubleshooting
- Plotting fails / no matplotlib: check `python -c "import matplotlib"`; set
  `JULIA_CONDAPKG_BACKEND=Null` to use a system Python, or let CondaPkg manage it locally.
- `UndefVarError` for a util constant/helper: it is not exported; import explicitly,
  e.g. `using Startup: μM, mV, mM, Hz, minute, hil, aec2atp`.
- Doc-only package not found locally: add it to `~/.julia/environments/v1.13`.
