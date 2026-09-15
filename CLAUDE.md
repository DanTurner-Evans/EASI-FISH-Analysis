# Claude Code instructions for data analysis

## Environment

This project uses **uv**. Dependencies live in `pyproject.toml`; `uv.lock` is the source of truth.

- Run everything through `uv run` — e.g. `uv run python src/analyze.py`, `uv run jupyter lab`.
- Never run bare `python`, `pytest`, or `jupyter`; they resolve to the miniconda base env, not this project.
- Add dependencies with `uv add <pkg>` (`uv add --dev <pkg>` for tooling). Do not hand-edit the dependency list and do not run `pip install`.
- Commit `uv.lock` alongside any dependency change.

## Layout

- `workflow/` — Snakefile and analysis scripts
- `data/raw/` — source data, read-only; never write here
- `data/processed/` — derived data, regenerable from `data/raw/`; git-ignored
- `results/` — final tables and figures; tracked in git

## System-wide tools (mamba)

uv owns this project's Python dependencies. **mamba** owns system-wide, non-Python tooling — CLI
binaries and image-processing utilities that are not importable from the analysis code.

- Install with `mamba install -c conda-forge <tool>` into the conda `base` env, so the binary is on `PATH` for every project.
- Keep `base` minimal: only genuinely global tools. A bloated `base` is how conda installs break.
- Never install Python libraries this way. If the analysis code `import`s it, it belongs in `pyproject.toml` via `uv add` — the uv venv cannot see conda-installed Python packages.
- If a required Python package exists only on conda-forge and not PyPI, flag it rather than mixing the two environments; that case needs a decision, not a default.

## Pipeline (Snakemake)

The analysis is a Snakemake workflow: `workflow/Snakefile`, parameters in `config/config.yaml`.

- Run with `uv run snakemake --cores <n>`. Always dry-run (`-n`) first.
- Samples are discovered by globbing `data/raw/*.czi`. The filename IS the metadata:
  `<line A-D><receptor set 1-3>_brain<n><optional a/b>`, e.g. `A2_brain2`, `C1_brain1a`.
  A file that does not match this pattern aborts the run rather than being skipped.
- Adding a sample means dropping a CZI into `data/raw/`. Never edit the Snakefile to add one.
- Driver lines and receptor pairs are declared once in `config/config.yaml`; never hardcode
  a line name, neuron type, or receptor into a script.
- Scripts in `workflow/scripts/` expect Snakemake's injected `snakemake` object — they are not
  runnable directly. Shared helpers go in `workflow/scripts/common.py`.
- To force recomputation use `--forcerun <rule>`, not deletion of `data/processed/`.
