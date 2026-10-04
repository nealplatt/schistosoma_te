# Agent instructions

## Project scope

This comparative bioinformatics project studies transposable elements (TEs) across
the genus *Schistosoma*. Its long-term goal is to discover, curate, and classify a
pan-*Schistosoma* TE library and use it to study genome TE content, lineage-specific
TE families and expansions, and TE evolutionary history.

Development occurs locally in WSL. The eventual workflow will use Snakemake,
Conda, Python, standard TE-analysis software, Git, and Slurm HPC. Do not treat
undecided tools, methods, or biological assumptions as established requirements.

## Scientific decisions and data

- Keep scientific decisions explicit and inspectable. Record assumptions,
  rationale, parameters, and supporting evidence in reviewable files.
- Do not silently make biological classification decisions. Ask for direction
  when a task requires an undecided classification or curation choice.
- Identify unresolved scientific or methodological decisions explicitly; do not
  invent defaults. Continue independent work when possible, but stop before
  implementing anything that would implicitly resolve them.
- Do not discard candidate TE families or other biological data without explicit
  instructions. Preserve intermediate evidence used for classification and curation.
- Do not overwrite raw input data. Write derived outputs separately and retain
  traceability to their inputs.
- Retain analysis provenance: input assembly/database versions, relevant software
  versions, parameters, and relationships between source and derived files.
- Keep generated data separate from source code; avoid committing large generated
  datasets to Git unless explicitly instructed.

## Implementation and execution

- Prefer reproducible scripts and workflow rules over manual commands. Capture
  inputs, outputs, parameters, and software versions where appropriate.
- Before substantial architectural or methodological changes, propose a plan
  for review before implementing it.
- Make small, focused, Git-friendly changes that are easy to review. Avoid
  unrelated edits and validate changes with checks appropriate to their scope.
- Successful execution alone does not establish biological or methodological
  validity. Distinguish software validation from biological interpretation.
- Use Python with lowercase variable names and `pathlib` for paths where appropriate.
- Use TSV for ordinary tabular bioinformatics outputs unless another format is
  clearly preferable.
- Make software environments reproducible with Conda environment files as the
  workflow develops; document dependencies rather than relying on local state.
- Use local WSL for development and small validation runs. Do not run large
  compute jobs locally unless explicitly requested; target Slurm HPC for them.
- Do not submit Slurm jobs or download large datasets unless explicitly requested.
- Do not git commit, push, merge, or otherwise modify remote repositories unless
  explicitly requested. Inspection commands and small local validation tests
  are acceptable.
