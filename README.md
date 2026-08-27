# Sanger Seek

Sanger Seek is a local macOS desktop application for reviewing forward and reverse Sanger
chromatograms against a user-assigned reference, reconciling strand evidence, and auditing
candidate variants against the raw traces. Sequence data stays on the Mac.

> Sanger Seek is research-review software. Its calling and QC rules have not been formally
> validated for clinical diagnostic use.

## Review workflow

1. Choose **Load Reference…** and assign the `.seq` reference for the region or gene of
   interest. FASTA and GenBank references remain supported; GenBank CDS annotations enable
   coding and protein consequence labels.
2. Choose **New Case…**, name the case without assigning it a special biological role, then
   select all files belonging to it.
3. Include the forward and reverse `.ab1` chromatograms. Same-stem `.seq`, `.phd`, and
   numbered `.phd.1` files are optional companions.
4. Review reference-aligned differences, strand support, chromatograms, and QC/artifact flags.

Cases are neutral: a case can represent a sample, control, WT, standard, or any other material.
The app does not privilege one case as a control or subtract one case from another.

The AB1 remains authoritative for the raw electropherogram, embedded calls, qualities, and
peak positions. Optional SEQ and PHD calls are compared with the AB1 calls and disagreements
are surfaced for review; they do not replace the trace. Common Phred naming such as
`read_F.ab1.phd.1` is recognized and paired with `read_F.ab1`.

Forward and reverse reads are aligned independently to the reference and displayed in the
same reference orientation. Dragging or zooming either chromatogram synchronizes the other
through reference coordinates. Use the mouse/trackpad, visible `+`, `−`, and **Fit** controls,
or `⌘+`, `⌘−`, and `⌘0`. Clicking a variant row centers the alignment and both chromatograms at
that reference position.

Review flags include missing forward/reverse traces, low mean quality, heavy trimming,
frequent secondary peaks, weak reference alignment, companion-call disagreement, noisy
variant evidence, and strand disagreement. These are prompts for visual review, not automatic
pass/fail decisions.

Use **Assign Reads to Cases…** (`⌘P`) to correct case membership or explicitly set a read as
Forward or Reverse when filenames and automatic alignment are ambiguous.

## Run from source

Python 3.11–3.13 is recommended.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[fast]'
.venv/bin/sanger-seek --demo
```

You can also import a folder containing a reference and conventionally named case files:

```bash
.venv/bin/sanger-seek /path/to/sanger-folder
.venv/bin/sanger-seek --project review.sanger-seek.json
```

## macOS application and DMG

On macOS, run:

```bash
./scripts/build_macos.sh
```

The build creates:

- `dist/Sanger Seek.app`
- `dist/Sanger-Seek-0.1.0.dmg`

The default build is ad-hoc signed and suitable for local testing. For distribution, provide
a Developer ID Application identity and an optional notarytool keychain profile:

```bash
APPLE_SIGN_IDENTITY="Developer ID Application: Example, Inc. (TEAMID)" \
APPLE_NOTARY_PROFILE="sanger-seek-notary" \
./scripts/build_macos.sh
```

To select a compatible Python installation, set `PYTHON_BIN`, for example:

```bash
PYTHON_BIN=/opt/homebrew/bin/python3.13 ./scripts/build_macos.sh
```

## Tests

```bash
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest -q
```
