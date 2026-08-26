# Sanger Seek

Sanger Seek is a local macOS desktop application for reviewing Sanger chromatograms,
reconciling forward/reverse evidence, and auditing candidate variants against the raw trace.
Sequence data stays on the Mac.

> Sanger Seek is research-review software. Its calling and QC rules have not been formally
> validated for clinical diagnostic use.

## Run from source

Python 3.11–3.13 is recommended.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[fast]'
.venv/bin/sanger-seek --demo
```

You can also open a folder or files directly:

```bash
.venv/bin/sanger-seek /path/to/sanger-folder
.venv/bin/sanger-seek --project review.sanger-seek.json
```

The application accepts `.ab1`/`.abi`, `.seq`, FASTA, and GenBank files. Use **Load Reference**
for the supplied reference, **Load WT Control** for one or more forward/reverse control contigs,
and **Load Assessed Samples** for the samples under review. Reads are aligned independently to
the reference and sample variants are labeled according to whether the same allele is present
in the WT control.

An explicitly selected `.seq` file can be the reference. GenBank records saved with a `.seq`
suffix are detected from their `LOCUS` header. Plain `.seq` references are sequence-only;
GenBank references retain CDS annotations and enable coding/protein consequences.

Chromatograms support trackpad/mouse-wheel zoom, double-click-to-zoom, visible `+`, `−`, and
**Fit** controls, and synchronized keyboard shortcuts: `⌘+`, `⌘−`, and `⌘0`. Up to four
chromatograms are shown at once—normally WT forward/reverse and assessed-sample
forward/reverse. Zooming or horizontally dragging any trace synchronizes every visible trace
through reference coordinates. The base pair at the center is the shared anchor, even when
the underlying AB1 sample positions and trace lengths differ.

Use **Pair / Assign Reads…** (`⌘P`) to override automatic grouping: give matching reads the
same sample name and explicitly choose Forward or Reverse when filenames or automatic
orientation are ambiguous. The assignments are saved in the project.

## macOS application and DMG

On macOS, run:

```bash
./scripts/build_macos.sh
```

The build creates:

- `dist/Sanger Seek.app`
- `dist/Sanger-Seek-0.1.0.dmg`

The default build is ad-hoc signed and is suitable for local testing. For distribution, provide
a Developer ID Application identity and an optional notarytool keychain profile:

```bash
APPLE_SIGN_IDENTITY="Developer ID Application: Example, Inc. (TEAMID)" \
APPLE_NOTARY_PROFILE="sanger-seek-notary" \
./scripts/build_macos.sh
```

Create the notary profile once with `xcrun notarytool store-credentials`. Apple notarization
requires an Apple Developer account and cannot be completed using an ad-hoc signature.

To select a particular compatible Python installation, set `PYTHON_BIN`, for example:

```bash
PYTHON_BIN=/opt/homebrew/bin/python3.13 ./scripts/build_macos.sh
```

## Tests

```bash
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest -q
```
