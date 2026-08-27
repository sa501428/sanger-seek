# Sanger Seek — reference-first case review

Sanger Seek is a local PySide6 desktop application for reviewing bidirectional Sanger
sequencing against a user-assigned region/gene reference. The app uses neutral **cases**:
control, sample, WT, standard, and other biological labels do not change analysis behavior.

## Inputs and ownership

```text
Reference (.seq; FASTA/GenBank also supported)
  └── Case
      ├── forward.ab1              required for forward chromatogram
      ├── forward.seq              optional companion calls
      ├── forward.phd.1            optional Phred calls/qualities/locations
      ├── reverse.ab1              required for reverse chromatogram
      ├── reverse.seq              optional companion calls
      └── reverse.phd.1            optional Phred calls/qualities/locations
```

The selected AB1 files are authoritative for the raw traces, embedded calls, peak locations,
and qualities. Optional companions are paired by logical read stem (including
`read.ab1.phd.1`) and compared with AB1 calls. Companion disagreements are review flags.

## Pipeline

For each read: parse inputs, measure peaks/noise, quality-trim, choose or honor orientation,
align to the reference, and map AB1 sample positions into reference coordinates. Candidate
SNVs, insertions, deletions, and mixed peaks are then reconciled across the case's forward and
reverse evidence and annotated against reference CDS features when available.

## Review UI

The main window contains a case list, reference/read/variant summary, reference-centered
alignment strip, synchronized forward and reverse chromatograms, and a clickable variant
table. Selecting a variant centers every view on its reference coordinate. Trace dragging,
wheel zoom, toolbar zoom, and Fit operate in shared reference coordinates even though the two
AB1 files have different raw sample axes.

Review flags cover missing strands, low mean Q, heavy trimming, frequent secondary peaks,
weak alignment, AB1/SEQ/PHD disagreement, noisy variant evidence, and forward/reverse
discordance. They prompt human inspection and are not clinical pass/fail decisions.

All processing is local and background analysis uses Qt's thread pool.
