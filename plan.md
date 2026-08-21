# Sanger Seek — desktop Sanger trace review app

A local desktop application for reviewing Sanger sequencing results: load `.ab1` chromatograms
and their exported `.seq` base calls, align reads against a reference, reconcile forward/reverse
evidence, and review candidate variants directly against the raw traces. The primary review
interface is a modernized, decluttered version of classic trace-review software: **keep the
information density available, but use progressive disclosure.** The default screen should
immediately answer *"where are the differences, do forward/reverse support them, and what do
the peaks actually look like?"* Everything else sits behind expandable QC/details panels.

Nothing ever leaves the machine.

## 1. Stack and architecture (desktop Python)

* **Python 3.11+** (developed against 3.13)
* **PySide6 / Qt** — desktop UI
* **pyqtgraph** — fast interactive chromatogram plotting, zoom/pan/cursors
* **Biopython** — `.ab1` (ABIF) parsing, FASTA, GenBank, sequence manipulation, translation
* **NumPy** — trace arrays / signal calculations
* **edlib** — fast read↔reference alignment (with a pure-Python/NumPy fallback aligner so the
  app still works if the compiled extension is unavailable)
* Plain Python code for:
  * `.seq` parsing
  * file pairing / forward-reverse matching
  * QC + quality trimming
  * consensus / strand reconciliation
  * SNV/indel/mixed-base detection
  * codon/protein consequence calculation
  * trace-quality metrics

Analysis runs on background threads (`QThreadPool`) so the UI never blocks; results are
delivered back to the UI via Qt signals.

```
Open Project
   │
   ├── reference.fasta / .gb
   │
   └── Sample
       ├── forward.ab1
       ├── forward.seq
       ├── reverse.ab1
       └── reverse.seq
                │
                ▼
         Parse + QC reads
                │
                ▼
       Align against reference
                │
                ▼
       Reconcile FWD + REV
                │
                ▼
       Detect / annotate variants
                │
                ▼
    ┌─────────────────────────────┐
    │ Modern chromatogram viewer  │
    │ Sequence alignment          │
    │ Variant navigator           │
    │ QC / noise inspection       │
    └─────────────────────────────┘
```

### Why desktop

The desktop removes the browser security model's friction. The app can directly open a folder
containing:

```
Patient001_F.ab1
Patient001_F.seq
Patient001_R.ab1
Patient001_R.seq
Patient002_F.ab1
...
```

and automatically pair the files, rather than requiring repeated file-selection operations.
It also gets straightforward drag-and-drop, folder scanning, local save/open projects, native
file dialogs, CSV exports, persistent settings, larger datasets, and real threading.

## 2. File handling

**`.ab1` — raw chromatogram source.** Parsed via Biopython's ABIF reader to extract:

* A/C/G/T trace arrays (`DATA9–12`, analyzed traces; channel order from `FWO_`)
* called bases (`PBAS`)
* peak/sample positions (`PLOC`)
* base quality values (`PCON`)
* run/instrument metadata where useful (`SMPL`, `MCHN`, `DySN`, run dates)

The `.ab1` is the authoritative source for the **actual electropherogram**, peak quality,
mixed peaks, noise, and manual visual review.

**`.seq` — exported called sequence.** Parsed as plain text (optionally FASTA-headered):
normalize whitespace/headers, associate with the corresponding `.ab1` (same basename), use as
the imported base-call sequence, and compare against the AB1's embedded calls — discrepancies
are counted and inspectable, never silently resolved. `.ab1` alone also works, since AB1
contains base calls itself.

```
sample_F.ab1   -> trace + peaks + qualities + embedded base calls
sample_F.seq   -> exported/called nucleotide sequence
```

They are two representations of the **same read**, not independent sequencing evidence.

**Reference.** FASTA for sequence-only; GenBank for CDS/exon/gene annotation and amino-acid
consequences (Biopython `SeqRecord.features`, including `join`/`complement` locations and
`codon_start`).

**Pairing heuristics.** Group files into samples by basename after stripping direction tokens
(`_F`, `_R`, `-FWD`, `_REV`, `_forward`, …). The token is only an orientation *hint*; final
orientation is decided by alignment. Users can re-assign reads between samples in the UI.

## 3. Analysis pipeline

For each sample/read:

```
Load files
   ↓
Parse AB1 / SEQ
   ↓
QC + low-quality trimming (Mott algorithm on PCON)
   ↓
Determine F/R orientation (best alignment of read vs revcomp)
   ↓
Reverse-complement when required
   ↓
Align read to reference (edlib infix alignment → op path)
   ↓
Map trace peaks → reference coordinates
   ↓
Call candidate differences
   ↓
Reconcile forward + reverse evidence
   ↓
Annotate variants
```

Detect:

* SNVs
* insertions
* deletions
* ambiguous/mixed bases (IUPAC calls and/or secondary-peak ratio)
* forward/reverse discordance

With an annotated coding reference, additionally calculate: synonymous / missense / nonsense /
frameshift / in-frame indel, codon change, and protein change, e.g. `c.944C>T / p.Thr315Ile`.

## 4. Main review interface

Single main window, top-to-bottom:

```
┌ Sample / Reference / QC summary ───────────────────────────────┐
│ Sample 014   Gene XYZ   2 reads   3 candidate variants        │
└────────────────────────────────────────────────────────────────┘

┌ Sequence alignment / amino-acid track ────────────────────────┐
│ Reference   ... C T A T C A C T G ...                         │
│ Forward     ... C T A T C A T T G ...          ↑ variant      │
│ Reverse     ... C T A T C A T T G ...                         │
│ Protein           I   T   E   F                                │
└────────────────────────────────────────────────────────────────┘

┌ Forward chromatogram ──────────────────────────────────────────┐
│       /\       /\  /\                                          │
│ ___A_/  \__G__/  \/  \_T____                                  │
│               ↑ c.944C>T                                      │
└────────────────────────────────────────────────────────────────┘

┌ Reverse chromatogram ──────────────────────────────────────────┐
│          /\     /\                                             │
│ ___G____/  \_A_/  \___                                        │
│               ↑ same genomic position                          │
└────────────────────────────────────────────────────────────────┘

┌ Variants ──────────────────────────────────────────────────────┐
│ c.944C>T │ p.Thr315Ile │ F ✓ │ R ✓ │ High confidence │ View │
└────────────────────────────────────────────────────────────────┘
```

Plus a samples list dock on the left (one row per sample with read counts and QC chips).

## 5. Chromatogram UX

Eliminate dense diagnostic-looking clutter. The trace viewer (pyqtgraph) has:

* smooth A/C/G/T traces
* synchronized base calls above peaks (rendered only for the visible range)
* reference-coordinate ruler (custom axis mapping trace samples → reference positions)
* quality shading
* selected-position crosshair
* variant flags
* zoom/pan, drag navigation
* keyboard left/right between bases; next/previous variant shortcuts

When a user clicks `c.944C>T`, every view instantly centers on that nucleotide. Forward and
reverse traces are displayed in the **same reference orientation** (reverse reads are drawn
mirrored with complemented channel assignment), so users never mentally reverse a strand.

## 6. Noise / mixed-peak visualization

For every called base calculate/display:

* primary peak height
* secondary peak height
* secondary/primary ratio
* local baseline/noise
* AB1 quality score
* distance from sequencing read ends

For candidate heterogeneity/mixed peaks, show a small indicator such as:

```
C 72%
T 28%
```

and let the user inspect the original trace rather than hiding this behind an automated call.

## 7. Variant review table

A polished fixed panel (QTableView):

| Variant    | Protein     | Fwd | Rev | Trace | Confidence |
| ---------- | ----------- | --: | --: | ----- | ---------- |
| c.944C>T   | p.Thr315Ile |   ✓ |   ✓ | Clean | High       |
| c.1012delA | p.X338fs    |   ✓ |   ✓ | Clean | High       |
| c.1173G>A  | synonymous  |   ✓ |   ? | Noisy | Review     |

Clicking any row jumps every view to that trace position.

Filters: all differences / coding only / missense+nonsense / indels / mixed peaks /
strand disagreement / low confidence. Export the (filtered) table to CSV.

## 8. Forward/reverse reconciliation

A first-class object, not two unrelated traces. For each reference position:

```
Reference: C
Forward:   T   Q=34
Reverse:   T   Q=37

Consensus: T
Status:    Supported on both strands
```

Versus:

```
Reference: C
Forward:   T   Q=31
Reverse:   C   Q=36

Consensus: ?
Status:    Strand disagreement — review
```

This drives the mutation flags and the confidence column.

## 9. Multi-read / multi-amplicon support

Data model:

```
Project
 └─ Samples[]
     ├─ Reference
     ├─ Reads[]
     │   ├─ AB1
     │   ├─ SEQ
     │   ├─ orientation
     │   ├─ trace
     │   ├─ quality
     │   └─ alignment
     ├─ Consensus
     └─ Variants[]
```

If multiple reads/amplicons exist, each is mapped independently onto the same reference and
merged into one reference-coordinate-based variant list.

## 10. Repository layout

```
sanger-seek/
  pyproject.toml
  sanger_seek/
    __main__.py          # python -m sanger_seek
    app.py               # QApplication bootstrap, CLI flags (--demo, --screenshot)
    core/
      model.py           # dataclasses: Project, Sample, Read, Variant, Reference
      abif.py            # AB1 → TraceData (traces, calls, ploc, quals, metadata)
      seqfile.py         # .seq parsing
      reference.py       # FASTA/GenBank loading (seq + CDS features)
      pairing.py         # folder scan, sample/read grouping, orientation hints
      trim.py            # Mott quality trimming
      align.py           # edlib alignment + read↔ref coordinate maps (NW fallback)
      peaks.py           # per-base peak metrics
      variants.py        # per-read candidate calls
      consensus.py       # strand reconciliation → sample variants
      consequence.py     # CDS/codon/protein annotation (c. / p.)
      pipeline.py        # per-sample orchestration
      projectio.py       # save/load project JSON
      export.py          # CSV export
    ui/
      main_window.py, sample_list.py, summary_bar.py, alignment_view.py,
      trace_view.py, variant_table.py, qc_panel.py, workers.py, theme.py
  scripts/
    make_demo_data.py    # synthetic AB1/SEQ/reference generator (ABIF writer)
  demo/                  # generated demo dataset
  tests/                 # pytest: parsers, trim, align, consensus, consequence
```

## 11. Implementation constraints

* Never send sequence data over the network.
* Perform alignment and analysis off the UI thread.
* Normalize everything to **reference coordinates**.
* Preserve the raw chromatogram separately from derived calls.
* Make every variant visually auditable against the trace.
* Do not silently overwrite imported `.seq` calls with re-derived AB1 calls; retain both and
  expose disagreements.
* Optimize primarily for **fast human review**, not maximal automated variant calling.

For anything intended to inform actual clinical decisions, the variant-calling/QC rules would
need formal validation beyond simply implementing the software correctly.
