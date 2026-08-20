Yes — I’d explicitly structure the app around the `.seq` + `.ab1` pairing and make the primary review interface feel like a modernized version of the screenshot you attached.

## Revised client-side Sanger viewer spec

### 1. Browser-only architecture

* **HTML + CSS + JavaScript/TypeScript only**
* No server, API, database, or sequence upload.
* Parse and analyze everything using:

  * `FileReader`
  * `ArrayBuffer` / `DataView`
  * Web Workers for alignment and analysis
  * Canvas/WebGL or SVG for chromatogram rendering
* Optional IndexedDB for locally saved projects.
* Deploy as static files and optionally support offline/PWA use.

### 2. File handling

**`.ab1` — raw chromatogram source**

Parse the ABIF binary directly in JavaScript to extract:

* A/C/G/T trace arrays (`DATA` records)
* called bases (`PBAS`)
* peak/sample positions (`PLOC`)
* base quality values (`PCON`)
* sequencing metadata where useful

The `.ab1` should be the authoritative source for the **actual electropherogram**, peak quality, mixed peaks, noise, and manual visual review.

**`.seq` — called sequence**

Parse as plain-text nucleotide sequence:

* normalize whitespace/header formatting
* associate it with the corresponding `.ab1`
* use it as an imported base-call sequence when present
* compare it against the calls extracted from the AB1 and flag discrepancies if useful

If both exist:

```text
sample_F.ab1   -> trace + peaks + qualities + embedded base calls
sample_F.seq   -> exported/called nucleotide sequence
```

Treat them as two representations of the **same read**, not independent sequencing evidence.

The app should also work with `.ab1` alone because AB1 generally contains the base calls itself.

**Reference**

* FASTA for sequence-only reference.
* GenBank optionally for CDS/exon/gene annotation and amino-acid consequences.

### 3. Analysis pipeline

For each sample/read:

```text
Load files
   ↓
Parse AB1 / SEQ
   ↓
QC + low-quality trimming
   ↓
Determine F/R orientation
   ↓
Reverse-complement when required
   ↓
Align read to reference
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
* ambiguous/mixed bases
* forward/reverse discordance

With an annotated coding reference, additionally calculate:

* synonymous
* missense
* nonsense
* frameshift
* in-frame indel
* codon change
* protein change, e.g. `c.944C>T / p.Thr315Ile`

### 4. Main review interface

Use the attached interface as the functional inspiration, but simplify it substantially.

A good screen hierarchy would be:

```text
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

### 5. Chromatogram UX

Compared with the screenshot, eliminate most of the dense diagnostic-looking clutter.

The trace viewer should have:

* smooth A/C/G/T traces
* synchronized base calls above peaks
* reference coordinate ruler
* quality shading
* selected-position crosshair
* variant flags
* zoom/pan
* drag navigation
* keyboard left/right between bases
* next/previous variant shortcuts

When a user clicks `c.944C>T`, every view should instantly center on that nucleotide.

Forward and reverse traces should be displayed in the **same reference orientation**, so users do not have to mentally reverse one strand.

### 6. Noise / mixed-peak visualization

The screenshot contains useful signal/noise information, but it can be communicated more cleanly.

For every called base calculate/display:

* primary peak height
* secondary peak height
* secondary/primary ratio
* local baseline/noise
* AB1 quality score
* distance from sequencing read ends

For candidate heterogeneity/mixed peaks, show a small indicator such as:

```text
C 72%
T 28%
```

and let the user inspect the original trace rather than hiding this behind an automated call.

### 7. Variant review table

Keep a polished fixed panel rather than the dense table in the screenshot:

| Variant    | Protein     | Fwd | Rev | Trace | Confidence |
| ---------- | ----------- | --: | --: | ----- | ---------- |
| c.944C>T   | p.Thr315Ile |   ✓ |   ✓ | Clean | High       |
| c.1012delA | p.X338fs    |   ✓ |   ✓ | Clean | High       |
| c.1173G>A  | synonymous  |   ✓ |   ? | Noisy | Review     |

Clicking any row jumps to that trace.

Filters:

* all differences
* coding only
* missense/nonsense
* indels
* mixed peaks
* strand disagreement
* low confidence

### 8. Forward/reverse reconciliation

Make this a first-class object rather than two unrelated traces.

For each reference position:

```text
Reference: C
Forward:   T   Q=34
Reverse:   T   Q=37

Consensus: T
Status:    Supported on both strands
```

Versus:

```text
Reference: C
Forward:   T   Q=31
Reverse:   C   Q=36

Consensus: ?
Status:    Strand disagreement — review
```

This should drive the mutation flags.

### 9. Multi-read / multi-amplicon support

Data model:

```text
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
     ├─ Amplicons[]
     ├─ Consensus
     └─ Variants[]
```

If multiple contigs/amplicons exist, map each independently onto the same reference and produce one reference-coordinate-based variant list.

### 10. Implementation decisions for an LLM

I’d give the implementing model these constraints:

* Use **TypeScript** rather than loose JavaScript.
* Separate modules:

  * `abif-parser`
  * `seq-parser`
  * `fasta-parser`
  * `genbank-parser`
  * `alignment-worker`
  * `variant-caller`
  * `consequence-annotator`
  * `trace-renderer`
  * `project-store`
* Never send sequence data over the network.
* Perform expensive alignment in Web Workers.
* Normalize everything to **reference coordinates**.
* Preserve the raw chromatogram separately from derived calls.
* Make every variant visually auditable against the trace.
* Do not silently overwrite imported `.seq` calls with re-derived AB1 calls; retain both and expose disagreements.
* Optimize primarily for **fast human review**, not maximal automated variant calling.

The key design change from the attached software would be: **keep its information density available, but use progressive disclosure**. The default screen should immediately answer *“where are the differences, do forward/reverse support them, and what do the peaks actually look like?”* Everything else can sit behind expandable QC/details panels.

For anything intended to inform actual clinical decisions, the variant-calling/QC rules would need formal validation beyond simply implementing the software correctly.

