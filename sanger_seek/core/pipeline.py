"""Per-sample analysis orchestration.

Load files -> parse AB1/SEQ -> QC + trimming -> orient + align -> peak
metrics -> per-read candidates -> strand reconciliation -> annotation.
"""

from __future__ import annotations

from pathlib import Path

from .abif import load_ab1
from .align import global_distance, orient_and_align
from .consensus import call_sample_variants
from .consequence import annotate_variants
from .model import Config, DiscrepancyReport, Project, Read, Reference, Sample
from .pairing import ScanResult
from .peaks import peak_metrics
from .seqfile import load_seq
from .trim import mott_trim


def build_samples_from_scan(scan: ScanResult) -> list[Sample]:
    samples: list[Sample] = []
    for skey in sorted(scan.samples):
        reads_map = scan.samples[skey]
        sample = Sample(key=skey, name=skey)
        for rkey in sorted(reads_map):
            rf = reads_map[rkey]
            sample.reads.append(
                Read(
                    id=f"{skey}/{rkey}",
                    label=rf.stem,
                    ab1_path=str(rf.ab1) if rf.ab1 else None,
                    seq_path=str(rf.seq) if rf.seq else None,
                    orientation_hint=rf.hint,
                )
            )
        samples.append(sample)
    return samples


def prepare_read(read: Read, cfg: Config) -> None:
    """Parse files and compute file-level derived data (idempotent)."""
    if read.trace is None and read.ab1_path:
        read.trace = load_ab1(read.ab1_path)
    if read.seq_imported is None and read.seq_path:
        read.seq_imported = load_seq(read.seq_path)

    if read.trace is not None:
        read.calls = read.trace.calls
        read.quals = read.trace.quals
        if read.peaks is None:
            read.peaks = peak_metrics(read.trace)
    elif read.seq_imported is not None:
        read.calls = read.seq_imported
        read.quals = None
    else:
        raise ValueError(f"{read.label}: no .ab1 or .seq data")

    # Imported .seq vs AB1-embedded calls: keep both, surface disagreements.
    if read.trace is not None and read.seq_imported is not None:
        pbas, imported = read.trace.calls, read.seq_imported
        if pbas == imported:
            read.discrepancies = DiscrepancyReport(0, [])
        elif len(pbas) == len(imported):
            pos = [i for i, (a, b) in enumerate(zip(pbas, imported)) if a != b]
            read.discrepancies = DiscrepancyReport(len(pos), pos)
        else:
            d = global_distance(pbas, imported)
            read.discrepancies = DiscrepancyReport(
                d, [], note=f"lengths differ (ab1 {len(pbas)}, seq {len(imported)})"
            )

    read.trim = mott_trim(read.quals, read.n, cfg.trim_cutoff)
    if read.trim == (0, 0) and read.quals is None:
        read.trim = (0, read.n)


def analyze_sample(sample: Sample, reference: Reference | None, cfg: Config) -> None:
    sample.error = ""
    for read in sample.reads:
        read.error = ""
        try:
            prepare_read(read, cfg)
        except Exception as e:  # keep other reads usable
            read.error = str(e)
            continue
        if reference is not None:
            read.alignment = orient_and_align(read, reference.seq, cfg)
            if read.alignment is not None:
                read.orientation = read.alignment.orientation
            else:
                read.orientation = read.orientation_hint or "?"
                read.error = "read did not align to the reference"
        else:
            read.orientation = read.orientation_hint or "?"

    if reference is not None:
        sample.variants = call_sample_variants(sample, reference.seq, cfg)
        annotate_variants(sample.variants, reference)
    else:
        sample.variants = []
    sample.analyzed = True


def load_project_inputs(paths: list[str | Path], project: Project) -> ScanResult:
    """Scan paths and merge resulting samples into the project (no analysis)."""
    from .pairing import scan_paths

    scan = scan_paths(list(paths))
    for new_sample in build_samples_from_scan(scan):
        existing = project.sample_by_key(new_sample.key)
        if existing is None:
            project.samples.append(new_sample)
        else:
            known = {r.id for r in existing.reads}
            for r in new_sample.reads:
                if r.id not in known:
                    existing.reads.append(r)
                    existing.analyzed = False
    project.samples.sort(key=lambda s: s.key.lower())
    return scan
