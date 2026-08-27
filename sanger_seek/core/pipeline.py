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
from .phd import load_phd
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
                    phd_path=str(rf.phd) if rf.phd else None,
                    orientation_hint=rf.hint,
                )
            )
        samples.append(sample)
    return samples


def prepare_read(read: Read, cfg: Config) -> None:
    """Parse files and compute file-level derived data (idempotent)."""
    if read.trace is None and read.ab1_path:
        read.trace = load_ab1(read.ab1_path)
    if read.seq_imported is None and read.seq_path and not any(
        item.startswith("Could not read SEQ") for item in read.companion_errors
    ):
        try:
            read.seq_imported = load_seq(read.seq_path)
        except Exception as exc:
            read.companion_errors.append(f"Could not read SEQ companion: {exc}")
    if read.phd_calls is None and read.phd_path and not any(
        item.startswith("Could not read PHD") for item in read.companion_errors
    ):
        try:
            phd = load_phd(read.phd_path)
            read.phd_calls = phd.calls
            read.phd_quals = phd.quals
            read.phd_ploc = phd.ploc
        except Exception as exc:
            read.companion_errors.append(f"Could not read PHD companion: {exc}")

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

    if read.trace is not None and read.phd_calls is not None:
        read.phd_discrepancies = _compare_calls(
            read.trace.calls, read.phd_calls, "PHD"
        )

    read.trim = mott_trim(read.quals, read.n, cfg.trim_cutoff)
    if read.trim == (0, 0) and read.quals is None:
        read.trim = (0, read.n)


def analyze_sample(sample: Sample, reference: Reference | None, cfg: Config) -> None:
    sample.error = ""
    sample.qc_flags = []
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
                read.orientation = read.orientation_override or read.orientation_hint or "?"
                read.error = "read did not align to the reference"
        else:
            read.orientation = read.orientation_override or read.orientation_hint or "?"
        read.qc_flags = _read_qc_flags(read, cfg)

    if reference is not None:
        sample.variants = call_sample_variants(sample, reference.seq, cfg)
        annotate_variants(sample.variants, reference)
    else:
        sample.variants = []
    orientations = {read.orientation for read in sample.reads if read.trace is not None}
    if "F" not in orientations:
        sample.qc_flags.append("Forward AB1 trace is missing or could not be oriented")
    if "R" not in orientations:
        sample.qc_flags.append("Reverse AB1 trace is missing or could not be oriented")
    flagged = sum(bool(read.qc_flags or read.error) for read in sample.reads)
    if flagged:
        sample.qc_flags.append(f"{flagged} read(s) have quality or artifact flags")
    sample.analyzed = True


def _compare_calls(ab1_calls: str, imported: str, label: str) -> DiscrepancyReport:
    if ab1_calls == imported:
        return DiscrepancyReport(0, [])
    if len(ab1_calls) == len(imported):
        positions = [i for i, (a, b) in enumerate(zip(ab1_calls, imported)) if a != b]
        return DiscrepancyReport(len(positions), positions)
    return DiscrepancyReport(
        global_distance(ab1_calls, imported), [],
        note=f"lengths differ (ab1 {len(ab1_calls)}, {label.lower()} {len(imported)})",
    )


def _read_qc_flags(read: Read, cfg: Config) -> list[str]:
    """Return conservative, review-oriented trace artifact flags."""
    flags: list[str] = list(read.companion_errors)
    if read.trace is None:
        flags.append("No AB1 chromatogram")
        return flags
    if read.quals is not None and len(read.quals):
        start, end = read.trim
        retained = max(end - start, 0) / len(read.quals)
        mean_q = float(read.quals[start:end].mean()) if end > start else 0.0
        if mean_q < cfg.low_qual:
            flags.append(f"Low mean base quality (Q{mean_q:.0f})")
        if retained < 0.50:
            flags.append(f"Heavy quality trimming ({retained:.0%} retained)")
    if read.peaks:
        noisy = sum(p.ratio >= 0.50 for p in read.peaks) / len(read.peaks)
        if noisy >= 0.15:
            flags.append(f"Frequent secondary peaks ({noisy:.0%} of calls)")
    if read.alignment is not None and read.alignment.identity < 0.90:
        flags.append(f"Low reference alignment identity ({read.alignment.identity:.1%})")
    if read.discrepancies is not None and read.discrepancies.count:
        flags.append(f"SEQ disagrees with AB1 at {read.discrepancies.count} call(s)")
    if read.phd_discrepancies is not None and read.phd_discrepancies.count:
        flags.append(f"PHD disagrees with AB1 at {read.phd_discrepancies.count} call(s)")
    return flags


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
