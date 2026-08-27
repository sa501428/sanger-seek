"""Project save/load: JSON referencing the original files on disk.

Raw data always stays in the source .ab1/.seq files; the project file only
records which files belong together plus analysis settings. Loading re-parses
and re-analyzes so results always reflect the current files.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

from .model import Config, Project, Read, Sample

FORMAT_VERSION = 2


def save_project(project: Project, path: str | Path) -> None:
    def sample_doc(s: Sample) -> dict:
        return {
            "key": s.key,
            "name": s.name,
            "reads": [
                {
                    "id": r.id,
                    "label": r.label,
                    "ab1": r.ab1_path,
                    "seq": r.seq_path,
                    "phd": r.phd_path,
                    "hint": r.orientation_hint,
                    "orientation_override": r.orientation_override,
                }
                for r in s.reads
            ],
        }

    doc = {
        "format": "sanger-seek-project",
        "version": FORMAT_VERSION,
        "reference": project.reference.path if project.reference else None,
        "config": dataclasses.asdict(project.config),
        "samples": [sample_doc(s) for s in project.samples],
    }
    p = Path(path)
    p.write_text(json.dumps(doc, indent=2))
    project.path = str(p)


def load_project(path: str | Path) -> tuple[Project, list[str]]:
    """Load a project skeleton. Returns (project, warnings). No analysis is run."""
    p = Path(path)
    doc = json.loads(p.read_text())
    if doc.get("format") != "sanger-seek-project":
        raise ValueError(f"{p.name}: not a sanger-seek project file")

    warnings: list[str] = []
    cfg = Config()
    for k, v in (doc.get("config") or {}).items():
        if hasattr(cfg, k):
            setattr(cfg, k, v)

    project = Project(config=cfg, path=str(p))

    ref_path = doc.get("reference")
    if ref_path:
        if Path(ref_path).exists():
            from .reference import load_reference

            project.reference = load_reference(ref_path)
        else:
            warnings.append(f"reference file missing: {ref_path}")

    def load_sample_doc(sdoc: dict) -> Sample | None:
        sample = Sample(key=sdoc["key"], name=sdoc.get("name", sdoc["key"]))
        for rdoc in sdoc.get("reads", []):
            ab1, seq, phd = rdoc.get("ab1"), rdoc.get("seq"), rdoc.get("phd")
            if ab1 and not Path(ab1).exists():
                warnings.append(f"missing file: {ab1}")
                ab1 = None
            if seq and not Path(seq).exists():
                warnings.append(f"missing file: {seq}")
                seq = None
            if phd and not Path(phd).exists():
                warnings.append(f"missing file: {phd}")
                phd = None
            if not ab1 and not seq and not phd:
                continue
            sample.reads.append(
                Read(
                    id=rdoc.get("id") or f"{sample.key}/{rdoc.get('label', '?')}",
                    label=rdoc.get("label", Path(ab1 or seq or phd).stem),
                    ab1_path=ab1,
                    seq_path=seq,
                    phd_path=phd,
                    orientation_hint=rdoc.get("hint"),
                    orientation_override=rdoc.get("orientation_override"),
                )
            )
        return sample if sample.reads else None

    control_doc = doc.get("wt_control")
    if control_doc:
        # Version-1 projects had a privileged WT lane.  In the case-centric
        # workflow it is migrated to an ordinary case so no data is lost.
        control = load_sample_doc(control_doc)
        if control:
            control.key = control.name
            project.samples.append(control)
    for sdoc in doc.get("samples", []):
        sample = load_sample_doc(sdoc)
        if sample:
            project.samples.append(sample)
    return project, warnings
