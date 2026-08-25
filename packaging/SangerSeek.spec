# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller recipe for the native macOS application bundle."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files

ROOT = Path(SPEC).resolve().parents[1]
ICON = ROOT / "build" / "macos" / "SangerSeek.icns"

datas = collect_data_files("pyqtgraph")
datas += [(str(ROOT / "sanger_seek" / "resources" / "app-icon.png"), "sanger_seek/resources")]

hiddenimports = [
    # SeqIO selects these parsers dynamically from the requested format.
    "Bio.SeqIO.AbiIO",
    "Bio.SeqIO.FastaIO",
    "Bio.SeqIO.InsdcIO",
    # The demo generator imports these inside its function body.
    "Bio.Seq",
    "Bio.SeqFeature",
    "Bio.SeqRecord",
]

a = Analysis(
    [str(ROOT / "packaging" / "launcher.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter", "pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Sanger Seek",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="Sanger Seek",
)
app = BUNDLE(
    coll,
    name="Sanger Seek.app",
    icon=str(ICON),
    bundle_identifier="org.sanger-seek.app",
    info_plist={
        "CFBundleDisplayName": "Sanger Seek",
        "CFBundleName": "Sanger Seek",
        "CFBundleShortVersionString": "0.1.0",
        "CFBundleVersion": "1",
        "LSMinimumSystemVersion": "12.0",
        "NSHighResolutionCapable": True,
        "NSHumanReadableCopyright": "Copyright © 2026 Sanger Seek contributors",
        "CFBundleDocumentTypes": [
            {
                "CFBundleTypeName": "Sanger Seek Project",
                "CFBundleTypeRole": "Editor",
                "LSHandlerRank": "Owner",
                "CFBundleTypeExtensions": ["json"],
            },
            {
                "CFBundleTypeName": "Sanger Chromatogram",
                "CFBundleTypeRole": "Viewer",
                "LSHandlerRank": "Alternate",
                "CFBundleTypeExtensions": ["ab1", "abi", "ab", "seq"],
            },
            {
                "CFBundleTypeName": "Sequence Reference",
                "CFBundleTypeRole": "Viewer",
                "LSHandlerRank": "Alternate",
                "CFBundleTypeExtensions": ["fa", "fasta", "fna", "gb", "gbk", "genbank"],
            },
        ],
    },
)
