# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['dna_translator.py'],
    pathex=[],
    binaries=[],
    datas=[('muscle-osx-arm64', '.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='DNAanalyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['app.icns'],
)
app = BUNDLE(
    exe,
    name='DNAanalyzer.app',
    icon='app.icns',
    bundle_identifier=None,
)
