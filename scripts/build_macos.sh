#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

APP_NAME="Sanger Seek"
VERSION="0.1.0"
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if command -v python3.13 >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python3.13)"
    else
        PYTHON_BIN="$(command -v python3)"
    fi
fi
VENV_DIR="${BUILD_VENV:-.venv-build}"
ICON_SOURCE="sanger_seek/resources/app-icon.png"
ICONSET="build/macos/SangerSeek.iconset"
ICNS="build/macos/SangerSeek.icns"
DMG="dist/Sanger-Seek-${VERSION}.dmg"

if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "The .app and DMG build must run on macOS." >&2
    exit 2
fi

"$PYTHON_BIN" -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install --upgrade pip
"$VENV_DIR/bin/python" -m pip install -e '.[build]'

mkdir -p build/macos dist
rm -rf "$ICONSET"
mkdir -p "$ICONSET"
for size in 16 32 128 256 512; do
    sips -z "$size" "$size" "$ICON_SOURCE" --out "$ICONSET/icon_${size}x${size}.png" >/dev/null
    double=$((size * 2))
    sips -z "$double" "$double" "$ICON_SOURCE" --out "$ICONSET/icon_${size}x${size}@2x.png" >/dev/null
done
iconutil -c icns "$ICONSET" -o "$ICNS"

"$VENV_DIR/bin/pyinstaller" --noconfirm --clean packaging/SangerSeek.spec

APP_PATH="dist/${APP_NAME}.app"
if [[ -n "${APPLE_SIGN_IDENTITY:-}" ]]; then
    codesign --force --deep --options runtime --timestamp \
        --entitlements packaging/entitlements.plist \
        --sign "$APPLE_SIGN_IDENTITY" "$APP_PATH"
else
    codesign --force --deep --sign - "$APP_PATH"
fi
codesign --verify --deep --strict --verbose=2 "$APP_PATH"

STAGING="$(mktemp -d "${TMPDIR:-/tmp}/sanger-seek-dmg.XXXXXX")"
cleanup() { rm -rf "$STAGING"; }
trap cleanup EXIT
ditto "$APP_PATH" "$STAGING/${APP_NAME}.app"
ln -s /Applications "$STAGING/Applications"
rm -f "$DMG"
hdiutil create -volname "$APP_NAME" -srcfolder "$STAGING" \
    -ov -format UDZO "$DMG"

if [[ -n "${APPLE_NOTARY_PROFILE:-}" ]]; then
    xcrun notarytool submit "$DMG" --keychain-profile "$APPLE_NOTARY_PROFILE" --wait
    xcrun stapler staple "$DMG"
fi

echo "Built $APP_PATH"
echo "Built $DMG"
