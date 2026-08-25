#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "$script_dir/.." && pwd)"
architecture="$(uname -m)"

if [[ "$architecture" == "arm64" ]]; then
  muscle_name="muscle-osx-arm64"
else
  muscle_name="muscle-osx-x86"
fi

binary_path="${DNA_ANALYZER_BINARY:-$project_dir/target/release/dna-analyzer}"
app_path="${DNA_ANALYZER_APP_OUTPUT:-$project_dir/target/release/bundle/macos/DNA Analyzer.app}"

if [[ -z "${DNA_ANALYZER_BINARY:-}" ]]; then
  cargo build --manifest-path "$project_dir/Cargo.toml" --locked --release
fi

if ! gcc_prefix="$(brew --prefix gcc@11 2>/dev/null)" || [[ ! -d "$gcc_prefix" ]]; then
  echo "gcc@11 is required only while packaging. Install it with: brew install gcc@11" >&2
  exit 1
fi

runtime_dir="$(dirname "$(find -L "$gcc_prefix" -type f -name 'libgomp.1.dylib' -print -quit)")"
if [[ ! -f "$runtime_dir/libgomp.1.dylib" || ! -f "$runtime_dir/libstdc++.6.dylib" ]]; then
  echo "Could not locate the gcc@11 runtime libraries under $gcc_prefix" >&2
  exit 1
fi

contents="$app_path/Contents"
macos_dir="$contents/MacOS"
resources_dir="$contents/Resources"
frameworks_dir="$contents/Frameworks"
rm -rf "$app_path"
mkdir -p "$macos_dir" "$resources_dir" "$frameworks_dir"

cp "$binary_path" "$macos_dir/dna-analyzer"
cp "$project_dir/$muscle_name" "$resources_dir/$muscle_name"
cp "$project_dir/app.icns" "$resources_dir/app.icns"
cp "$project_dir/packaging/macos/Info.plist" "$contents/Info.plist"
cp "$project_dir/THIRD_PARTY_NOTICES.md" "$resources_dir/THIRD_PARTY_NOTICES.md"

for library in libgomp.1.dylib libstdc++.6.dylib libgcc_s.1.1.dylib; do
  if [[ -f "$runtime_dir/$library" ]]; then
    cp -L "$runtime_dir/$library" "$frameworks_dir/$library"
  fi
done

while IFS= read -r dependency; do
  library="$(basename "$dependency")"
  if [[ -f "$frameworks_dir/$library" ]]; then
    install_name_tool \
      -change "$dependency" "@executable_path/../Frameworks/$library" \
      "$resources_dir/$muscle_name"
  fi
done < <(otool -L "$resources_dir/$muscle_name" | tail -n +2 | awk '{print $1}')

for bundled_library in "$frameworks_dir"/*.dylib; do
  library_name="$(basename "$bundled_library")"
  install_name_tool -id "@rpath/$library_name" "$bundled_library"
  while IFS= read -r dependency; do
    dependency_name="$(basename "$dependency")"
    if [[ -f "$frameworks_dir/$dependency_name" ]]; then
      install_name_tool \
        -change "$dependency" "@loader_path/$dependency_name" \
        "$bundled_library"
    fi
  done < <(otool -L "$bundled_library" | tail -n +2 | awk '{print $1}')
  codesign --force --sign - "$bundled_library"
done

chmod +x "$macos_dir/dna-analyzer" "$resources_dir/$muscle_name"
codesign --force --sign - "$resources_dir/$muscle_name"
codesign --force --sign - "$macos_dir/dna-analyzer"
codesign --deep --force --sign - "$app_path"

echo "Created $app_path"
