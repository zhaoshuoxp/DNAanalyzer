# Architecture

DNA Analyzer 2.0 separates the native UI from a testable analysis library.

## Layers

1. `src/lib.rs` owns parsing, complements, translation, search, and ORIGIN export. It has no GUI or operating-system dependencies.
2. `src/restriction.rs` loads the compile-time enzyme table once, represents IUPAC symbols as four-bit masks, and distributes enzyme scans across the Rayon thread pool.
3. `src/alignment.rs` tries every discovered MUSCLE executable, isolates each run in a temporary directory, checks process status, and validates the returned FASTA alignment. If no candidate can run, a native center-star progressive aligner performs Needleman–Wunsch pairwise alignment and merges the center gap profiles.
4. `src/app.rs` renders the original single-page egui layout. Restriction scanning and MUSCLE run on workers connected through channels, so repaint and editing remain responsive. Search and Multiple Alignment launch the same signed executable in dedicated window modes; this gives each tool a true native window and avoids unreliable eframe child-viewport teardown on Intel macOS.

## Portability

egui/eframe supplies a native window and GPU-backed renderer on Windows, macOS, and Linux. The Windows target disables the console subsystem. Release packages include the matching MUSCLE executable. On macOS, native window persistence is disabled to avoid the known AppKit Touch Bar KVO teardown crash on affected Intel Macs. Packaging copies the required GCC 11 runtime dylibs into `Contents/Frameworks`, rewrites MUSCLE and dylib load commands to bundle-relative paths, then ad-hoc signs the complete app. The destination Mac therefore does not need Homebrew. The Rust fallback keeps MSA available if any external engine still fails.

## Data conventions

- Internal nucleic-acid strings use uppercase DNA (`U` is normalized to `T`).
- Unambiguously RNA input is rendered with `U` in derived sequences.
- Restriction positions are one-based top-strand cleavage coordinates.
- Search offsets are stored zero-based internally and shown one-based.
- Export deliberately uses the project's numbered ORIGIN-style layout.
