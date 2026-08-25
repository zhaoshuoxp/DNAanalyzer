# Third-party notices

DNA Analyzer itself is released under the MIT License. Packaged distributions may include the following separate components.

## MUSCLE 5

MUSCLE is Copyright © Robert C. Edgar and is distributed under the GNU General Public License, version 3.

- Project and source: <https://github.com/rcedgar/muscle>
- Releases: <https://github.com/rcedgar/muscle/releases>
- License: <https://github.com/rcedgar/muscle/blob/main/LICENSE>

MUSCLE is executed as a separate process. It is not linked into the DNA Analyzer Rust executable.

## GCC runtime libraries (macOS package)

The macOS application bundle includes `libgomp`, `libstdc++`, and, when required by the architecture, `libgcc_s` from GCC 11 so that the bundled MUSCLE executable does not require Homebrew on the destination computer.

GCC runtime libraries are covered by the GNU General Public License, version 3, with the GCC Runtime Library Exception, version 3.1.

- GCC source: <https://gcc.gnu.org/git/gcc.git>
- GCC licenses: <https://gcc.gnu.org/onlinedocs/libstdc++/manual/license.html>
- Runtime Library Exception: <https://www.gnu.org/licenses/gcc-exception-3.1.html>
