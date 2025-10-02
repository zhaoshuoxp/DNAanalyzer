# DNA analyzer

## Overview

This is a Python-based GUI tool built with **PyQt5** and **Biopython** for sequence analysis. It allows you to:

- Input DNA or RNA sequences.
- Automatically detect RNA and warn about A:U pairing.
- Remove invalid characters and display warnings in red.
- Display **complement**, **reverse**, and **reverse-complement** sequences.
- Perform **six-frame translation** (optional display).
- Analyze **restriction enzyme sites** with keyword search/filter.
- Save sequences in **standard FASTA/ORIGIN format**.
- Search subsequences and highlight matches in a formatted display.

![GUI](https://raw.githubusercontent.com/zhaoshuoxp/DNAanalyzer/refs/heads/main//screenshot1.png)

![GUI](https://raw.githubusercontent.com/zhaoshuoxp/DNAanalyzer/refs/heads/main//screenshot2.png)

------

## Requirements

- Python 3.10+ (installed via Homebrew recommended)
- PyQt5
- Biopython

Install dependencies:

```
pip install pyqt5 biopython
```

------

## Usage

Run the GUI:

```
python dna_translator.py
```

### Features

1. **Input Sequence**: Paste or type DNA/RNA sequence.
2. **Complement / Reverse / Reverse-Complement**: Displayed in separate areas.
3. **Six-frame Translation**: Click the "Translate 6-Frames" button to show.
4. **Restriction Enzyme Analysis**: Click "Show Restriction Enzymes" to display sites. Filter by enzyme name if needed.
5. **Search Subsequence**: Click the button to open a search dialog. Subsequence matches are highlighted in **ORIGIN/FASTA format**.
6. **Save as FASTA**: Save the sequence in a standard, formatted FASTA file.

------

## macOS App Packaging

You can create a standalone macOS `.app` using **pyinstaller**.

### 1. Install PyInstaller

```
pip install pyinstaller
```

### 2. Package the app

```
pyinstaller --name "DNAanalyzer" --onefile --windowed dna_translator.py
```

- `--onefile` generates a single executable.
- `--windowed` prevents a terminal window from opening.

After packaging, the `.app` will be inside the `dist` folder.

### 3. Gatekeeper (macOS security)

If macOS blocks the app:

```
xattr -d com.apple.quarantine /path/to/DNAanalyzer.app
```

- Or right-click → Open → confirm.

------

## Notes

- Input sequences are automatically cleaned (invalid characters removed).
- RNA sequences are automatically detected, and A:U pairing is used.
- Restriction enzyme analysis uses **Biopython AllEnzymes**.

------

## License

MIT License