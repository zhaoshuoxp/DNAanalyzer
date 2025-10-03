# DNA analyzer

## Overview

This is a Python-based GUI tool built with **PyQt5** and **Biopython** for sequence analysis. It allows you to:

- Input DNA or RNA sequences.
- Display **complement**, **reverse**, and **reverse-complement** sequences.
- Perform **six-frame translation** (optional display).
- Analyze **restriction enzyme sites** with keyword search/filter.
- Save sequences in **standard FASTA/ORIGIN format**.
- **Search** subsequences and **highlight** matches in a formatted display.
- **BLAST** search (remote via NCBI)
- Multiple Sequence **Alignment** (MSA) using MUSCLE

![GUI](https://raw.githubusercontent.com/zhaoshuoxp/DNAanalyzer/refs/heads/main//screenshot1.png)

![GUI](https://raw.githubusercontent.com/zhaoshuoxp/DNAanalyzer/refs/heads/main/screenshot2.png)

![GUI](https://raw.githubusercontent.com/zhaoshuoxp/DNAanalyzer/refs/heads/main/screenshot3.png)

------

## Requirements

- PyQt5
- Biopython
- Python 3.8+
- MUSCLE executable (place in the same directory as the script if using MSA)
- Internet access if using NCBI BLAST remote search

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
7. **Run BLAST**: re-direct to NCBI blast webpage.
- **Run Multiple Sequence Alignment (MSA)**:
  - Add/remove sequences dynamically
  - MUSCLE executable auto-detection by OS:
    - `muscle-osx-x86` for macOS
    - `muscle-linux-x86` for Linux
    - `muscle-win64.exe` for Windows
  - Aligned output shown with mismatches in **red**

------

## Package into Single App (macOS example)

You can create a standalone macOS `.app` using **pyinstaller**.

### 1. Install PyInstaller

```
pip install pyinstaller
```

### 2. Package the app

```
pyinstaller --name "DNAanalyzer" --add-data "muscle-osx-x86:." --onefile --windowed dna_translator.py
```

- `--onefile` generates a single executable.
- `--windowed` prevents a terminal window from opening.
- `--add-data` includes the MUSCLE binary for aligning.

After packaging, the `.app` will be inside the `dist` folder.

### 3. Gatekeeper (macOS security)

If macOS blocks the app:

```
xattr -d com.apple.quarantine /path/to/DNAanalyzer.app
```

- Or allow from **System Preferences → Security & Privacy**.

------

## Notes

- Input sequences are automatically cleaned (invalid characters removed).
- RNA sequences are automatically detected, and A:U pairing is used.
- Restriction enzyme analysis uses **Biopython AllEnzymes**.

------

## License

MIT License