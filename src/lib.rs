//! High-performance, UI-independent sequence analysis for DNA Analyzer.

pub mod alignment;
pub mod restriction;

use std::collections::BTreeSet;

/// The molecule type inferred from the input alphabet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleculeKind {
    Dna,
    Rna,
    Mixed,
}

impl MoleculeKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Mixed => "Mixed DNA/RNA",
        }
    }
}

/// The text format recognized around the nucleotide sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Plain,
    Fasta,
    GenBankOrigin,
    Embl,
}

impl InputFormat {
    pub fn label(self) -> &'static str {
        match self {
            Self::Plain => "plain or numbered text",
            Self::Fasta => "FASTA",
            Self::GenBankOrigin => "GenBank ORIGIN",
            Self::Embl => "EMBL SQ",
        }
    }
}

/// Cleaned input and its immediately useful derived values.
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceAnalysis {
    /// Canonical internal representation. RNA uracil is normalized to thymine.
    pub dna: String,
    /// Canonical display representation (uses U for unambiguously RNA input).
    pub display: String,
    pub complement: String,
    pub reverse: String,
    pub reverse_complement: String,
    pub kind: MoleculeKind,
    pub input_format: InputFormat,
    pub invalid_letters: BTreeSet<char>,
    pub gc_percent: f64,
}

impl Default for SequenceAnalysis {
    fn default() -> Self {
        Self {
            dna: String::new(),
            display: String::new(),
            complement: String::new(),
            reverse: String::new(),
            reverse_complement: String::new(),
            kind: MoleculeKind::Dna,
            input_format: InputFormat::Plain,
            invalid_letters: BTreeSet::new(),
            gc_percent: 0.0,
        }
    }
}

/// Parse plain/numbered sequences, FASTA, GenBank ORIGIN, or EMBL SQ text.
///
/// Section markers may have trailing whitespace, and an ORIGIN marker may share
/// its line with the first numbered sequence row. Unsupported alphabetic
/// characters inside sequence rows are reported to the caller.
pub fn analyze_input(raw: &str) -> SequenceAnalysis {
    let mut dna = String::with_capacity(raw.len());
    let mut invalid_letters = BTreeSet::new();
    let mut saw_t = false;
    let mut saw_u = false;
    let has_origin_section = raw.lines().any(|line| origin_remainder(line).is_some());
    let has_embl_section = !has_origin_section && raw.lines().any(is_embl_sequence_header);
    let input_format = if has_origin_section {
        InputFormat::GenBankOrigin
    } else if has_embl_section {
        InputFormat::Embl
    } else if raw.lines().any(|line| line.trim_start().starts_with('>')) {
        InputFormat::Fasta
    } else {
        InputFormat::Plain
    };
    let has_sequence_section = has_origin_section || has_embl_section;
    let mut inside_sequence_section = !has_sequence_section;

    for line in raw.lines() {
        let trimmed = line.trim_start();
        let sequence_text = if let Some(remainder) = origin_remainder(line) {
            inside_sequence_section = true;
            remainder
        } else if is_embl_sequence_header(line) {
            inside_sequence_section = true;
            continue;
        } else if !inside_sequence_section || trimmed.starts_with('>') || trimmed.starts_with(';') {
            continue;
        } else {
            line
        };

        let (sequence_text, reached_terminator) = sequence_text
            .split_once("//")
            .map_or((sequence_text, false), |(sequence, _)| (sequence, true));

        for ch in sequence_text.chars() {
            if !ch.is_alphabetic() {
                continue;
            }
            let upper = ch.to_ascii_uppercase();
            if is_iupac_nucleotide(upper) {
                saw_t |= upper == 'T';
                saw_u |= upper == 'U';
                dna.push(if upper == 'U' { 'T' } else { upper });
            } else {
                invalid_letters.insert(upper);
            }
        }

        if reached_terminator && has_sequence_section {
            break;
        }
    }

    let kind = match (saw_t, saw_u) {
        (false, true) => MoleculeKind::Rna,
        (true, true) => MoleculeKind::Mixed,
        _ => MoleculeKind::Dna,
    };
    let display = display_alphabet(&dna, kind);
    let complement_dna: String = dna.chars().map(complement_base).collect();
    let reverse_dna: String = dna.chars().rev().collect();
    let reverse_complement_dna: String = complement_dna.chars().rev().collect();
    let gc_count = dna
        .bytes()
        .filter(|base| matches!(base, b'G' | b'C'))
        .count();
    let gc_percent = if dna.is_empty() {
        0.0
    } else {
        gc_count as f64 * 100.0 / dna.len() as f64
    };

    SequenceAnalysis {
        dna,
        display,
        complement: display_alphabet(&complement_dna, kind),
        reverse: display_alphabet(&reverse_dna, kind),
        reverse_complement: display_alphabet(&reverse_complement_dna, kind),
        kind,
        input_format,
        invalid_letters,
        gc_percent,
    }
}

fn origin_remainder(line: &str) -> Option<&str> {
    marker_remainder(line, "ORIGIN")
}

fn is_embl_sequence_header(line: &str) -> bool {
    marker_remainder(line, "SQ").is_some()
}

fn marker_remainder<'a>(line: &'a str, marker: &str) -> Option<&'a str> {
    let trimmed = line.trim_start();
    let keyword = trimmed.get(..marker.len())?;
    if !keyword.eq_ignore_ascii_case(marker) {
        return None;
    }

    let remainder = trimmed.get(marker.len()..)?;
    match remainder.chars().next() {
        None => Some(remainder),
        Some(character) if character.is_whitespace() => Some(remainder.trim_start()),
        Some(':') => Some(remainder.get(1..)?.trim_start()),
        Some(_) => None,
    }
}

fn display_alphabet(dna: &str, kind: MoleculeKind) -> String {
    if kind == MoleculeKind::Rna {
        dna.replace('T', "U")
    } else {
        dna.to_owned()
    }
}

pub fn is_iupac_nucleotide(base: char) -> bool {
    matches!(
        base,
        'A' | 'C'
            | 'G'
            | 'T'
            | 'U'
            | 'R'
            | 'Y'
            | 'S'
            | 'W'
            | 'K'
            | 'M'
            | 'B'
            | 'D'
            | 'H'
            | 'V'
            | 'N'
    )
}

pub fn complement_base(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' | 'U' => 'A',
        'R' => 'Y',
        'Y' => 'R',
        'S' => 'S',
        'W' => 'W',
        'K' => 'M',
        'M' => 'K',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N',
        other => other,
    }
}

pub fn reverse_complement(dna: &str) -> String {
    dna.chars().rev().map(complement_base).collect()
}

/// Translate all forward and reverse-complement reading frames.
pub fn six_frame_translation(dna: &str) -> [String; 6] {
    let reverse = reverse_complement(dna);
    std::array::from_fn(|index| {
        if index < 3 {
            translate_frame(dna, index)
        } else {
            translate_frame(&reverse, index - 3)
        }
    })
}

pub fn translate_frame(dna: &str, offset: usize) -> String {
    if offset >= dna.len() {
        return String::new();
    }
    let bytes = dna.as_bytes();
    let codon_count = (bytes.len() - offset) / 3;
    let mut protein = String::with_capacity(codon_count);
    for index in (offset..bytes.len().saturating_sub(2)).step_by(3) {
        protein.push(translate_codon([
            bytes[index],
            bytes[index + 1],
            bytes[index + 2],
        ]));
    }
    protein
}

fn translate_codon(codon: [u8; 3]) -> char {
    match &codon {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

/// Return overlapping, zero-based match offsets.
pub fn find_overlapping(haystack: &str, needle: &str) -> Vec<usize> {
    let query = analyze_input(needle).dna;
    if query.is_empty() || query.len() > haystack.len() {
        return Vec::new();
    }
    haystack
        .as_bytes()
        .windows(query.len())
        .enumerate()
        .filter_map(|(index, window)| (window == query.as_bytes()).then_some(index))
        .collect()
}

/// FASTA header followed by GenBank-style ORIGIN rows.
///
/// This intentionally preserves the desktop app's established export format:
/// one-based row offsets, groups of ten bases, sixty bases per row, and a `//`
/// terminator.
pub fn format_origin_export(name: &str, sequence: &str) -> String {
    let safe_name = name.trim().trim_start_matches('>');
    let safe_name = if safe_name.is_empty() {
        "sequence_1"
    } else {
        safe_name
    };
    let mut out = format!(">{safe_name}\n");
    for (line_index, line) in sequence.as_bytes().chunks(60).enumerate() {
        out.push_str(&format!("{:>9} ", line_index * 60 + 1));
        for (group_index, group) in line.chunks(10).enumerate() {
            if group_index > 0 {
                out.push(' ');
            }
            out.push_str(std::str::from_utf8(group).expect("IUPAC sequence is ASCII"));
        }
        out.push('\n');
    }
    out.push_str("//\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cleans_fasta_and_origin_without_polluting_sequence() {
        let parsed = analyze_input(">sample human gene\n  1 acgu ry n 60\n//");
        assert_eq!(parsed.dna, "ACGTRYN");
        assert_eq!(parsed.kind, MoleculeKind::Rna);
        assert_eq!(parsed.input_format, InputFormat::Fasta);
        assert!(parsed.invalid_letters.is_empty());
    }

    #[test]
    fn reads_only_the_genbank_origin_section() {
        let parsed = analyze_input(
            "LOCUS       SCU49845 12 bp DNA\nDEFINITION  fake header\nORIGIN\n        1 acgt acgtac gt\n//\n",
        );
        assert_eq!(parsed.dna, "ACGTACGTACGT");
        assert_eq!(parsed.input_format, InputFormat::GenBankOrigin);
    }

    #[test]
    fn reads_origin_marker_with_trailing_spaces() {
        let parsed = analyze_input("ORIGIN               \n        1 acgt acgt\n//\n");
        assert_eq!(parsed.dna, "ACGTACGT");
        assert_eq!(parsed.input_format, InputFormat::GenBankOrigin);
        assert!(parsed.invalid_letters.is_empty());
    }

    #[test]
    fn reads_sequence_on_the_origin_marker_line() {
        let parsed = analyze_input(
            "ORIGIN               1 agcgacacat cacacgggct\n\
                    21 caacagtgta gttggtgttc //\n\
             FEATURES             ignored after terminator\n",
        );
        assert_eq!(parsed.dna, "AGCGACACATCACACGGGCTCAACAGTGTAGTTGGTGTTC");
        assert_eq!(parsed.input_format, InputFormat::GenBankOrigin);
        assert!(parsed.invalid_letters.is_empty());
    }

    #[test]
    fn reads_embl_sq_sequence_rows() {
        let parsed = analyze_input(
            "ID   TEST; SV 1; linear; genomic DNA; STD; UNC; 12 BP.\n\
             SQ   Sequence 12 BP; 3 A; 3 C; 3 G; 3 T;\n\
                  acgt acgtac gt 12\n\
             //\n",
        );
        assert_eq!(parsed.dna, "ACGTACGTACGT");
        assert_eq!(parsed.input_format, InputFormat::Embl);
        assert!(parsed.invalid_letters.is_empty());
    }

    #[test]
    fn computes_iupac_reverse_complement() {
        assert_eq!(reverse_complement("ACGTRYN"), "NRYACGT");
    }

    #[test]
    fn translates_forward_and_reverse_frames() {
        assert_eq!(translate_frame("ATGGCTTAA", 0), "MA*");
        let frames = six_frame_translation("ATGGCTTAA");
        assert_eq!(frames[0], "MA*");
        assert_eq!(frames.len(), 6);
    }

    #[test]
    fn reports_overlapping_matches() {
        assert_eq!(find_overlapping("AAAA", "AAA"), vec![0, 1]);
    }

    #[test]
    fn emits_numbered_origin_rows() {
        assert_eq!(
            format_origin_export("demo", "ACGTAC"),
            ">demo\n        1 ACGTAC\n//\n"
        );
    }
}
