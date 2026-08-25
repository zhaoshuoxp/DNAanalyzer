//! MUSCLE discovery, execution, and FASTA result parsing.

use std::collections::HashSet;
use std::env;
use std::ffi::OsString;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentResult {
    pub labels: Vec<String>,
    pub sequences: Vec<String>,
    pub engine: String,
}

/// Prefer MUSCLE, but guarantee cross-platform functionality with the native
/// Rust progressive aligner when no runnable external binary is available.
pub fn run_alignment(sequences: &[String]) -> Result<AlignmentResult, String> {
    let cleaned: Vec<String> = sequences
        .iter()
        .map(|sequence| clean_alignment_sequence(sequence))
        .filter(|sequence| !sequence.is_empty())
        .collect();
    if cleaned.len() < 2 {
        return Err("Enter at least two valid DNA/RNA or protein sequences.".to_owned());
    }

    for muscle in muscle_candidates() {
        if let Ok(result) = run_muscle_binary(&cleaned, &muscle) {
            return Ok(result);
        }
    }

    progressive_alignment(&cleaned).map(|mut result| {
        result.engine = "Built-in Rust fallback".to_owned();
        result
    })
}

fn run_muscle_binary(
    cleaned: &[String],
    muscle: &std::path::Path,
) -> Result<AlignmentResult, String> {
    let directory = tempfile::tempdir()
        .map_err(|error| format!("Could not create a temporary folder: {error}"))?;
    let input_path = directory.path().join("input.fasta");
    let output_path = directory.path().join("output.fasta");
    let mut fasta = String::new();
    for (index, sequence) in cleaned.iter().enumerate() {
        fasta.push_str(&format!(">Seq{}\n{}\n", index + 1, sequence));
    }
    fs::write(&input_path, fasta)
        .map_err(|error| format!("Could not write MUSCLE input: {error}"))?;

    let output = Command::new(muscle)
        .arg("-align")
        .arg(&input_path)
        .arg("-output")
        .arg(&output_path)
        .output()
        .map_err(|error| format!("Could not start {}: {error}", muscle.display()))?;
    if !output.status.success() {
        let details = String::from_utf8_lossy(&output.stderr);
        return Err(format!(
            "MUSCLE failed ({}): {}",
            output.status,
            details.trim()
        ));
    }

    let aligned = fs::read_to_string(&output_path)
        .map_err(|error| format!("Could not read MUSCLE output: {error}"))?;
    let mut result = parse_fasta_alignment(&aligned)?;
    result.engine = format!("MUSCLE ({})", muscle.display());
    Ok(result)
}

/// Find a bundled binary first, then a user override, then the system PATH.
pub fn detect_muscle() -> Option<PathBuf> {
    muscle_candidates().into_iter().next()
}

fn muscle_candidates() -> Vec<PathBuf> {
    let binary_name = platform_binary_name();
    let mut candidates = Vec::new();

    if let Some(override_path) = env::var_os("MUSCLE_PATH") {
        candidates.push(PathBuf::from(override_path));
    }
    if let Ok(executable) = env::current_exe()
        && let Some(directory) = executable.parent()
    {
        candidates.push(directory.join(binary_name));
        candidates.push(directory.join("muscle"));
        candidates.push(directory.join("resources").join(binary_name));
        candidates.push(directory.join("..").join("Resources").join(binary_name));
    }
    if let Ok(directory) = env::current_dir() {
        candidates.push(directory.join(binary_name));
        candidates.push(directory.join("muscle"));
    }
    candidates.extend(path_candidates(["muscle", "muscle5", "muscle.exe"]));

    let mut seen = HashSet::new();
    candidates
        .into_iter()
        .filter(|candidate| candidate.is_file() && seen.insert(candidate.clone()))
        .collect()
}

fn platform_binary_name() -> &'static str {
    if cfg!(target_os = "windows") {
        "muscle-win64.exe"
    } else if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
        "muscle-osx-arm64"
    } else if cfg!(target_os = "macos") {
        "muscle-osx-x86"
    } else {
        "muscle-linux-x86"
    }
}

fn path_candidates<const N: usize>(names: [&str; N]) -> Vec<PathBuf> {
    let Some(path) = env::var_os("PATH") else {
        return Vec::new();
    };
    env::split_paths(&path)
        .flat_map(|directory| names.iter().map(move |name| directory.join(name)))
        .collect()
}

pub fn clean_alignment_sequence(raw: &str) -> String {
    let letters: String = raw
        .lines()
        .filter(|line| !line.trim_start().starts_with('>'))
        .flat_map(|line| line.chars())
        .filter(|character| {
            character.is_ascii_alphabetic() || *character == '-' || *character == '*'
        })
        .map(|character| character.to_ascii_uppercase())
        .collect();

    let nucleotide = letters
        .chars()
        .all(|character| "ACGTURYSWKMBDHVN-".contains(character));
    let alphabet = if nucleotide {
        "ACGTURYSWKMBDHVN-"
    } else {
        "ABCDEFGHIKLMNPQRSTVWXYZ*-"
    };
    letters
        .chars()
        .filter(|character| alphabet.contains(*character))
        .map(|character| {
            if character == 'U' && nucleotide {
                'T'
            } else {
                character
            }
        })
        .collect()
}

pub fn parse_fasta_alignment(fasta: &str) -> Result<AlignmentResult, String> {
    let mut labels = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    for line in fasta.lines() {
        if let Some(label) = line.strip_prefix('>') {
            labels.push(label.trim().to_owned());
            sequences.push(String::new());
        } else if let Some(sequence) = sequences.last_mut() {
            sequence.push_str(line.trim());
        }
    }
    if labels.len() < 2 || labels.len() != sequences.len() {
        return Err("MUSCLE returned an invalid alignment.".to_owned());
    }
    let width = sequences.first().map(String::len).unwrap_or(0);
    if width == 0 || sequences.iter().any(|sequence| sequence.len() != width) {
        return Err("MUSCLE returned rows of different lengths.".to_owned());
    }
    Ok(AlignmentResult {
        labels,
        sequences,
        engine: "MUSCLE".to_owned(),
    })
}

/// Center-star progressive MSA using Needleman–Wunsch pairwise alignments.
/// This is deliberately dependency-free and acts as a reliable fallback.
fn progressive_alignment(sequences: &[String]) -> Result<AlignmentResult, String> {
    let center_index = sequences
        .iter()
        .enumerate()
        .max_by_key(|(_, sequence)| sequence.len())
        .map(|(index, _)| index)
        .ok_or_else(|| "No sequences to align.".to_owned())?;
    let center = &sequences[center_index];
    let mut master_center = center.clone();
    let mut rows = vec![(center_index, center.clone())];

    for (index, sequence) in sequences.iter().enumerate() {
        if index == center_index {
            continue;
        }
        let (new_center, new_sequence) = needleman_wunsch(center, sequence)?;
        let existing: Vec<String> = rows.iter().map(|(_, row)| row.clone()).collect();
        let (merged_center, merged_rows, merged_new) =
            merge_center_alignments(&master_center, &existing, &new_center, &new_sequence);
        master_center = merged_center;
        for ((_, row), merged) in rows.iter_mut().zip(merged_rows) {
            *row = merged;
        }
        rows.push((index, merged_new));
    }

    rows.sort_unstable_by_key(|(index, _)| *index);
    Ok(AlignmentResult {
        labels: (1..=rows.len())
            .map(|index| format!("Seq{index}"))
            .collect(),
        sequences: rows.into_iter().map(|(_, sequence)| sequence).collect(),
        engine: "Built-in Rust fallback".to_owned(),
    })
}

fn needleman_wunsch(left: &str, right: &str) -> Result<(String, String), String> {
    let left = left.as_bytes();
    let right = right.as_bytes();
    let columns = right.len() + 1;
    let cells = (left.len() + 1)
        .checked_mul(columns)
        .ok_or_else(|| "Sequences are too large for built-in alignment.".to_owned())?;
    if cells > 25_000_000 {
        return Err(
            "Sequences are too large for the built-in aligner; install MUSCLE or set MUSCLE_PATH."
                .to_owned(),
        );
    }

    const DIAGONAL: u8 = 0;
    const UP: u8 = 1;
    const LEFT: u8 = 2;
    const GAP: i32 = -2;
    let mut scores = vec![0_i32; cells];
    let mut trace = vec![DIAGONAL; cells];
    for row in 1..=left.len() {
        scores[row * columns] = row as i32 * GAP;
        trace[row * columns] = UP;
    }
    for column in 1..=right.len() {
        scores[column] = column as i32 * GAP;
        trace[column] = LEFT;
    }

    for row in 1..=left.len() {
        for column in 1..=right.len() {
            let index = row * columns + column;
            let diagonal = scores[(row - 1) * columns + column - 1]
                + if left[row - 1] == right[column - 1] {
                    2
                } else {
                    -1
                };
            let up = scores[(row - 1) * columns + column] + GAP;
            let left_score = scores[row * columns + column - 1] + GAP;
            let (score, direction) = if diagonal >= up && diagonal >= left_score {
                (diagonal, DIAGONAL)
            } else if up >= left_score {
                (up, UP)
            } else {
                (left_score, LEFT)
            };
            scores[index] = score;
            trace[index] = direction;
        }
    }

    let mut aligned_left = Vec::with_capacity(left.len().max(right.len()));
    let mut aligned_right = Vec::with_capacity(left.len().max(right.len()));
    let (mut row, mut column) = (left.len(), right.len());
    while row > 0 || column > 0 {
        let direction = trace[row * columns + column];
        if row > 0 && column > 0 && direction == DIAGONAL {
            aligned_left.push(left[row - 1]);
            aligned_right.push(right[column - 1]);
            row -= 1;
            column -= 1;
        } else if row > 0 && (column == 0 || direction == UP) {
            aligned_left.push(left[row - 1]);
            aligned_right.push(b'-');
            row -= 1;
        } else {
            aligned_left.push(b'-');
            aligned_right.push(right[column - 1]);
            column -= 1;
        }
    }
    aligned_left.reverse();
    aligned_right.reverse();
    Ok((
        String::from_utf8(aligned_left).expect("sequence alphabet is ASCII"),
        String::from_utf8(aligned_right).expect("sequence alphabet is ASCII"),
    ))
}

fn merge_center_alignments(
    master_center: &str,
    existing_rows: &[String],
    new_center: &str,
    new_row: &str,
) -> (String, Vec<String>, String) {
    let master = master_center.as_bytes();
    let new = new_center.as_bytes();
    let new_sequence = new_row.as_bytes();
    let existing: Vec<&[u8]> = existing_rows.iter().map(String::as_bytes).collect();
    let mut merged_center = Vec::new();
    let mut merged_existing = vec![Vec::new(); existing.len()];
    let mut merged_new = Vec::new();
    let (mut master_index, mut new_index) = (0, 0);

    while master_index < master.len() || new_index < new.len() {
        let master_base = master.get(master_index).copied();
        let new_base = new.get(new_index).copied();
        match (master_base, new_base) {
            (Some(b'-'), Some(b'-')) => {
                merged_center.push(b'-');
                for (output, row) in merged_existing.iter_mut().zip(&existing) {
                    output.push(row[master_index]);
                }
                merged_new.push(new_sequence[new_index]);
                master_index += 1;
                new_index += 1;
            }
            (Some(b'-'), _) => {
                merged_center.push(b'-');
                for (output, row) in merged_existing.iter_mut().zip(&existing) {
                    output.push(row[master_index]);
                }
                merged_new.push(b'-');
                master_index += 1;
            }
            (_, Some(b'-')) => {
                merged_center.push(b'-');
                for output in &mut merged_existing {
                    output.push(b'-');
                }
                merged_new.push(new_sequence[new_index]);
                new_index += 1;
            }
            (Some(master_residue), Some(new_residue)) => {
                debug_assert_eq!(master_residue, new_residue);
                merged_center.push(master_residue);
                for (output, row) in merged_existing.iter_mut().zip(&existing) {
                    output.push(row[master_index]);
                }
                merged_new.push(new_sequence[new_index]);
                master_index += 1;
                new_index += 1;
            }
            (Some(master_residue), None) => {
                merged_center.push(master_residue);
                for (output, row) in merged_existing.iter_mut().zip(&existing) {
                    output.push(row[master_index]);
                }
                merged_new.push(b'-');
                master_index += 1;
            }
            (None, Some(new_residue)) => {
                merged_center.push(new_residue);
                for output in &mut merged_existing {
                    output.push(b'-');
                }
                merged_new.push(new_sequence[new_index]);
                new_index += 1;
            }
            (None, None) => break,
        }
    }

    (
        String::from_utf8(merged_center).expect("sequence alphabet is ASCII"),
        merged_existing
            .into_iter()
            .map(|row| String::from_utf8(row).expect("sequence alphabet is ASCII"))
            .collect(),
        String::from_utf8(merged_new).expect("sequence alphabet is ASCII"),
    )
}

/// Exposed for a concise status message in the UI.
pub fn muscle_location() -> Option<OsString> {
    detect_muscle().map(PathBuf::into_os_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cleans_dna_and_protein_inputs() {
        assert_eq!(clean_alignment_sequence(">dna\nacgu n 12"), "ACGTN");
        assert_eq!(clean_alignment_sequence("MEEPQSDPSV"), "MEEPQSDPSV");
    }

    #[test]
    fn parses_wrapped_fasta_alignment() {
        let parsed = parse_fasta_alignment(">Seq1\nAC-\nGT\n>Seq2\nACT\nGT\n").unwrap();
        assert_eq!(parsed.labels, ["Seq1", "Seq2"]);
        assert_eq!(parsed.sequences, ["AC-GT", "ACTGT"]);
    }

    #[test]
    fn built_in_alignment_keeps_rows_equal_width() {
        let result = progressive_alignment(&[
            "ACGCTCGCT".to_owned(),
            "ACGTCGCT".to_owned(),
            "ACGCTTAGCT".to_owned(),
        ])
        .unwrap();
        assert_eq!(result.sequences.len(), 3);
        assert!(
            result
                .sequences
                .iter()
                .all(|row| row.len() == result.sequences[0].len())
        );
        for (source, aligned) in ["ACGCTCGCT", "ACGTCGCT", "ACGCTTAGCT"]
            .iter()
            .zip(&result.sequences)
        {
            assert_eq!(&aligned.replace('-', ""), source);
        }
    }
}
