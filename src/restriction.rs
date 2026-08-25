//! Restriction enzyme recognition and cleavage-site analysis.

use rayon::prelude::*;
use std::collections::BTreeSet;
use std::sync::OnceLock;

const DATABASE: &str = include_str!("../assets/restriction_enzymes.tsv");

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RestrictionEnzyme {
    pub name: String,
    pub site: String,
    pub first_forward_cut: Option<i32>,
    pub first_reverse_cut: Option<i32>,
    pub second_forward_cut: Option<i32>,
    pub second_reverse_cut: Option<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RestrictionHit {
    pub name: String,
    pub site: String,
    /// One-based top-strand cleavage coordinates, matching Biopython semantics.
    pub positions: Vec<usize>,
}

pub fn enzymes() -> &'static [RestrictionEnzyme] {
    static ENZYMES: OnceLock<Vec<RestrictionEnzyme>> = OnceLock::new();
    ENZYMES.get_or_init(|| {
        DATABASE
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .filter_map(parse_enzyme)
            .collect()
    })
}

fn parse_enzyme(line: &str) -> Option<RestrictionEnzyme> {
    let mut fields = line.split('\t');
    Some(RestrictionEnzyme {
        name: fields.next()?.to_owned(),
        site: fields.next()?.to_owned(),
        first_forward_cut: parse_cut(fields.next()?),
        first_reverse_cut: parse_cut(fields.next()?),
        second_forward_cut: parse_cut(fields.next()?),
        second_reverse_cut: parse_cut(fields.next()?),
    })
}

fn parse_cut(value: &str) -> Option<i32> {
    (value != "-").then(|| value.parse().ok()).flatten()
}

/// Scan every known enzyme in parallel. Results include only enzymes with at
/// least one valid cleavage inside this linear sequence.
pub fn analyze_restriction_sites(sequence: &str) -> Vec<RestrictionHit> {
    let mut hits: Vec<_> = enzymes()
        .par_iter()
        .filter_map(|enzyme| scan_enzyme(sequence, enzyme))
        .collect();
    hits.sort_unstable_by(|left, right| left.name.cmp(&right.name));
    hits
}

fn scan_enzyme(sequence: &str, enzyme: &RestrictionEnzyme) -> Option<RestrictionHit> {
    let site = enzyme.site.as_bytes();
    if site.is_empty() || site.len() > sequence.len() || enzyme.first_forward_cut.is_none() {
        return None;
    }

    let mut positions = BTreeSet::new();
    scan_orientation(
        sequence.as_bytes(),
        site,
        &[enzyme.first_forward_cut, enzyme.second_forward_cut],
        false,
        &mut positions,
    );

    let reverse_site = reverse_complement_site(&enzyme.site);
    if reverse_site != enzyme.site {
        scan_orientation(
            sequence.as_bytes(),
            reverse_site.as_bytes(),
            &[enzyme.first_reverse_cut, enzyme.second_reverse_cut],
            true,
            &mut positions,
        );
    }

    (!positions.is_empty()).then(|| RestrictionHit {
        name: enzyme.name.clone(),
        site: enzyme.site.clone(),
        positions: positions.into_iter().collect(),
    })
}

fn scan_orientation(
    sequence: &[u8],
    site: &[u8],
    cuts: &[Option<i32>],
    reverse: bool,
    positions: &mut BTreeSet<usize>,
) {
    for (start, window) in sequence.windows(site.len()).enumerate() {
        if !window.iter().zip(site).all(|(&sequence_base, &site_base)| {
            iupac_mask(sequence_base) & iupac_mask(site_base) != 0
        }) {
            continue;
        }
        for cut in cuts.iter().flatten() {
            let coordinate = if reverse {
                start as i32 - cut + 1
            } else {
                start as i32 + cut + 1
            };
            if (1..=sequence.len() as i32).contains(&coordinate) {
                positions.insert(coordinate as usize);
            }
        }
    }
}

fn reverse_complement_site(site: &str) -> String {
    site.bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' | b'U' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'D' => 'H',
            b'H' => 'D',
            b'V' => 'B',
            _ => 'N',
        })
        .collect()
}

fn iupac_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enzyme(name: &str) -> &'static RestrictionEnzyme {
        enzymes().iter().find(|enzyme| enzyme.name == name).unwrap()
    }

    #[test]
    fn loads_full_biopython_database() {
        assert!(enzymes().len() >= 1_000);
        assert_eq!(enzyme("EcoRI").site, "GAATTC");
    }

    #[test]
    fn reports_correct_one_based_ecori_cut() {
        let hit = scan_enzyme("AAAAGAATTCTTT", enzyme("EcoRI")).unwrap();
        assert_eq!(hit.positions, vec![6]);
    }

    #[test]
    fn handles_type_iis_on_both_orientations() {
        let forward = scan_enzyme("AAAAAGGTCTCAAAAAAAAAAAA", enzyme("BsaI")).unwrap();
        assert_eq!(forward.positions, vec![13]);

        let reverse = scan_enzyme("AAAAAAAAAAGAGACCAAAAAAAAAAAAAAAAAAAA", enzyme("BsaI")).unwrap();
        assert_eq!(reverse.positions, vec![6]);
    }

    #[test]
    fn recognizes_iupac_site_symbols() {
        let synthetic = RestrictionEnzyme {
            name: "test".into(),
            site: "TCNGA".into(),
            first_forward_cut: Some(3),
            first_reverse_cut: Some(-3),
            second_forward_cut: None,
            second_reverse_cut: None,
        };
        assert!(scan_enzyme("TCCGA", &synthetic).is_some());
    }
}
