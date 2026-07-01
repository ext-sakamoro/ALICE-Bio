//! `FASTA` (`RFC` de-facto) sequence file parser + minimal biological utilities.
//!
//! `FASTA` records begin with a `>` header line followed by one or more
//! sequence lines. This module returns one [`FastaRecord`] per header, with
//! whitespace stripped from the concatenated sequence body.
//!
//! Convenience helpers:
//!
//! - [`FastaRecord::gc_content`] — proportion of `G` + `C` bases.
//! - [`reverse_complement`] — Watson-Crick reverse complement for `DNA`.

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// `FASTA` parser errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastaError {
    /// A sequence line appeared before any header.
    OrphanSequence { line_number: usize },
    /// The header carried no identifier after the `>` marker.
    EmptyHeader { line_number: usize },
}

// ---------------------------------------------------------------------------
// FastaRecord
// ---------------------------------------------------------------------------

/// A single record in a `FASTA` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// The identifier (everything on the header line up to the first
    /// whitespace, without the leading `>`).
    pub id: String,
    /// The remainder of the header line after the identifier.
    pub description: String,
    /// The concatenated sequence with whitespace removed.
    pub sequence: String,
}

impl FastaRecord {
    /// Number of residues (nucleotides or amino acids) in the sequence.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the sequence is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// `GC` content — proportion of `G` and `C` bases, case-insensitive.
    /// Returns `0.0` for an empty sequence.
    #[must_use]
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc = self
            .sequence
            .chars()
            .filter(|c| matches!(c.to_ascii_uppercase(), 'G' | 'C'))
            .count();
        gc as f64 / self.sequence.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a `FASTA` text buffer into a vector of records.
///
/// # Errors
///
/// See [`FastaError`].
pub fn parse(text: &str) -> Result<Vec<FastaRecord>, FastaError> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_description = String::new();
    let mut current_sequence = String::new();

    for (index, raw_line) in text.lines().enumerate() {
        let line_number = index + 1;
        let line = raw_line.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(FastaRecord {
                    id,
                    description: std::mem::take(&mut current_description),
                    sequence: std::mem::take(&mut current_sequence),
                });
            }
            let (id_part, desc_part) = split_header(header);
            if id_part.is_empty() {
                return Err(FastaError::EmptyHeader { line_number });
            }
            current_id = Some(id_part.to_owned());
            current_description = desc_part.to_owned();
        } else {
            if current_id.is_none() {
                return Err(FastaError::OrphanSequence { line_number });
            }
            current_sequence.extend(line.chars().filter(|c| !c.is_whitespace()));
        }
    }
    if let Some(id) = current_id.take() {
        records.push(FastaRecord {
            id,
            description: current_description,
            sequence: current_sequence,
        });
    }
    Ok(records)
}

fn split_header(header: &str) -> (&str, &str) {
    let header = header.trim_start();
    header
        .split_once(char::is_whitespace)
        .map_or_else(|| (header, ""), |(id, rest)| (id, rest.trim()))
}

// ---------------------------------------------------------------------------
// DNA utilities
// ---------------------------------------------------------------------------

/// Watson-Crick reverse complement of a `DNA` sequence.
///
/// Unknown characters are preserved as-is; case is retained.
#[must_use]
pub fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_single_record() {
        let input = ">seq1 description\nACGT\n";
        let records = parse(input).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, "description");
        assert_eq!(records[0].sequence, "ACGT");
    }

    #[test]
    fn parses_multi_line_sequence() {
        let input = ">seq1\nACGT\nGGGG\nCCCC\n";
        let records = parse(input).unwrap();
        assert_eq!(records[0].sequence, "ACGTGGGGCCCC");
    }

    #[test]
    fn parses_multiple_records() {
        let input = ">a\nAAAA\n>b\nCCCC\n";
        let records = parse(input).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "a");
        assert_eq!(records[1].id, "b");
    }

    #[test]
    fn skips_blank_lines_between_records() {
        let input = ">a\nAAAA\n\n\n>b\nCCCC\n";
        let records = parse(input).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn empty_header_is_rejected() {
        let err = parse(">\nAAAA\n").unwrap_err();
        assert!(matches!(err, FastaError::EmptyHeader { .. }));
    }

    #[test]
    fn orphan_sequence_is_rejected() {
        let err = parse("AAAA\n>seq1\nCCCC\n").unwrap_err();
        assert!(matches!(err, FastaError::OrphanSequence { .. }));
    }

    #[test]
    fn gc_content_matches_expectation() {
        let rec = FastaRecord {
            id: "s".into(),
            description: String::new(),
            sequence: "GGCCAT".into(),
        };
        assert!((rec.gc_content() - 4.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn gc_content_is_zero_for_empty_sequence() {
        let rec = FastaRecord {
            id: "s".into(),
            description: String::new(),
            sequence: String::new(),
        };
        assert!(rec.gc_content().abs() < 1e-12);
    }

    #[test]
    fn reverse_complement_covers_all_bases() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("acgtN"), "Nacgt");
    }

    #[test]
    fn header_without_description_yields_empty_description() {
        let input = ">only_id\nACGT\n";
        let records = parse(input).unwrap();
        assert_eq!(records[0].id, "only_id");
        assert!(records[0].description.is_empty());
    }
}
