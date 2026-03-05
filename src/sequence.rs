//! 配列解析 — FASTA解析、Needleman-Wunsch グローバルアラインメント、配列同一性
//!
//! タンパク質一次配列の入出力と比較機能。

use crate::amino::AminoAcid;

// ============================================================================
// FASTA 解析
// ============================================================================

/// FASTA レコード。
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// ヘッダー行 (先頭の '>' を除く)。
    pub header: String,
    /// アミノ酸配列。
    pub sequence: Vec<AminoAcid>,
}

/// 1文字コードから `AminoAcid` への変換。
#[must_use]
pub const fn aa_from_char(c: char) -> Option<AminoAcid> {
    match c {
        'A' => Some(AminoAcid::Ala),
        'R' => Some(AminoAcid::Arg),
        'N' => Some(AminoAcid::Asn),
        'D' => Some(AminoAcid::Asp),
        'C' => Some(AminoAcid::Cys),
        'Q' => Some(AminoAcid::Gln),
        'E' => Some(AminoAcid::Glu),
        'G' => Some(AminoAcid::Gly),
        'H' => Some(AminoAcid::His),
        'I' => Some(AminoAcid::Ile),
        'L' => Some(AminoAcid::Leu),
        'K' => Some(AminoAcid::Lys),
        'M' => Some(AminoAcid::Met),
        'F' => Some(AminoAcid::Phe),
        'P' => Some(AminoAcid::Pro),
        'S' => Some(AminoAcid::Ser),
        'T' => Some(AminoAcid::Thr),
        'W' => Some(AminoAcid::Trp),
        'Y' => Some(AminoAcid::Tyr),
        'V' => Some(AminoAcid::Val),
        _ => None,
    }
}

/// FASTA 形式のテキストを解析。
///
/// 不明な文字はスキップ。複数レコード対応。
#[must_use]
pub fn parse_fasta(input: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq: Vec<AminoAcid> = Vec::new();

    for line in input.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(hdr) = line.strip_prefix('>') {
            // 前のレコードを保存
            if let Some(header) = current_header.take() {
                records.push(FastaRecord {
                    header,
                    sequence: std::mem::take(&mut current_seq),
                });
            }
            current_header = Some(hdr.trim().to_string());
            current_seq.clear();
        } else if current_header.is_some() {
            for c in line.chars() {
                if let Some(aa) = aa_from_char(c.to_ascii_uppercase()) {
                    current_seq.push(aa);
                }
            }
        }
    }

    // 最後のレコード
    if let Some(header) = current_header {
        records.push(FastaRecord {
            header,
            sequence: current_seq,
        });
    }

    records
}

/// アミノ酸配列をFASTA文字列に変換。
#[must_use]
pub fn to_fasta(header: &str, sequence: &[AminoAcid]) -> String {
    let mut out = String::with_capacity(header.len() + sequence.len() + 10);
    out.push('>');
    out.push_str(header);
    out.push('\n');
    for (i, aa) in sequence.iter().enumerate() {
        out.push(aa.one_letter());
        // 80文字で改行 (FASTA慣例)
        if (i + 1) % 80 == 0 {
            out.push('\n');
        }
    }
    if !sequence.is_empty() && !sequence.len().is_multiple_of(80) {
        out.push('\n');
    }
    out
}

// ============================================================================
// Needleman-Wunsch グローバルアラインメント
// ============================================================================

/// アラインメント設定。
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    /// 一致スコア。
    pub match_score: i32,
    /// 不一致ペナルティ (負値)。
    pub mismatch_penalty: i32,
    /// ギャップペナルティ (負値)。
    pub gap_penalty: i32,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -1,
            gap_penalty: -2,
        }
    }
}

/// アラインメント結果。
#[derive(Debug, Clone)]
pub struct AlignResult {
    /// アラインメントスコア。
    pub score: i32,
    /// 一致数。
    pub matches: usize,
    /// 配列同一性 (0.0–1.0)。
    pub identity: f64,
    /// アラインメント長。
    pub alignment_length: usize,
}

/// Needleman-Wunsch グローバルアラインメント。
///
/// O(mn) 時間・空間。スコアと一致率を返す。
#[must_use]
pub fn needleman_wunsch(
    seq_a: &[AminoAcid],
    seq_b: &[AminoAcid],
    config: &AlignConfig,
) -> AlignResult {
    let m = seq_a.len();
    let n = seq_b.len();

    if m == 0 && n == 0 {
        return AlignResult {
            score: 0,
            matches: 0,
            identity: 1.0,
            alignment_length: 0,
        };
    }
    if m == 0 {
        return AlignResult {
            score: config.gap_penalty * n as i32,
            matches: 0,
            identity: 0.0,
            alignment_length: n,
        };
    }
    if n == 0 {
        return AlignResult {
            score: config.gap_penalty * m as i32,
            matches: 0,
            identity: 0.0,
            alignment_length: m,
        };
    }

    // DP テーブル
    let cols = n + 1;
    let mut dp = vec![0i32; (m + 1) * cols];

    // 初期化
    for i in 0..=m {
        dp[i * cols] = config.gap_penalty * i as i32;
    }
    for (j, val) in dp.iter_mut().enumerate().take(n + 1) {
        *val = config.gap_penalty * j as i32;
    }

    // 充填
    for i in 1..=m {
        for j in 1..=n {
            let s = if seq_a[i - 1] == seq_b[j - 1] {
                config.match_score
            } else {
                config.mismatch_penalty
            };
            let diag = dp[(i - 1) * cols + (j - 1)] + s;
            let up = dp[(i - 1) * cols + j] + config.gap_penalty;
            let left = dp[i * cols + (j - 1)] + config.gap_penalty;
            dp[i * cols + j] = diag.max(up).max(left);
        }
    }

    let score = dp[m * cols + n];

    // トレースバックで一致数・アラインメント長を算出
    let mut matches = 0usize;
    let mut align_len = 0usize;
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        align_len += 1;
        if i > 0 && j > 0 {
            let s = if seq_a[i - 1] == seq_b[j - 1] {
                config.match_score
            } else {
                config.mismatch_penalty
            };
            if dp[i * cols + j] == dp[(i - 1) * cols + (j - 1)] + s {
                if seq_a[i - 1] == seq_b[j - 1] {
                    matches += 1;
                }
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i * cols + j] == dp[(i - 1) * cols + j] + config.gap_penalty {
            i -= 1;
        } else {
            j -= 1;
        }
    }

    let identity = if align_len > 0 {
        matches as f64 / align_len as f64
    } else {
        1.0
    };

    AlignResult {
        score,
        matches,
        identity,
        alignment_length: align_len,
    }
}

/// 二配列の配列同一性 (0.0–1.0)。
///
/// 長さが同じ場合は位置ごとの比較。異なる場合はNW アラインメントで算出。
#[must_use]
pub fn sequence_identity(seq_a: &[AminoAcid], seq_b: &[AminoAcid]) -> f64 {
    if seq_a.is_empty() && seq_b.is_empty() {
        return 1.0;
    }
    if seq_a.len() == seq_b.len() {
        let matches = seq_a
            .iter()
            .zip(seq_b.iter())
            .filter(|(a, b)| a == b)
            .count();
        return matches as f64 / seq_a.len() as f64;
    }
    let result = needleman_wunsch(seq_a, seq_b, &AlignConfig::default());
    result.identity
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aa_from_char_all_20() {
        let codes = "ARNDCQEGHILKMFPSTWYV";
        for c in codes.chars() {
            assert!(aa_from_char(c).is_some(), "Failed for '{c}'");
        }
    }

    #[test]
    fn aa_from_char_invalid() {
        assert!(aa_from_char('X').is_none());
        assert!(aa_from_char('0').is_none());
    }

    #[test]
    fn parse_fasta_single_record() {
        let input = ">test protein\nARNDCQ\nEGHILK\n";
        let records = parse_fasta(input);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, "test protein");
        assert_eq!(records[0].sequence.len(), 12);
        assert_eq!(records[0].sequence[0], AminoAcid::Ala);
        assert_eq!(records[0].sequence[11], AminoAcid::Lys);
    }

    #[test]
    fn parse_fasta_multiple_records() {
        let input = ">seq1\nARND\n>seq2\nGHIL\n";
        let records = parse_fasta(input);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence.len(), 4);
        assert_eq!(records[1].sequence.len(), 4);
    }

    #[test]
    fn parse_fasta_empty() {
        let records = parse_fasta("");
        assert!(records.is_empty());
    }

    #[test]
    fn parse_fasta_ignores_unknown_chars() {
        let input = ">test\nAR-ND.CQ\n";
        let records = parse_fasta(input);
        assert_eq!(records[0].sequence.len(), 6);
    }

    #[test]
    fn to_fasta_roundtrip() {
        let seq = vec![AminoAcid::Ala, AminoAcid::Arg, AminoAcid::Asn];
        let fasta = to_fasta("test", &seq);
        let records = parse_fasta(&fasta);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, seq);
    }

    #[test]
    fn nw_identical_sequences() {
        let seq: Vec<AminoAcid> = vec![AminoAcid::Ala, AminoAcid::Arg, AminoAcid::Gly];
        let result = needleman_wunsch(&seq, &seq, &AlignConfig::default());
        assert_eq!(result.matches, 3);
        assert!((result.identity - 1.0).abs() < 1e-10);
        assert_eq!(result.score, 6); // 3 * match_score(2)
    }

    #[test]
    fn nw_completely_different() {
        let a = vec![AminoAcid::Ala, AminoAcid::Ala, AminoAcid::Ala];
        let b = vec![AminoAcid::Gly, AminoAcid::Gly, AminoAcid::Gly];
        let result = needleman_wunsch(&a, &b, &AlignConfig::default());
        assert_eq!(result.matches, 0);
        assert!((result.identity - 0.0).abs() < 1e-10);
    }

    #[test]
    fn nw_empty_sequences() {
        let result = needleman_wunsch(&[], &[], &AlignConfig::default());
        assert_eq!(result.score, 0);
        assert_eq!(result.matches, 0);
        assert!((result.identity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn nw_one_empty() {
        let a = vec![AminoAcid::Ala, AminoAcid::Arg];
        let result = needleman_wunsch(&a, &[], &AlignConfig::default());
        assert_eq!(result.score, -4); // 2 * gap_penalty(-2)
        assert_eq!(result.matches, 0);
    }

    #[test]
    fn nw_with_gaps() {
        // ARND vs ARD → 挿入ギャップ1つ
        let a = vec![
            AminoAcid::Ala,
            AminoAcid::Arg,
            AminoAcid::Asn,
            AminoAcid::Asp,
        ];
        let b = vec![AminoAcid::Ala, AminoAcid::Arg, AminoAcid::Asp];
        let result = needleman_wunsch(&a, &b, &AlignConfig::default());
        assert!(result.matches >= 3);
    }

    #[test]
    fn sequence_identity_identical() {
        let seq = vec![AminoAcid::Ala, AminoAcid::Arg];
        assert!((sequence_identity(&seq, &seq) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sequence_identity_half() {
        let a = vec![AminoAcid::Ala, AminoAcid::Arg];
        let b = vec![AminoAcid::Ala, AminoAcid::Gly];
        assert!((sequence_identity(&a, &b) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn sequence_identity_empty() {
        assert!((sequence_identity(&[], &[]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sequence_identity_different_lengths() {
        let a = vec![AminoAcid::Ala, AminoAcid::Arg, AminoAcid::Gly];
        let b = vec![AminoAcid::Ala, AminoAcid::Gly];
        let id = sequence_identity(&a, &b);
        assert!((0.0..=1.0).contains(&id));
    }

    #[test]
    fn align_config_default() {
        let cfg = AlignConfig::default();
        assert_eq!(cfg.match_score, 2);
        assert_eq!(cfg.mismatch_penalty, -1);
        assert_eq!(cfg.gap_penalty, -2);
    }

    #[test]
    fn parse_fasta_lowercase_input() {
        let input = ">test\narnd\n";
        let records = parse_fasta(input);
        assert_eq!(records[0].sequence.len(), 4);
    }
}
