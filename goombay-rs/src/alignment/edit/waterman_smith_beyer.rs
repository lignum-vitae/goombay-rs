use crate::align::global_base::{GlobalAlgorithm, GlobalAlignmentModel, Metric};
use crate::align::scoring::ExtendedGapScoring;
use crate::align::{AlignmentData, GlobalAlignmentMatrix, PointerValues, Scoring};

pub struct WatermanSmithBeyer<S: Scoring + Clone> {
    pub scores: S,
}

impl Default for WatermanSmithBeyer<ExtendedGapScoring> {
    fn default() -> Self {
        let scores = ExtendedGapScoring {
            identity: 2,
            mismatch: 1,
            gap: 3,
            extended_gap: 1,
        };
        Self { scores }
    }
}

/// Affine gap cost for a gap of length k:
/// cost = gap_open + (k - 1) * gap_extend
fn gap_cost(gap_open: i32, gap_extend: i32, k: usize) -> i32 {
    gap_open + (k as i32 - 1) * gap_extend
}

impl GlobalAlignmentMatrix<ExtendedGapScoring> for WatermanSmithBeyer<ExtendedGapScoring> {
    fn compute(query: &str, subject: &str) -> GlobalAlignmentModel {
        let wsb_default = WatermanSmithBeyer::default();
        wsb_default.calculate_matrix(query, subject)
    }

    fn set_scores(scores: &ExtendedGapScoring) -> Self {
        Self {
            scores: scores.clone(),
        }
    }

    fn calculate_matrix(&self, query: &str, subject: &str) -> GlobalAlignmentModel {
        // Use new_gotoh to get 3 pointer matrices:
        // pointer_matrix[0] = direction pointers (Match/Up/Left)
        // pointer_matrix[1] = i_step (optimal up gap length)
        // pointer_matrix[2] = j_step (optimal left gap length)
        let mut alignments = AlignmentData::new_gotoh(query, subject);
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;

        let gap_open = self.scores.gap as i32;
        let gap_extend = self.scores.extended_gap as i32;

        // Initialise first column (gaps in subject)
        alignments.pointer_matrix[0][0][0] = PointerValues::Left as i32;
        for i in 1..query_len {
            alignments.score_matrix[0][i][0] = -gap_cost(gap_open, gap_extend, i);
            alignments.pointer_matrix[0][i][0] = PointerValues::Up as i32;
            alignments.pointer_matrix[1][i][0] = i as i32; // i_step
        }

        // Initialise first row (gaps in query)
        for j in 1..subject_len {
            alignments.score_matrix[0][0][j] = -gap_cost(gap_open, gap_extend, j);
            alignments.pointer_matrix[0][0][j] = PointerValues::Left as i32;
            alignments.pointer_matrix[2][0][j] = j as i32; // j_step
        }

        // Fill scoring matrix
        for i in 1..query_len {
            for j in 1..subject_len {
                // Diagonal (match/mismatch)
                let match_score = if alignments.query[i - 1] == alignments.subject[j - 1] {
                    alignments.score_matrix[0][i - 1][j - 1] + self.scores.identity as i32
                } else {
                    alignments.score_matrix[0][i - 1][j - 1] - self.scores.mismatch as i32
                };

                // Up gap: consider all gap lengths k from 1 to i
                let mut ugap_score = i32::MIN;
                let mut u_step: usize = 1;
                for k in 1..=i {
                    let candidate =
                        alignments.score_matrix[0][i - k][j] - gap_cost(gap_open, gap_extend, k);
                    if candidate > ugap_score {
                        ugap_score = candidate;
                        u_step = k;
                    }
                }

                // Left gap: consider all gap lengths k from 1 to j
                let mut lgap_score = i32::MIN;
                let mut l_step: usize = 1;
                for k in 1..=j {
                    let candidate =
                        alignments.score_matrix[0][i][j - k] - gap_cost(gap_open, gap_extend, k);
                    if candidate > lgap_score {
                        lgap_score = candidate;
                        l_step = k;
                    }
                }

                let tmax = *[match_score, ugap_score, lgap_score].iter().max().unwrap();
                alignments.score_matrix[0][i][j] = tmax;

                // Set pointer values (cumulative, matching NW pattern)
                if tmax == match_score {
                    alignments.pointer_matrix[0][i][j] += PointerValues::Match as i32;
                }
                if tmax == ugap_score {
                    alignments.pointer_matrix[0][i][j] += PointerValues::Up as i32;
                    alignments.pointer_matrix[1][i][j] = u_step as i32;
                }
                if tmax == lgap_score {
                    alignments.pointer_matrix[0][i][j] += PointerValues::Left as i32;
                    alignments.pointer_matrix[2][i][j] = l_step as i32;
                }
            }
        }

        GlobalAlignmentModel {
            data: alignments,
            aligner: GlobalAlgorithm::WatermanSmithBeyer,
            metric: Metric::Similarity,
            identity: self.scores.identity,
            mismatch: self.scores.mismatch,
            gap: self.scores.gap,
            all_alignments: false,
        }
    }
}
