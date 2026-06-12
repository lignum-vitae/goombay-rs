use crate::align::global_base::{GlobalAlgorithm, GlobalAlignmentModel, Metric};
use crate::align::scoring::ExtendedGapScoring;
use crate::align::{AlignmentData, GlobalAlignmentMatrix, PointerValues, Scoring};

pub struct WatermanSmithBeyer<S: Scoring + Clone> {
    pub scores: S,
}

impl Default for WatermanSmithBeyer<ExtendedGapScoring> {
    fn default() -> Self {
        let scores = ExtendedGapScoring {
            identity: 1,
            mismatch: 1,
            gap: 3,
            extended_gap: 1,
        };
        Self { scores }
    }
}

/// Affine gap cost for a gap of length k:
/// cost = gap_open + (k * gap_extend)
fn gap_cost(gap_open: i32, gap_extend: i32, k: usize) -> i32 {
    gap_open + (k as i32 * gap_extend)
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
        // Use new_wsb to get 3 pointer matrices:
        // pointer_matrix[0] = direction pointers (Match/Up/Left)
        // pointer_matrix[1] = i_step (optimal up gap length)
        // pointer_matrix[2] = j_step (optimal left gap length)
        let mut alignments = AlignmentData::new_wsb(query, subject);
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;
        let score_matrix = &mut alignments.score_matrix[0];
        let (pointer_first, pointer_rest) = alignments.pointer_matrix.split_at_mut(1);
        let (i_step_matrix, j_step_matrix) = pointer_rest.split_at_mut(1);
        let pointer = &mut pointer_first[0];
        let i_step_pointer = &mut i_step_matrix[0];
        let j_step_pointer = &mut j_step_matrix[0];

        let gap_open = self.scores.gap as i32;
        let gap_extend = self.scores.extended_gap as i32;

        // Initialise first column (gaps in subject)
        pointer[0][0] = PointerValues::Left as i32;
        for i in 1..query_len {
            score_matrix[i][0] = -gap_open - (i as i32 * gap_extend);
            pointer[i][0] = PointerValues::Up as i32;
            i_step_pointer[i][0] = 1_i32; // i_step
        }

        // Initialise first row (gaps in query)
        for j in 1..subject_len {
            score_matrix[0][j] = -gap_open - (j as i32 * gap_extend);
            pointer[0][j] = PointerValues::Left as i32;
            j_step_pointer[0][j] = 1_i32; // j_step
        }

        // Fill scoring matrix
        for i in 1..query_len {
            for j in 1..subject_len {
                // Diagonal (match/mismatch)
                let match_score = if alignments.query[i - 1] == alignments.subject[j - 1] {
                    score_matrix[i - 1][j - 1] + self.scores.identity as i32
                } else {
                    score_matrix[i - 1][j - 1] - self.scores.mismatch as i32
                };

                // Up gap: consider all gap lengths k from 1 to i
                let mut ugap_score = i32::MIN;
                let mut u_step: usize = 1;
                for k in 1..=i {
                    let candidate = score_matrix[i - k][j] - gap_cost(gap_open, gap_extend, k);
                    if candidate > ugap_score {
                        ugap_score = candidate;
                        u_step = k;
                    }
                }

                // Left gap: consider all gap lengths k from 1 to j
                let mut lgap_score = i32::MIN;
                let mut l_step: usize = 1;
                for k in 1..=j {
                    let candidate = score_matrix[i][j - k] - gap_cost(gap_open, gap_extend, k);
                    if candidate > lgap_score {
                        lgap_score = candidate;
                        l_step = k;
                    }
                }

                let tmax = *[match_score, ugap_score, lgap_score].iter().max().unwrap();
                score_matrix[i][j] = tmax;

                // Set pointer values (cumulative, matching NW pattern)
                if tmax == match_score {
                    pointer[i][j] += PointerValues::Match as i32;
                }
                if tmax == ugap_score {
                    pointer[i][j] += PointerValues::Up as i32;
                    i_step_pointer[i][j] = u_step as i32;
                }
                if tmax == lgap_score {
                    pointer[i][j] += PointerValues::Left as i32;
                    j_step_pointer[i][j] = l_step as i32;
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
