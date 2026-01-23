use crate::align::global_base::{GlobalAlgorithm, GlobalAlignmentModel, Metric};
use crate::align::scoring::LevenshteinScoring;
use crate::align::{AlignmentData, GlobalAlignmentMatrix, PointerValues, Scoring};
use std::cmp::min;

pub struct WagnerFischer<S: Scoring + Clone> {
    pub scores: S,
}

impl Default for WagnerFischer<LevenshteinScoring> {
    fn default() -> Self {
        let scores = LevenshteinScoring {
            substitution: 1,
            gap: 1,
        };
        Self { scores }
    }
}

impl GlobalAlignmentMatrix<LevenshteinScoring> for WagnerFischer<LevenshteinScoring> {
    fn compute(query: &str, subject: &str) -> GlobalAlignmentModel {
        // Use default scores to calculate scoring and pointer matrices
        let wf_default = WagnerFischer::default();
        wf_default.calculate_matrix(query, subject)
    }
    fn set_scores(scores: &LevenshteinScoring) -> Self {
        // Set custom scores before manually calculating matrices
        Self {
            scores: scores.clone(),
        }
    }
    fn calculate_matrix(&self, query: &str, subject: &str) -> GlobalAlignmentModel {
        let mut alignments = AlignmentData::new(query, subject);
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;
        let score_matrix = &mut alignments.score_matrix[0];
        let pointer_matrix = &mut alignments.pointer_matrix[0];

        // initialise score and pointer matrices
        pointer_matrix[0][0] = PointerValues::Left as i32;
        for i in 1..query_len {
            score_matrix[i][0] = i as i32 * self.scores.gap as i32;
            pointer_matrix[i][0] = PointerValues::Up as i32;
        }
        for j in 1..subject_len {
            score_matrix[0][j] = j as i32 * self.scores.gap as i32;
            pointer_matrix[0][j] = PointerValues::Left as i32;
        }

        // Build pointer and score matrix
        for i in 1..query_len {
            for j in 1..subject_len {
                let identity = {
                    if alignments.query[i - 1] != alignments.subject[j - 1] {
                        score_matrix[i - 1][j - 1] + self.scores.substitution as i32
                    } else {
                        score_matrix[i - 1][j - 1] // Score unchanged for matching letters
                    }
                };
                let ugap = score_matrix[i - 1][j] + self.scores.gap as i32;
                let lgap = score_matrix[i][j - 1] + self.scores.gap as i32;

                let tmax = min(min(identity, ugap), lgap);
                score_matrix[i][j] = tmax;

                if tmax == identity {
                    pointer_matrix[i][j] += PointerValues::Match as i32;
                }
                if tmax == ugap {
                    pointer_matrix[i][j] += PointerValues::Up as i32;
                }
                if tmax == lgap {
                    pointer_matrix[i][j] += PointerValues::Left as i32;
                }
            }
        }

        GlobalAlignmentModel {
            data: alignments,
            aligner: GlobalAlgorithm::WagnerFischer,
            metric: Metric::Distance,
            identity: 0,
            mismatch: self.scores.substitution,
            gap: self.scores.gap,
            all_alignments: false,
        }
    }
}
