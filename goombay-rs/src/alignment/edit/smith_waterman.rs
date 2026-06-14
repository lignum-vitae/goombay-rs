use crate::align::local_alignment_model::{LocalAlgorithm, LocalAlignmentModel, LocalMetric};
use crate::align::scoring::GeneralScoring;
use crate::align::{AlignmentData, LocalAlignmentMatrix, PointerValues, Scoring};

pub struct SmithWaterman<S: Scoring + Clone> {
    pub scores: S,
}

impl Default for SmithWaterman<GeneralScoring> {
    fn default() -> Self {
        let scores = GeneralScoring {
            identity: 2,
            mismatch: 1,
            gap: 2,
        };
        Self { scores }
    }
}

impl LocalAlignmentMatrix<GeneralScoring> for SmithWaterman<GeneralScoring> {
    fn compute(query: &str, subject: &str) -> LocalAlignmentModel {
        let sw_default = SmithWaterman::default();
        sw_default.calculate_matrix(query, subject)
    }

    fn set_scores(scores: &GeneralScoring) -> Self {
        Self {
            scores: scores.clone(),
        }
    }

    fn calculate_matrix(&self, query: &str, subject: &str) -> LocalAlignmentModel {
        let mut alignments = AlignmentData::new(query, subject);
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;
        let score_matrix = &mut alignments.score_matrix[0];
        let pointer_matrix = &mut alignments.pointer_matrix[0];

        let mut max_score = 0;
        let mut start_indices = Vec::new();

        // initialise score and pointer matrices (first row and column are 0 for SW)
        // AlignmentData::new already initialises to 0, so we just set pointers if needed.
        // For local alignment, we don't strictly need to set pointers for the first row/col as we stop at 0.

        // Build pointer and score matrix
        for i in 1..query_len {
            for j in 1..subject_len {
                let identity = {
                    if alignments.query[i - 1] == alignments.subject[j - 1] {
                        score_matrix[i - 1][j - 1] + self.scores.identity as i32
                    } else {
                        score_matrix[i - 1][j - 1] - self.scores.mismatch as i32
                    }
                };
                let ugap = score_matrix[i - 1][j] - self.scores.gap as i32;
                let lgap = score_matrix[i][j - 1] - self.scores.gap as i32;

                // Smith-Waterman: score_matrix[i][j] = max(0, identity, ugap, lgap)
                let mut current_max = 0;
                let mut choices = Vec::new();

                if identity > current_max {
                    current_max = identity;
                    choices = vec![PointerValues::Match as i32];
                } else if identity == current_max && current_max > 0 {
                    choices.push(PointerValues::Match as i32);
                }

                if ugap > current_max {
                    current_max = ugap;
                    choices = vec![PointerValues::Up as i32];
                } else if ugap == current_max && current_max > 0 {
                    choices.push(PointerValues::Up as i32);
                }

                if lgap > current_max {
                    current_max = lgap;
                    choices = vec![PointerValues::Left as i32];
                } else if lgap == current_max && current_max > 0 {
                    choices.push(PointerValues::Left as i32);
                }

                score_matrix[i][j] = current_max;
                for choice in choices {
                    pointer_matrix[i][j] += choice;
                }

                if current_max > max_score {
                    max_score = current_max;
                    start_indices = vec![(i, j)];
                } else if current_max == max_score && max_score > 0 {
                    start_indices.push((i, j));
                }
            }
        }

        LocalAlignmentModel {
            data: alignments,
            aligner: LocalAlgorithm::SmithWaterman,
            metric: LocalMetric::Similarity,
            identity: self.scores.identity,
            mismatch: self.scores.mismatch,
            gap: self.scores.gap,
            all_alignments: false,
            max_score,
            start_indices,
        }
    }
}
