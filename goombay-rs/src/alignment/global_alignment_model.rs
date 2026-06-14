use crate::align::{AlignmentData, TracebackState};
use crate::aligners::global::{GlobalAligner, GotohAligner, WsbAligner};

#[derive(Clone)]
pub enum GlobalAlgorithm {
    NeedlemanWunsch,
    WagnerFischer,
    WatermanSmithBeyer,
    Gotoh,
}

// Handles matrices that store similarity score vs distance score
#[derive(Clone)]
pub enum Metric {
    Similarity,
    Distance,
}

// This struct holds the user-facing alignment and scoring functions
pub struct GlobalAlignmentModel {
    pub data: AlignmentData,
    pub aligner: GlobalAlgorithm,
    pub metric: Metric,
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub all_alignments: bool,
}

impl GlobalAlignmentModel {
    pub fn all_alignments(&self, value: bool) -> Self {
        Self {
            data: self.data.clone(),
            aligner: self.aligner.clone(),
            metric: self.metric.clone(),
            identity: self.identity,
            mismatch: self.mismatch,
            gap: self.gap,
            all_alignments: value,
        }
    }

    fn select_aligner(&self) -> Box<dyn Iterator<Item = (String, String)> + '_> {
        let i = self.data.query.len();
        let j = self.data.subject.len();
        match self.aligner {
            GlobalAlgorithm::NeedlemanWunsch | GlobalAlgorithm::WagnerFischer => {
                let global_aligner = GlobalAligner {
                    query_chars: &self.data.query,
                    subject_chars: &self.data.subject,
                    pointer_matrix: self.data.single_pointer_matrix(),
                    stack: vec![(Vec::new(), Vec::new(), i, j)],
                    all_alignments: self.all_alignments,
                };
                Box::new(global_aligner)
            }
            GlobalAlgorithm::WatermanSmithBeyer => {
                let wsb_aligner = WsbAligner {
                    query_chars: &self.data.query,
                    subject_chars: &self.data.subject,
                    pointer_matrix: &self.data.pointer_matrix[0],
                    i_step_matrix: &self.data.pointer_matrix[1],
                    j_step_matrix: &self.data.pointer_matrix[2],
                    stack: vec![(Vec::new(), Vec::new(), i, j)],
                    all_alignments: self.all_alignments,
                };
                Box::new(wsb_aligner)
            }
            GlobalAlgorithm::Gotoh => {
                let gotoh_aligner = GotohAligner {
                    query_chars: &self.data.query,
                    subject_chars: &self.data.subject,
                    pointer_matrix: &self.data.pointer_matrix,
                    stack: vec![TracebackState {
                        query_seq: Vec::new(),
                        subject_seq: Vec::new(),
                        row: i,
                        col: j,
                        active_ptr_matrix: 0,
                    }],
                    all_alignments: self.all_alignments,
                };
                // Turns struct into dynamically dispatched iterator
                Box::new(gotoh_aligner)
            }
        }
    }

    pub fn align(&self) -> Vec<String> {
        let iterator = self.select_aligner();
        let aligned_results: Vec<String> = iterator.map(|(qs, ss)| format!("{qs}\n{ss}")).collect();
        aligned_results
    }

    pub fn similarity(&self) -> i32 {
        if self.data.query.is_empty() && self.data.subject.is_empty() {
            return 1;
        }
        match self.metric {
            Metric::Similarity => {
                if self.data.query.is_empty() && self.data.subject.is_empty() {
                    return 1;
                }
                let score_matrix = self.data.single_score_matrix();
                let i = self.data.query.len();
                let j = self.data.subject.len();
                score_matrix[i][j]
            }
            Metric::Distance => self
                .data
                .query
                .len()
                .max(self.data.subject.len())
                .saturating_sub(self.distance() as usize) as i32,
            // converting self.distance to usize is fine because
            // Lowrance Wagner and Wagner Fischer only return
            // positive values for `self.distance`.
        }
    }

    pub fn distance(&self) -> i32 {
        if self.data.query.is_empty() && self.data.subject.is_empty() {
            return 0;
        }
        match self.metric {
            Metric::Similarity => {
                if self.data.query.is_empty() && self.data.subject.is_empty() {
                    return 0;
                }
                let similarity = self.similarity();
                let max_possible =
                    self.data.query.len().max(self.data.subject.len()) * self.identity;
                if similarity > 0 {
                    return max_possible as i32 - similarity.abs();
                }
                similarity.abs()
            }
            Metric::Distance => {
                let score_matrix = self.data.single_score_matrix();
                let i = self.data.query.len();
                let j = self.data.subject.len();
                score_matrix[i][j]
            }
        }
    }

    pub fn normalized_similarity(&self) -> f64 {
        match self.metric {
            Metric::Similarity => {
                if self.data.query.is_empty() && self.data.subject.is_empty() {
                    return 1.0;
                }
                if self.data.query.is_empty() || self.data.subject.is_empty() {
                    return 0.0;
                }
                let raw_sim = (self.similarity()) as f64;
                let max_length = self.data.query.len().max(self.data.subject.len());
                let min_length = self.data.query.len().min(self.data.subject.len());
                let max_possible = (max_length * self.identity) as f64;
                let min_possible =
                    -((min_length * self.mismatch + (max_length - min_length) * self.gap) as f64);

                let score_range = max_possible - min_possible;
                if score_range.abs() < f64::EPSILON {
                    return 1.0;
                }
                (raw_sim - min_possible) / score_range
            }
            Metric::Distance => 1_f64 - self.normalized_distance(),
        }
    }

    pub fn normalized_distance(&self) -> f64 {
        match self.metric {
            Metric::Similarity => 1_f64 - self.normalized_similarity(),
            Metric::Distance => {
                // This distance metric is valid for Wagner-Fischer
                // May need to update if other algorithms are added
                let max_poss_dist = self.data.query.len().max(self.data.subject.len());
                self.distance() as f64 / max_poss_dist as f64
            }
        }
    }
}
