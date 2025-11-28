use crate::align::base::GlobalAlignment;
use crate::align::{AlignmentData, PointerValues, Scoring};

pub struct NeedlemanWunsch {
    pub scores: Scoring,
}

impl NeedlemanWunsch {
    pub fn new(scores: &Scoring) -> Self {
        Self {
            scores: scores.clone(),
        }
    }

    pub fn build(&self, alignments: &mut AlignmentData) {
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;
        let score_matrix = &mut alignments.score_matrix[0];
        let pointer_matrix = &mut alignments.pointer_matrix[0];
        // initialise matrices
        pointer_matrix[0][0] = PointerValues::Left as i32;
        for i in 1..query_len {
            score_matrix[i][0] = -(i as i32 * self.scores.gap as i32);
            pointer_matrix[i][0] = PointerValues::Up as i32;
        }
        for j in 1..subject_len {
            score_matrix[0][j] = -(j as i32 * self.scores.gap as i32);
            pointer_matrix[0][j] = PointerValues::Left as i32;
        }

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

                let tmax = [identity, ugap, lgap].iter().max().copied().unwrap();
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
        alignments.is_built = true;
    }

    pub fn align(&self, alignments: &mut AlignmentData, all_alignments: bool) -> Vec<String> {
        if !alignments.is_built {
            self.build(alignments);
        }
        let i = alignments.query.len();
        let j = alignments.subject.len();

        let iterator = GlobalAlignment {
            query_chars: &alignments.query,
            subject_chars: &alignments.subject,
            pointer_matrix: alignments.pointer_matrix(),
            stack: vec![(Vec::new(), Vec::new(), i, j)],
            all_alignments,
            match_val: PointerValues::Match as i32,
            up_val: PointerValues::Up as i32,
            left_val: PointerValues::Left as i32,
        };
        let aligned_results: Vec<String> = iterator.map(|(qs, ss)| format!("{qs}\n{ss}")).collect();
        aligned_results
    }

    pub fn similarity(&self, alignments: &mut AlignmentData) -> i32 {
        if alignments.query.is_empty() && alignments.subject.is_empty() {
            return 1;
        }
        if !alignments.is_built {
            self.build(alignments);
        }
        let score_matrix = alignments.score_matrix();
        let i = alignments.query.len();
        let j = alignments.subject.len();
        score_matrix[i][j]
    }

    pub fn distance(&self, alignments: &mut AlignmentData) -> i32 {
        if alignments.query.is_empty() && alignments.subject.is_empty() {
            return 0;
        }
        if alignments.query.is_empty() || alignments.subject.is_empty() {
            let max_len = [alignments.query.len(), alignments.subject.len()]
                .iter()
                .max()
                .copied()
                .unwrap_or(0_usize);
            return (max_len * self.scores.mismatch) as i32;
        }
        if !alignments.is_built {
            self.build(alignments);
        }
        let similarity = self.similarity(alignments);
        let max_possible = [alignments.query.len(), alignments.subject.len()]
            .iter()
            .max()
            .copied()
            .unwrap()
            * self.scores.identity;
        max_possible as i32 - similarity.abs()
    }

    pub fn normalized_similarity(&self, alignments: &mut AlignmentData) -> f64 {
        let raw_sim = (self.similarity(alignments)) as f64;
        let max_length = [alignments.query.len(), alignments.subject.len()]
            .iter()
            .max()
            .copied()
            .unwrap();
        let max_possible = (max_length * self.scores.identity) as f64;
        let min_possible = (max_length * self.scores.mismatch) as f64;

        let score_range = max_possible + min_possible.abs();

        (raw_sim + min_possible.abs()) / score_range
    }

    pub fn normalized_distance(&self, alignments: &mut AlignmentData) -> f64 {
        1_f64 - self.normalized_similarity(alignments)
    }
}
