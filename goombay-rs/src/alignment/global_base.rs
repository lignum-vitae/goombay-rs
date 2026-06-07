use crate::align::{AlignmentData, PointerValues, TracebackState};
use jedvek::Matrix2D;

#[derive(Clone)]
pub enum GlobalAlgorithm {
    NeedlemanWunsch,
    WagnerFischer,
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
        match self.aligner {
            GlobalAlgorithm::NeedlemanWunsch | GlobalAlgorithm::WagnerFischer => {
                let i = self.data.query.len();
                let j = self.data.subject.len();
                let global_aligner = GlobalAligner {
                    query_chars: &self.data.query,
                    subject_chars: &self.data.subject,
                    pointer_matrix: self.data.single_pointer_matrix(),
                    stack: vec![(Vec::new(), Vec::new(), i, j)],
                    all_alignments: self.all_alignments,
                    match_val: PointerValues::Match as i32,
                    up_val: PointerValues::Up as i32,
                    left_val: PointerValues::Left as i32,
                };
                // Turns struct into dynamically dispatched iterator
                Box::new(global_aligner)
            }
            GlobalAlgorithm::Gotoh => {
                let i = self.data.query.len();
                let j = self.data.subject.len();
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
                    match_val: PointerValues::Match as i32,
                    up_val: PointerValues::Up as i32,
                    left_val: PointerValues::Left as i32,
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

// This struct does the actual alignment
pub struct GlobalAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Matrix2D<i32>,
    pub stack: Vec<(Vec<char>, Vec<char>, usize, usize)>,
    pub all_alignments: bool,
    pub match_val: i32,
    pub up_val: i32,
    pub left_val: i32,
}

impl<'a> Iterator for GlobalAligner<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        let identity = PointerValues::Match as i32; // 2
        let up = PointerValues::Up as i32; // 3
        let left = PointerValues::Left as i32; // 4

        let identity_array = [
            identity,
            identity + up,
            identity + left,
            identity + up + left,
        ];
        let left_array = [left, left + identity, left + up, left + identity + up];
        let up_array = [up, up + identity, up + left, up + identity + left];

        while let Some((qs_align, ss_align, i, j)) = self.stack.pop() {
            if i == 0 && j == 0 {
                let mut qs_align = qs_align;
                let mut ss_align = ss_align;
                qs_align.reverse();
                ss_align.reverse();
                let qs_aligned = qs_align.into_iter().collect::<String>();
                let ss_aligned = ss_align.into_iter().collect::<String>();

                if !self.all_alignments {
                    self.stack.clear();
                }
                return Some((qs_aligned, ss_aligned));
            }

            if identity_array.contains(&self.pointer_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push(self.query_chars[i - 1]);
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push(self.subject_chars[j - 1]);
                self.stack.push((new_qs_align, new_ss_align, i - 1, j - 1));
                if !self.all_alignments {
                    continue;
                }
            }

            if up_array.contains(&self.pointer_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push(self.query_chars[i - 1]);
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push('-');

                self.stack.push((new_qs_align, new_ss_align, i - 1, j));
                if !self.all_alignments {
                    continue;
                }
            }

            if left_array.contains(&self.pointer_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push('-');
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push(self.subject_chars[j - 1]);
                self.stack.push((new_qs_align, new_ss_align, i, j - 1));
                if !self.all_alignments {
                    continue;
                }
            }
        }
        None
    }
}

pub struct GotohAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Vec<Matrix2D<i32>>,
    pub stack: Vec<TracebackState>,
    pub all_alignments: bool,
    pub match_val: i32,
    pub up_val: i32,
    pub left_val: i32,
}

impl<'a> Iterator for GotohAligner<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        let identity = PointerValues::Match as i32; // 2
        let up = PointerValues::Up as i32; // 3
        let left = PointerValues::Left as i32; // 4

        let identity_array = [
            identity,
            identity + up,
            identity + left,
            identity + up + left,
        ];
        let left_array = [left, left + identity, left + up, left + identity + up];
        let up_array = [up, up + identity, up + left, up + identity + left];

        let (d_ptr_idx, p_ptr_idx, q_ptr_idx) = (0, 1, 2);
        while let Some(state) = self.stack.pop() {
            let TracebackState {
                query_seq: qs_align,
                subject_seq: ss_align,
                row: i,
                col: j,
                active_ptr_matrix: ptr_matrix_idx,
            } = state;
            let active_matrix = &self.pointer_matrix[ptr_matrix_idx];
            if i == 0 && j == 0 {
                let mut qs_align = qs_align;
                let mut ss_align = ss_align;
                qs_align.reverse();
                ss_align.reverse();
                let qs_aligned = qs_align.into_iter().collect::<String>();
                let ss_aligned = ss_align.into_iter().collect::<String>();

                if !self.all_alignments {
                    self.stack.clear();
                }
                return Some((qs_aligned, ss_aligned));
            }

            if identity_array.contains(&active_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push(self.query_chars[i - 1]);
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push(self.subject_chars[j - 1]);
                self.stack.push(TracebackState {
                    query_seq: new_qs_align,
                    subject_seq: new_ss_align,
                    row: i - 1,
                    col: j - 1,
                    active_ptr_matrix: d_ptr_idx,
                });
                if !self.all_alignments {
                    continue;
                }
            }

            if up_array.contains(&active_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push(self.query_chars[i - 1]);
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push('-');
                self.stack.push(TracebackState {
                    query_seq: new_qs_align,
                    subject_seq: new_ss_align,
                    row: i - 1,
                    col: j,
                    active_ptr_matrix: p_ptr_idx,
                });
                if !self.all_alignments {
                    continue;
                }
            }

            if left_array.contains(&active_matrix[i][j]) {
                let mut new_qs_align = qs_align.clone();
                new_qs_align.push('-');
                let mut new_ss_align = ss_align.clone();
                new_ss_align.push(self.subject_chars[j - 1]);
                self.stack.push(TracebackState {
                    query_seq: new_qs_align,
                    subject_seq: new_ss_align,
                    row: i,
                    col: j - 1,
                    active_ptr_matrix: q_ptr_idx,
                });
                if !self.all_alignments {
                    continue;
                }
            }
        }
        None
    }
}
