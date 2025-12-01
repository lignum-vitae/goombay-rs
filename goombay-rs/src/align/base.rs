use crate::align::{AlignmentData, PointerValues};
use spindalis::utils::Arr2D;

pub struct GloablAlignmentModel {
    pub data: AlignmentData,
    pub identity: usize,
    pub mismatch: usize,
    pub all_alignments: bool,
}

impl GloablAlignmentModel {
    pub fn all_alignments(&self, value: bool) -> Self {
        Self {
            data: self.data.clone(),
            identity: self.identity,
            mismatch: self.mismatch,
            all_alignments: value,
        }
    }
    pub fn align(&self) -> Vec<String> {
        let i = self.data.query.len();
        let j = self.data.subject.len();

        let iterator = GlobalAligner {
            query_chars: &self.data.query,
            subject_chars: &self.data.subject,
            pointer_matrix: self.data.pointer_matrix(),
            stack: vec![(Vec::new(), Vec::new(), i, j)],
            all_alignments: self.all_alignments,
            match_val: PointerValues::Match as i32,
            up_val: PointerValues::Up as i32,
            left_val: PointerValues::Left as i32,
        };
        let aligned_results: Vec<String> = iterator.map(|(qs, ss)| format!("{qs}\n{ss}")).collect();
        aligned_results
    }

    pub fn similarity(&self) -> i32 {
        if self.data.query.is_empty() && self.data.subject.is_empty() {
            return 1;
        }
        let score_matrix = self.data.score_matrix();
        let i = self.data.query.len();
        let j = self.data.subject.len();
        score_matrix[i][j]
    }

    pub fn distance(&self) -> i32 {
        if self.data.query.is_empty() && self.data.subject.is_empty() {
            return 0;
        }
        if self.data.query.is_empty() || self.data.subject.is_empty() {
            let max_len = [self.data.query.len(), self.data.subject.len()]
                .iter()
                .max()
                .copied()
                .unwrap_or(0_usize);
            return (max_len * self.mismatch) as i32;
        }
        let similarity = self.similarity();
        let max_possible = [self.data.query.len(), self.data.subject.len()]
            .iter()
            .max()
            .copied()
            .unwrap()
            * self.identity;
        max_possible as i32 - similarity.abs()
    }

    pub fn normalized_similarity(&self) -> f64 {
        let raw_sim = (self.similarity()) as f64;
        let max_length = [self.data.query.len(), self.data.subject.len()]
            .iter()
            .max()
            .copied()
            .unwrap();
        let max_possible = (max_length * self.identity) as f64;
        let min_possible = (max_length * self.mismatch) as f64;

        let score_range = max_possible + min_possible.abs();

        (raw_sim + min_possible.abs()) / score_range
    }

    pub fn normalized_distance(&self) -> f64 {
        1_f64 - self.normalized_similarity()
    }
}

pub struct GlobalAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Arr2D<i32>,
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
