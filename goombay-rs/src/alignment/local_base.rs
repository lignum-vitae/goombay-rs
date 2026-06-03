use crate::align::{AlignmentData, PointerValues};
use jedvek::Matrix2D;

#[derive(Clone)]
pub enum LocalAlgorithm {
    SmithWaterman,
}

// Handles matrices that store similarity score vs distance score
// For local alignment, we primarily use Similarity
#[derive(Clone)]
pub enum LocalMetric {
    Similarity,
}

pub struct LocalAlignmentModel {
    pub data: AlignmentData,
    pub aligner: LocalAlgorithm,
    pub metric: LocalMetric,
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub all_alignments: bool,
    pub max_score: i32,
    pub start_indices: Vec<(usize, usize)>, // Locations of max_score in the matrix
}

impl LocalAlignmentModel {
    pub fn all_alignments(&self, value: bool) -> Self {
        Self {
            data: self.data.clone(),
            aligner: self.aligner.clone(),
            metric: self.metric.clone(),
            identity: self.identity,
            mismatch: self.mismatch,
            gap: self.gap,
            all_alignments: value,
            max_score: self.max_score,
            start_indices: self.start_indices.clone(),
        }
    }

    fn select_aligner(&self) -> Box<dyn Iterator<Item = (String, String)> + '_> {
        match self.aligner {
            LocalAlgorithm::SmithWaterman => {
                let local_aligner = LocalAligner {
                    query_chars: &self.data.query,
                    subject_chars: &self.data.subject,
                    pointer_matrix: self.data.pointer_matrix(),
                    score_matrix: self.data.score_matrix(),
                    stack: self
                        .start_indices
                        .iter()
                        .map(|&(i, j)| (Vec::new(), Vec::new(), i, j))
                        .collect(),
                    all_alignments: self.all_alignments,
                };
                // Turns struct into dynamically dispatched iterator
                Box::new(local_aligner)
            }
        }
    }

    pub fn align(&self) -> Vec<String> {
        let iterator = self.select_aligner();
        let aligned_results: Vec<String> = iterator.map(|(qs, ss)| format!("{qs}\n{ss}")).collect();
        aligned_results
    }

    pub fn similarity(&self) -> i32 {
        self.max_score
    }
}

pub struct LocalAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Matrix2D<i32>,
    pub score_matrix: &'a Matrix2D<i32>,
    pub stack: Vec<(Vec<char>, Vec<char>, usize, usize)>,
    pub all_alignments: bool,
}

impl<'a> Iterator for LocalAligner<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        let identity = PointerValues::Match as i32;
        let up = PointerValues::Up as i32;
        let left = PointerValues::Left as i32;

        let identity_array = [
            identity,
            identity + up,
            identity + left,
            identity + up + left,
        ];
        let up_array = [up, up + identity, up + left, up + identity + left];
        let left_array = [left, left + identity, left + up, left + identity + left];

        while let Some((qs_align, ss_align, i, j)) = self.stack.pop() {
            // Local alignment stops when score reaches 0
            if self.score_matrix[i][j] == 0 {
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
