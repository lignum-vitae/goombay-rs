use crate::align::PointerValues;
use jedvek::Matrix2D;

/// Traceback aligner for WatermanSmithBeyer that handles variable-length gaps.
///
/// Unlike GlobalAligner which always steps by 1, WSB gaps can skip multiple
/// positions. The step sizes are stored in separate matrices.
pub struct WsbAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Matrix2D<i32>,
    pub i_step_matrix: &'a Matrix2D<i32>,
    pub j_step_matrix: &'a Matrix2D<i32>,
    pub stack: Vec<(Vec<char>, Vec<char>, usize, usize)>,
    pub all_alignments: bool,
}

impl<'a> Iterator for WsbAligner<'a> {
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
        let left_array = [left, left + identity, left + up, left + identity + up];

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

            // Diagonal: match/mismatch (step by 1 in both dimensions)
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

            // Up gap: step by i_step positions (gap in subject)
            if up_array.contains(&self.pointer_matrix[i][j]) {
                let step = self.i_step_matrix[i][j] as usize;
                let mut new_qs_align = qs_align.clone();
                let mut new_ss_align = ss_align.clone();
                for s in 0..step {
                    new_qs_align.push(self.query_chars[i - 1 - s]);
                    new_ss_align.push('-');
                }
                self.stack.push((new_qs_align, new_ss_align, i - step, j));
                if !self.all_alignments {
                    continue;
                }
            }

            // Left gap: step by j_step positions (gap in query)
            if left_array.contains(&self.pointer_matrix[i][j]) {
                let step = self.j_step_matrix[i][j] as usize;
                let mut new_qs_align = qs_align.clone();
                let mut new_ss_align = ss_align.clone();
                for s in 0..step {
                    new_qs_align.push('-');
                    new_ss_align.push(self.subject_chars[j - 1 - s]);
                }
                self.stack.push((new_qs_align, new_ss_align, i, j - step));
                if !self.all_alignments {
                    continue;
                }
            }
        }
        None
    }
}
