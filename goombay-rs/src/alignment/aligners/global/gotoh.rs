use crate::align::{PointerValues, TracebackState};
use jedvek::Matrix2D;

pub struct GotohAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Vec<Matrix2D<i32>>,
    pub stack: Vec<TracebackState>,
    pub all_alignments: bool,
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
