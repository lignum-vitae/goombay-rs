use crate::align::PointerValues;
use jedvek::Matrix2D;

// This struct does the actual alignment
pub struct GlobalAligner<'a> {
    pub query_chars: &'a [char],
    pub subject_chars: &'a [char],
    pub pointer_matrix: &'a Matrix2D<i32>,
    pub stack: Vec<(Vec<char>, Vec<char>, usize, usize)>,
    pub all_alignments: bool,
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
