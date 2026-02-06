use crate::align::global_base::{GlobalAlgorithm, GlobalAlignmentModel, Metric};
use crate::align::scoring::ExtendedGapScoring;
use crate::align::{AlignmentData, GlobalAlignmentMatrix, PointerValues, Scoring};

pub struct Gotoh<S: Scoring + Clone> {
    pub scores: S,
}

impl Default for Gotoh<ExtendedGapScoring> {
    fn default() -> Self {
        let scores = ExtendedGapScoring {
            identity: 2,
            mismatch: 1,
            gap: 3,
            extended_gap: 1,
        };
        Self { scores }
    }
}

impl GlobalAlignmentMatrix<ExtendedGapScoring> for Gotoh<ExtendedGapScoring> {
    fn compute(query: &str, subject: &str) -> GlobalAlignmentModel {
        // Use default scores to calculate scoring and pointer matrices
        let gotoh_default = Gotoh::default();
        gotoh_default.calculate_matrix(query, subject)
    }
    
    fn set_scores(scores: &ExtendedGapScoring) -> Self {
        // Set custom scores before manually calculating matrices
        Self {
            scores: scores.clone(),
        }
    }
    
    fn calculate_matrix(&self, query: &str, subject: &str) -> GlobalAlignmentModel {
        let mut alignments = AlignmentData::new_gotoh(query, subject);
        let query_len = alignments.query.len() + 1;
        let subject_len = alignments.subject.len() + 1;
        
        // D, P, Q matrices correspond to score_matrix[0], score_matrix[1], score_matrix[2]
        let (score_first, score_rest) = alignments.score_matrix.split_at_mut(1);
        let d = &mut score_first[0];
        let (score_second, score_third) = score_rest.split_at_mut(1);
        let p = &mut score_second[0];
        let q = &mut score_third[0];
        
        let (ptr_first, ptr_rest) = alignments.pointer_matrix.split_at_mut(1);
        let d_ptr = &mut ptr_first[0];
        let (ptr_second, ptr_third) = ptr_rest.split_at_mut(1);
        let p_ptr = &mut ptr_second[0];
        let q_ptr = &mut ptr_third[0];
        
        // Initialize D matrix
        d[0][0] = 0;
        for i in 1..query_len {
            d[i][0] = -(self.scores.gap as i32 + (i as i32 * self.scores.extended_gap as i32));
        }
        for j in 1..subject_len {
            d[0][j] = -(self.scores.gap as i32 + (j as i32 * self.scores.extended_gap as i32));
        }
        
        // Initialize P matrix (gap extension in query)
        for i in 0..query_len {
            p[i][0] = 0;
        }
        
        // Initialize Q matrix (gap extension in subject)
        for j in 0..subject_len {
            q[0][j] = 0;
        }
        
        // Initialize pointer matrices
        for i in 1..query_len {
            d_ptr[i][0] = PointerValues::Up as i32;
        }
        for j in 1..subject_len {
            d_ptr[0][j] = PointerValues::Left as i32;
        }
        
        // Fill the matrices
        for i in 1..query_len {
            for j in 1..subject_len {
                // Calculate substitution score
                let substitution = if alignments.query[i - 1] == alignments.subject[j - 1] {
                    self.scores.identity as i32
                } else {
                    -(self.scores.mismatch as i32)
                };
                
                // Calculate D matrix (match/mismatch)
                let d_match = d[i-1][j-1] + substitution;
                let d_from_p = p[i-1][j-1] + substitution;
                let d_from_q = q[i-1][j-1] + substitution;
                d[i][j] = [d_match, d_from_p, d_from_q].iter().max().copied().unwrap_or(i32::MIN);
                
                // Calculate P matrix (gap extension in query)
                let p_gap_open = d[i-1][j] - (self.scores.gap as i32 + self.scores.extended_gap as i32);
                let p_gap_extend = p[i-1][j] - (self.scores.extended_gap as i32);
                p[i][j] = [p_gap_open, p_gap_extend].iter().max().copied().unwrap_or(i32::MIN);
                
                // Calculate Q matrix (gap extension in subject)
                let q_gap_open = d[i][j-1] - (self.scores.gap as i32 + self.scores.extended_gap as i32);
                let q_gap_extend = q[i][j-1] - (self.scores.extended_gap as i32);
                q[i][j] = [q_gap_open, q_gap_extend].iter().max().copied().unwrap_or(i32::MIN);
                
                // Set pointers for P matrix
                if p[i][j] == p_gap_open {
                    p_ptr[i][j] = PointerValues::Up as i32;
                    p_ptr[i-1][j] = PointerValues::Match as i32;
                } else if p[i][j] == p_gap_extend {
                    p_ptr[i][j] = PointerValues::Up as i32;
                    p_ptr[i-1][j] = PointerValues::Up as i32;
                }
                
                // Set pointers for Q matrix
                if q[i][j] == q_gap_open {
                    q_ptr[i][j] = PointerValues::Left as i32;
                    q_ptr[i][j-1] = PointerValues::Match as i32;
                } else if q[i][j] == q_gap_extend {
                    q_ptr[i][j] = PointerValues::Left as i32;
                    q_ptr[i][j-1] = PointerValues::Left as i32;
                }
                
                // Set pointers for D matrix based on which matrix gave the maximum score
                let max_score = [d[i][j], p[i][j], q[i][j]].iter().max().copied().unwrap_or(i32::MIN);
                
                if max_score == d[i][j] {
                    d_ptr[i][j] += PointerValues::Match as i32;
                }
                if max_score == p[i][j] {
                    d_ptr[i][j] += PointerValues::Up as i32;
                }
                if max_score == q[i][j] {
                    d_ptr[i][j] += PointerValues::Left as i32;
                }
            }
        }
        
        GlobalAlignmentModel {
            data: alignments,
            aligner: GlobalAlgorithm::Gotoh,
            metric: Metric::Similarity,
            identity: self.scores.identity,
            mismatch: self.scores.mismatch,
            gap: self.scores.gap,
            all_alignments: false,
        }
    }
}
