use crate::align::global_alignment_model::GlobalAlignmentModel;
use crate::align::local_alignment_model::LocalAlignmentModel;
use jedvek::Matrix2D;
pub mod scoring;

pub mod aligners;
pub mod edit;
pub mod global_alignment_model;
pub mod local_alignment_model;

pub use scoring::Scoring;

pub enum PointerValues {
    Match = 2,
    Up = 3,
    Left = 4,
    Transpose = 8,
}

pub struct TracebackState {
    pub query_seq: Vec<char>,
    pub subject_seq: Vec<char>,
    pub row: usize,
    pub col: usize,
    pub active_ptr_matrix: usize, // e.g. for Gotoh , 0 for D, 1 for P, 2 for Q
}

pub trait GlobalAlignmentMatrix<S: Scoring + Clone> {
    fn compute(query: &str, subject: &str) -> GlobalAlignmentModel;
    fn set_scores(scores: &S) -> Self;
    fn calculate_matrix(&self, query: &str, subject: &str) -> GlobalAlignmentModel;
}

pub trait LocalAlignmentMatrix<S: Scoring + Clone> {
    fn compute(query: &str, subject: &str) -> LocalAlignmentModel;
    fn set_scores(scores: &S) -> Self;
    fn calculate_matrix(&self, query: &str, subject: &str) -> LocalAlignmentModel;
}

#[derive(Clone)]
pub struct AlignmentData {
    pub query: Vec<char>,
    pub subject: Vec<char>,
    pub score_matrix: Vec<Matrix2D<i32>>,
    pub pointer_matrix: Vec<Matrix2D<i32>>,
}

impl AlignmentData {
    pub fn new(query: &str, subject: &str) -> AlignmentData {
        let query: Vec<char> = query.to_uppercase().chars().collect();
        let subject: Vec<char> = subject.to_uppercase().chars().collect();
        let score_matrix = vec![Matrix2D::full(0, query.len() + 1, subject.len() + 1)];
        let pointer_matrix = vec![Matrix2D::full(0, query.len() + 1, subject.len() + 1)];
        AlignmentData {
            query,
            subject,
            score_matrix,
            pointer_matrix,
        }
    }

    pub fn new_gotoh(query: &str, subject: &str) -> AlignmentData {
        let query: Vec<char> = query.to_uppercase().chars().collect();
        let subject: Vec<char> = subject.to_uppercase().chars().collect();
        let score_matrix = vec![
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
        ];
        let pointer_matrix = vec![
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
            Matrix2D::full(3, query.len() + 1, subject.len() + 1),
            Matrix2D::full(4, query.len() + 1, subject.len() + 1),
        ];
        AlignmentData {
            query,
            subject,
            score_matrix,
            pointer_matrix,
        }
    }

    pub fn new_wsb(query: &str, subject: &str) -> AlignmentData {
        let query: Vec<char> = query.to_uppercase().chars().collect();
        let subject: Vec<char> = subject.to_uppercase().chars().collect();
        let score_matrix = vec![Matrix2D::full(0, query.len() + 1, subject.len() + 1)];
        let pointer_matrix = vec![
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
            Matrix2D::full(0, query.len() + 1, subject.len() + 1),
        ];
        AlignmentData {
            query,
            subject,
            score_matrix,
            pointer_matrix,
        }
    }

    pub fn single_score_matrix(&self) -> &Matrix2D<i32> {
        &self.score_matrix[0]
    }

    pub fn single_pointer_matrix(&self) -> &Matrix2D<i32> {
        &self.pointer_matrix[0]
    }
}
