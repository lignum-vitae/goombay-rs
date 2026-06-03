use crate::align::global_base::GlobalAlignmentModel;
use crate::align::local_base::LocalAlignmentModel;
use jedvek::Matrix2D;
pub mod scoring;

pub mod edit;
pub mod global_base;
pub mod local_base;

pub use scoring::Scoring;

pub enum PointerValues {
    Match = 2,
    Up = 3,
    Left = 4,
    Transpose = 8,
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

    pub fn score_matrix(&self) -> &Matrix2D<i32> {
        &self.score_matrix[0]
    }

    pub fn pointer_matrix(&self) -> &Matrix2D<i32> {
        &self.pointer_matrix[0]
    }
}
