use spindalis::utils::Arr2D;

pub mod base;
pub mod edit;

pub use edit::NeedlemanWunsch;

pub enum PointerValues {
    Match = 2,
    Up = 3,
    Left = 4,
    Transpose = 8,
}

#[derive(Clone)]
pub struct Scoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub transpose: Option<usize>,
    pub extended_gap: Option<usize>,
}

pub struct AlignmentData {
    pub query: Vec<char>,
    pub subject: Vec<char>,
    pub score_matrix: Vec<Arr2D<i32>>,
    pub pointer_matrix: Vec<Arr2D<i32>>,
    pub is_built: bool,
}

impl AlignmentData {
    pub fn new(query: &str, subject: &str) -> AlignmentData {
        let query: Vec<char> = query.to_uppercase().chars().collect();
        let subject: Vec<char> = subject.to_uppercase().chars().collect();
        let score_matrix = vec![Arr2D::full(0, query.len() + 1, subject.len() + 1)];
        let pointer_matrix = vec![Arr2D::full(0, query.len() + 1, subject.len() + 1)];
        AlignmentData {
            query,
            subject,
            score_matrix,
            pointer_matrix,
            is_built: false,
        }
    }

    pub fn new_gotoh(query: &str, subject: &str) -> AlignmentData {
        let query: Vec<char> = query.to_uppercase().chars().collect();
        let subject: Vec<char> = subject.to_uppercase().chars().collect();
        let score_matrix = vec![
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
        ];
        let pointer_matrix = vec![
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
            Arr2D::full(0, query.len() + 1, subject.len() + 1),
        ];
        AlignmentData {
            query,
            subject,
            score_matrix,
            pointer_matrix,
            is_built: false,
        }
    }

    pub fn score_matrix(&self) -> &Arr2D<i32> {
        &self.score_matrix[0]
    }

    pub fn pointer_matrix(&self) -> &Arr2D<i32> {
        &self.pointer_matrix[0]
    }
}
