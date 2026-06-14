use crate::align::AlignmentData;
use crate::aligners::local::LocalAligner;

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
                    pointer_matrix: self.data.single_pointer_matrix(),
                    score_matrix: self.data.single_score_matrix(),
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
