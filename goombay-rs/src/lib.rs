pub mod alignment;
pub mod phylo;

pub mod scoring {
    pub use crate::alignment::scoring;

    pub use scoring::ExtendedGapScoring;
    pub use scoring::GeneralScoring;
    pub use scoring::LevenshteinScoring;
    pub use scoring::TransposeScoring;
}

pub mod align {
    // Re-exports everything from alignment folder as align
    pub use crate::alignment::*;

    pub use edit::gotoh::Gotoh;
    pub use edit::needleman_wunsch::NeedlemanWunsch;
    pub use edit::smith_waterman::SmithWaterman;
    pub use edit::wagner_fischer::WagnerFischer;
    pub use edit::waterman_smith_beyer::WatermanSmithBeyer;
}
