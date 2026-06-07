use goombay_rs::align::{LocalAlignmentMatrix, SmithWaterman};
use goombay_rs::scoring::GeneralScoring;

fn main() {
    // Sequences to be aligned
    let query = "TGTTACGGAAAAAAAAAAGTGAC";
    let subject = "GGTTGACTA";

    println!("Default Scoring and all alignments");
    // Use Default parameters
    let sw_default = SmithWaterman::compute(query, subject);
    // Align the sequences based on the pointer matrix
    let aligned = sw_default.all_alignments(true).align();
    println!("{}", sw_default.data.single_score_matrix());
    println!("{}", sw_default.data.single_pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = sw_default.similarity();
    println!("Similarity: {sim}");

    println!("Custom Scoring and single alignment");
    // Set custom scoring parameters for Smith Waterman
    let scores = GeneralScoring {
        identity: 5,
        mismatch: 3,
        gap: 2,
    };
    let sw_custom_scores = SmithWaterman::set_scores(&scores);
    let sw_custom = sw_custom_scores.calculate_matrix(query, subject);

    // Align the sequences based on the pointer matrix
    let aligned = sw_custom.align(); // One alignment returned by default
    println!("{}", sw_custom.data.single_score_matrix());
    println!("{}", sw_custom.data.single_pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = sw_custom.similarity();
    println!("Similarity: {sim}");
}
