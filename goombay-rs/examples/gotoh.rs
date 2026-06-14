use goombay_rs::align::{GlobalAlignmentMatrix, Gotoh};
use goombay_rs::scoring::ExtendedGapScoring;

fn main() {
    // Sequences to be aligned
    let query = "ATGA";
    let subject = "ACCCCCCCCTGCCCCCCCCA";

    println!("Default Scoring and all alignments");
    // Use Default parameters
    let gotoh_default = Gotoh::compute(query, subject);
    // Align the sequences based on the pointer matrix
    let aligned = gotoh_default.all_alignments(true).align();
    println!("{}", gotoh_default.data.single_score_matrix());
    println!("{}", gotoh_default.data.single_pointer_matrix());
    println!("\nSingle alignment");
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = gotoh_default.similarity();
    let dist = gotoh_default.distance();
    let norm_sim = gotoh_default.normalized_similarity();
    let norm_dist = gotoh_default.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}\n"
    );

    let query = "ATCA";
    let subject = "ACCCCCCCCTGCCCCCCCCA";
    let gotoh_default = Gotoh::compute(query, subject);
    let aligned = gotoh_default.all_alignments(true).align();
    println!("\n Multiple alignments");
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    println!("\nCustom Scoring and single alignment");
    // Set custom scoring parameters for Needleman Wunsch
    let scores = ExtendedGapScoring {
        identity: 5,
        mismatch: 3,
        gap: 2,
        extended_gap: 1,
    };
    let gotoh_custom_scores = Gotoh::set_scores(&scores);
    let gotoh_custom = gotoh_custom_scores.calculate_matrix(query, subject);

    // Align the sequences based on the pointer matrix
    let aligned = gotoh_custom.align(); // One alignment returned by default
    println!("{}", gotoh_custom.data.single_score_matrix());
    println!("{}", gotoh_custom.data.single_pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = gotoh_custom.similarity();
    let dist = gotoh_custom.distance();
    let norm_sim = gotoh_custom.normalized_similarity();
    let norm_dist = gotoh_custom.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}"
    );
}
