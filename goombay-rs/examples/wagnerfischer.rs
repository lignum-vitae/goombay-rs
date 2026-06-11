use goombay_rs::align::{GlobalAlignmentMatrix, WagnerFischer};
use goombay_rs::scoring::LevenshteinScoring;

fn main() {
    // Sequences to be aligned
    let query = "attain";
    let subject = "atin";

    println!("Default Scoring and all alignments");
    // Use Default parameters
    let wf_default = WagnerFischer::compute(query, subject);
    // Align the sequences based on the pointer matrix
    let aligned = wf_default.all_alignments(true).align();
    println!("{}", wf_default.data.single_score_matrix());
    println!("{}", wf_default.data.single_pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = wf_default.similarity();
    let dist = wf_default.distance();
    let norm_sim = wf_default.normalized_similarity();
    let norm_dist = wf_default.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}\n"
    );

    println!("Custom Scoring and single alignment");
    // Set custom scoring parameters for Wagner Fischer
    // NOTE: Wagner Fischer is not meant to be used with custom scoring
    let scores = LevenshteinScoring {
        substitution: 3,
        gap: 2,
    };
    let nw_custom_scores = WagnerFischer::set_scores(&scores);
    let wf_custom = nw_custom_scores.calculate_matrix(query, subject);

    // Align the sequences based on the pointer matrix
    let aligned = wf_custom.align(); // One alignment returned by default
    println!("{}", wf_custom.data.single_score_matrix());
    println!("{}", wf_custom.data.single_pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = wf_custom.similarity();
    let dist = wf_custom.distance();
    let norm_sim = wf_custom.normalized_similarity();
    let norm_dist = wf_custom.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}"
    );
}
