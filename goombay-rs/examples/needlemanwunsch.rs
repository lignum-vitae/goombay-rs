use goombay_rs::align::{NeedlemanWunsch, Scoring};

fn main() {
    // Sequences to be aligned
    let query = "attain";
    let subject = "atin";

    // Use Default parameters
    let nw_default = NeedlemanWunsch::compute(query, subject);
    // Align the sequences based on the pointer matrix
    let aligned = nw_default.all_alignments(true).align();
    println!("{}", nw_default.data.score_matrix());
    println!("{}", nw_default.data.pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = nw_default.similarity();
    let dist = nw_default.distance();
    let norm_sim = nw_default.normalized_similarity();
    let norm_dist = nw_default.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}"
    );

    // Set custom scoring parameters for Needleman Wunsch
    let scores = Scoring {
        identity: 5,
        mismatch: 3,
        gap: 2,
        transpose: None,    // This is an Optional parameter
        extended_gap: None, // This is an Optional parameter
    };
    let nw_custom_scores = NeedlemanWunsch::set_scores(&scores);
    let nw_custom = nw_custom_scores.calculate_matrix(query, subject);

    // Align the sequences based on the pointer matrix
    let aligned = nw_custom.align(); // One alignment returned by default
    println!("{}", nw_custom.data.score_matrix());
    println!("{}", nw_custom.data.pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    // Note: Alignment scores can be calculated independently from alignment
    let sim = nw_custom.similarity();
    let dist = nw_custom.distance();
    let norm_sim = nw_custom.normalized_similarity();
    let norm_dist = nw_custom.normalized_distance();
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}"
    );
}
