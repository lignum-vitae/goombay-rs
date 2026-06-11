use goombay_rs::align::{GlobalAlignmentMatrix, Gotoh};
use goombay_rs::scoring::ExtendedGapScoring;

#[test]
fn test_identical_sequences() {
    let gotoh = Gotoh::compute("ACTG", "ACTG");

    let aligned = gotoh.align();
    let sim = gotoh.similarity();
    let dist = gotoh.distance();
    let norm_sim = gotoh.normalized_similarity();
    let norm_dist = gotoh.normalized_distance();

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * gotoh.identity) as i32);
    assert_eq!(dist, 0);
    assert_eq!(norm_sim, 1.0);
    assert_eq!(norm_dist, 0.0);
}

#[test]
fn test_completely_different() {
    let gotoh = Gotoh::compute("AAAA", "TTTT");
    let aligned = gotoh.align();
    let sim = gotoh.similarity();
    let dist = gotoh.distance();
    let norm_sim = gotoh.normalized_similarity();
    let norm_dist = gotoh.normalized_distance();

    assert_eq!(aligned[0], "AAAA\nTTTT");
    assert_eq!(sim, -4_i32 * gotoh.mismatch as i32);
    assert_eq!(dist, (4 * gotoh.mismatch) as i32);
    assert_eq!(norm_sim, 0.0);
    assert_eq!(norm_dist, 1.0);
}

#[test]
fn test_different_length() {
    let test_cases = vec![
        ("ACTG", "ACT", "ACTG\nACT-"),
        ("ACT", "ACTG", "ACT-\nACTG"),
        ("ACGTAGTC", "ACAGTC", "ACGTAGTC\nAC--AGTC"),
    ];
    for (query, subject, expected) in test_cases {
        let gotoh = Gotoh::compute(query, subject);

        let aligned = gotoh.align();
        assert_eq!(aligned[0], expected);
    }
}

#[test]
fn test_normalisation() {
    let sequences = [
        ("ACTG", "BBBB"),
        ("ACTG", "ABBB"),
        ("ACTG", "ACBB"),
        ("ACTG", "ACTB"),
        ("ACTG", "ACTG"),
    ];
    let expected = [
        (0.0, 1.0),
        (0.25, 0.75),
        (0.5, 0.5),
        (0.75, 0.25),
        (1.0, 0.0),
    ];
    for ((query, subject), (expected_sim, expected_dist)) in sequences.iter().zip(expected) {
        let gotoh = Gotoh::compute(query, subject);
        let norm_sim = gotoh.normalized_similarity();
        let norm_dist = gotoh.normalized_distance();

        assert_eq!(norm_sim, expected_sim);
        assert_eq!(norm_dist, expected_dist);
    }
}

#[test]
fn test_empty_sequences() {
    let custom_scores = ExtendedGapScoring {
        identity: 2,
        mismatch: 1,
        gap: 1,          // Gap opening
        extended_gap: 1, // Gap extension
    };
    let custom_gotoh = Gotoh::set_scores(&custom_scores);

    let gap_open_score = custom_scores.gap as i32;
    let gap_ext_score = custom_scores.extended_gap as i32;

    // A gap of length 4 under affine scoring costs: gap_open + (len * gap_ext)
    let total_gap_penalty = gap_open_score + (4 * gap_ext_score);

    let test_cases = vec![
        (
            "",
            "ACTG",
            "----\nACTG",
            -total_gap_penalty,
            total_gap_penalty,
        ),
        (
            "ACTG",
            "",
            "ACTG\n----",
            -total_gap_penalty,
            total_gap_penalty,
        ),
        ("", "", "\n", 1, 0),
    ];

    for (query, subject, expected_align, expected_sim, expected_dist) in test_cases {
        let gotoh = custom_gotoh.calculate_matrix(query, subject);

        let aligned = gotoh.align();
        let sim = gotoh.similarity();
        let dist = gotoh.distance();

        assert_eq!(aligned[0], expected_align);
        assert_eq!(sim, expected_sim);
        assert_eq!(dist, expected_dist);
    }
}

#[test]
fn test_single_character() {
    let gotoh_match = Gotoh::compute("A", "A");
    assert_eq!(gotoh_match.align()[0], "A\nA");
    assert_eq!(gotoh_match.similarity(), gotoh_match.identity as i32);
    assert_eq!(gotoh_match.distance(), 0);

    let gotoh_mismatch = Gotoh::compute("A", "T");
    assert_eq!(gotoh_mismatch.align()[0], "A\nT");
    assert_eq!(
        gotoh_mismatch.similarity(),
        -(gotoh_mismatch.mismatch as i32)
    );
    assert_eq!(gotoh_mismatch.distance(), gotoh_mismatch.mismatch as i32);
}

#[test]
fn test_case_sensitivity() {
    let test_cases = vec![("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")];

    for (query, subject) in test_cases {
        let gotoh_mixed = Gotoh::compute(query, subject);
        let aligned_mixed = gotoh_mixed.align();

        let gotoh_upper = Gotoh::compute(
            query.to_uppercase().as_str(),
            subject.to_uppercase().as_str(),
        );
        let aligned_upper = gotoh_upper.align();

        assert_eq!(aligned_mixed[0], aligned_upper[0]);

        let sim_mixed = gotoh_mixed.similarity();
        let sim_upper = gotoh_upper.similarity();
        assert_eq!((sim_mixed - sim_upper).abs(), 0);
    }
}

#[test]
fn test_all_alignments() {
    let (query, subject) = ("ACCG", "ACG");
    let gotoh = Gotoh::compute(query, subject);
    let all_aligned = gotoh.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ACCG\nAC-G".to_string()));
    assert!(all_aligned.contains(&"ACCG\nA-CG".to_string()));

    let (query, subject) = ("ATGTGTA", "ATA");
    let gotoh = Gotoh::compute(query, subject);
    let all_aligned = gotoh.all_alignments(true).align();

    // Gotoh tracks affine blocks natively, optimizing choices cleanly
    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ATGTGTA\nAT----A".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA----TA".to_string()));
}
