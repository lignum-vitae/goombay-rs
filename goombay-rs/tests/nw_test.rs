use goombay_rs::align::{NeedlemanWunsch, Scoring};
use spindalis::utils::Arr2D;

#[test]
fn test_identical_sequences() {
    let nw = NeedlemanWunsch::compute("ACTG", "ACTG");

    let aligned = nw.align();
    let sim = nw.similarity();
    let dist = nw.distance();
    let norm_sim = nw.normalized_similarity();
    let norm_dist = nw.normalized_distance();

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * nw.identity) as i32); // alignment length = 4, identity score = 2
    assert_eq!(dist, 0);
    assert_eq!(norm_sim, 1.0);
    assert_eq!(norm_dist, 0.0);
}

#[test]
fn test_completely_different() {
    let nw = NeedlemanWunsch::compute("AAAA", "TTTT");
    let aligned = nw.align();
    let sim = nw.similarity();
    let dist = nw.distance();
    let norm_sim = nw.normalized_similarity();
    let norm_dist = nw.normalized_distance();

    assert_eq!(aligned[0], "AAAA\nTTTT");
    assert_eq!(sim, -4_i32 * nw.mismatch as i32); // mismatch score = 1
    assert_eq!(dist, (4 * nw.mismatch) as i32);
    assert_eq!(norm_sim, 0.0);
    assert_eq!(norm_dist, 1.0);
}

#[test]
fn test_different_length() {
    let test_cases = vec![
        ("ACTG", "ACT", "ACTG\nACT-"), // Longer query
        ("ACT", "ACTG", "ACT-\nACTG"), // Longer subject
        ("ACGT", "AGT", "ACGT\nA-GT"), // Internal gap
    ];
    for (query, subject, expected) in test_cases {
        let nw = NeedlemanWunsch::compute(query, subject);

        let aligned = nw.align();
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
        let nw = NeedlemanWunsch::compute(query, subject);
        let norm_sim = nw.normalized_similarity();
        let norm_dist = nw.normalized_distance();

        assert_eq!(norm_sim, expected_sim);
        assert_eq!(norm_dist, expected_dist);
    }
}

#[test]
fn test_empty_sequences() {
    let custom_scores = Scoring {
        identity: 2,
        mismatch: 1,
        gap: 1,
        transpose: None,
        extended_gap: None,
    };
    let custom_nw = NeedlemanWunsch::set_scores(&custom_scores);

    let gap_score = custom_scores.gap as i32;

    let test_cases = vec![
        ("", "ACTG", "----\nACTG", -4 * gap_score, 4 * gap_score),
        ("ACTG", "", "ACTG\n----", -4 * gap_score, 4 * gap_score),
        ("", "", "\n", 1, 0),
    ];

    for (query, subject, expected_align, expected_sim, expected_dist) in test_cases {
        let nw = custom_nw.calculate_matrix(query, subject);

        let aligned = nw.align();
        let sim = nw.similarity();
        let dist = nw.distance();

        assert_eq!(aligned[0], expected_align);
        assert_eq!(sim, expected_sim);
        assert_eq!(dist, expected_dist);
    }
}

#[test]
fn test_single_character() {
    let nw_match = NeedlemanWunsch::compute("A", "A");
    assert_eq!(nw_match.align()[0], "A\nA");
    assert_eq!(nw_match.similarity(), nw_match.identity as i32);
    assert_eq!(nw_match.distance(), 0);

    let nw_mismatch = NeedlemanWunsch::compute("A", "T");
    assert_eq!(nw_mismatch.align()[0], "A\nT");
    assert_eq!(nw_mismatch.similarity(), -(nw_mismatch.mismatch as i32));
    assert_eq!(nw_mismatch.distance(), nw_mismatch.mismatch as i32);
}

#[test]
fn test_case_sensitivity() {
    let test_cases = vec![("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")];

    for (query, subject) in test_cases {
        let nw_mixed = NeedlemanWunsch::compute(query, subject);

        let aligned_mixed = nw_mixed.align();

        let nw_upper = NeedlemanWunsch::compute(
            query.to_uppercase().as_str(),
            subject.to_uppercase().as_str(),
        );
        let aligned_upper = nw_upper.align();

        assert_eq!(aligned_mixed[0], aligned_upper[0]);

        let sim_mixed = nw_mixed.similarity();
        let sim_upper = nw_upper.similarity();
        assert!((sim_mixed - sim_upper).abs() == 0);
    }
}

#[test]
fn test_scoring_parameters() {
    let custom_scores = Scoring {
        identity: 1,
        mismatch: 2,
        gap: 3,
        transpose: None,
        extended_gap: None,
    };
    let custom_nw = NeedlemanWunsch::set_scores(&custom_scores);

    let nw_alignment = custom_nw.calculate_matrix("ACGT", "AGT");
    assert_eq!(nw_alignment.align()[0], "ACGT\nA-GT");

    let query = "AC";
    let subject = "AT";
    let nw_alignment_matrix = custom_nw.calculate_matrix(query, subject);

    let expected_score: Arr2D<i32> = Arr2D::from(&[[0, -3, -6], [-3, 1, -2], [-6, -2, -1]]);

    let score_matrix = &nw_alignment_matrix.data.score_matrix[0];
    for r in 0..=query.len() {
        for c in 0..=subject.len() {
            assert_eq!(
                score_matrix[r][c], expected_score[r][c],
                "Custom score matrix value mismatch at ({r}, {c})"
            );
        }
    }
}

#[test]
fn test_all_alignments() {
    let (query, subject) = ("ACCG", "ACG");
    let nw = NeedlemanWunsch::compute(query, subject);
    let all_aligned = nw.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ACCG\nAC-G".to_string()));
    assert!(all_aligned.contains(&"ACCG\nA-CG".to_string()));

    let (query, subject) = ("ATGTGTA", "ATA");
    let nw = NeedlemanWunsch::compute(query, subject);
    let all_aligned = nw.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 3);
    assert!(all_aligned.contains(&"ATGTGTA\nAT----A".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA----TA".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA--T--A".to_string()));
}
