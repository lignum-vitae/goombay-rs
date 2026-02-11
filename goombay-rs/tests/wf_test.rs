use goombay_rs::align::{GlobalAlignmentMatrix, WagnerFischer};
use goombay_rs::scoring::LevenshteinScoring;
use jedvek::Matrix2D;

#[test]
fn test_identical_sequences() {
    let wf = WagnerFischer::compute("ACTG", "ACTG");

    let aligned = wf.align();
    let sim = wf.similarity();
    let dist = wf.distance();
    let norm_sim = wf.normalized_similarity();
    let norm_dist = wf.normalized_distance();

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, 4_i32);
    assert_eq!(dist, 0);
    assert_eq!(norm_sim, 1.0);
    assert_eq!(norm_dist, 0.0);
}

#[test]
fn test_completely_different() {
    let wf = WagnerFischer::compute("AAAA", "TTTT");
    let aligned = wf.align();
    let sim = wf.similarity();
    let dist = wf.distance();
    let norm_sim = wf.normalized_similarity();
    let norm_dist = wf.normalized_distance();

    assert_eq!(aligned[0], "AAAA\nTTTT");
    assert_eq!(sim, 0);
    assert_eq!(dist, (4 * wf.mismatch) as i32);
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
        let wf = WagnerFischer::compute(query, subject);

        let aligned = wf.align();
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
        let wf = WagnerFischer::compute(query, subject);
        let norm_sim = wf.normalized_similarity();
        let norm_dist = wf.normalized_distance();

        assert_eq!(norm_sim, expected_sim);
        assert_eq!(norm_dist, expected_dist);
    }
}

#[test]
fn test_empty_sequences() {
    let custom_scores = LevenshteinScoring {
        substitution: 1,
        gap: 5,
    };
    let custom_wf = WagnerFischer::set_scores(&custom_scores);

    let gap_score = custom_scores.gap as i32;

    let test_cases = vec![
        ("", "ACTG", "----\nACTG", 0, 4 * gap_score),
        ("ACTG", "", "ACTG\n----", 0, 4 * gap_score),
        ("", "", "\n", 1, 0),
    ];

    for (query, subject, expected_align, expected_sim, expected_dist) in test_cases {
        let wf = custom_wf.calculate_matrix(query, subject);

        let aligned = wf.align();
        let sim = wf.similarity();
        let dist = wf.distance();

        assert_eq!(aligned[0], expected_align);
        assert_eq!(sim, expected_sim);
        assert_eq!(dist, expected_dist);
    }
}

#[test]
fn test_single_character() {
    let nw_match = WagnerFischer::compute("A", "A");
    assert_eq!(nw_match.align()[0], "A\nA");
    assert_eq!(nw_match.similarity(), 1);
    assert_eq!(nw_match.distance(), 0);

    let nw_mismatch = WagnerFischer::compute("A", "T");
    assert_eq!(nw_mismatch.align()[0], "A\nT");
    assert_eq!(nw_mismatch.similarity(), 0);
    assert_eq!(nw_mismatch.distance(), nw_mismatch.mismatch as i32);
}

#[test]
fn test_case_sensitivity() {
    let test_cases = vec![("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")];

    for (query, subject) in test_cases {
        let nw_mixed = WagnerFischer::compute(query, subject);

        let aligned_mixed = nw_mixed.align();

        let nw_upper = WagnerFischer::compute(
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
    let custom_scores = LevenshteinScoring {
        substitution: 2,
        gap: 3,
    };
    let custom_wf = WagnerFischer::set_scores(&custom_scores);

    let wf_alignment = custom_wf.calculate_matrix("ACGT", "AGT");
    assert_eq!(wf_alignment.align()[0], "ACGT\nA-GT");

    let query = "AC";
    let subject = "AT";
    let wf_alignment_matrix = custom_wf.calculate_matrix(query, subject);

    let expected_score: Matrix2D<i32> = Matrix2D::from(&[[0, 3, 6], [3, 0, 3], [6, 3, 2]]);

    let score_matrix = &wf_alignment_matrix.data.score_matrix();
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
    let wf = WagnerFischer::compute(query, subject);
    let all_aligned = wf.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ACCG\nAC-G".to_string()));
    assert!(all_aligned.contains(&"ACCG\nA-CG".to_string()));

    let (query, subject) = ("ATGTGTA", "ATA");
    let wf = WagnerFischer::compute(query, subject);
    let all_aligned = wf.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 3);
    assert!(all_aligned.contains(&"ATGTGTA\nAT----A".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA----TA".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA--T--A".to_string()));
}
