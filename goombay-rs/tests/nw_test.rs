use goombay_rs::align::{AlignmentData, NeedlemanWunsch, Scoring};
use spindalis::utils::Arr2D;

fn setup_nw() -> NeedlemanWunsch {
    let scores = Scoring {
        identity: 2,
        mismatch: 1,
        gap: 1,
        transpose: None,
        extended_gap: None,
    };
    NeedlemanWunsch::new(&scores)
}

#[test]
fn test_identical_sequences() {
    let nw = setup_nw();
    let mut nw_alignment = AlignmentData::new("ACTG", "ACTG");

    let aligned = nw.align(&mut nw_alignment, false);
    let sim = nw.similarity(&mut nw_alignment);
    let dist = nw.distance(&mut nw_alignment);
    let norm_sim = nw.normalized_similarity(&mut nw_alignment);
    let norm_dist = nw.normalized_distance(&mut nw_alignment);

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * nw.scores.identity) as i32); // alignment length = 4, identity score = 2
    assert_eq!(dist, 0);
    assert_eq!(norm_sim, 1.0);
    assert_eq!(norm_dist, 0.0);
}

#[test]
fn test_completely_different() {
    let nw = setup_nw();
    let mut nw_alignment = AlignmentData::new("AAAA", "TTTT");
    nw.build(&mut nw_alignment);
    let aligned = nw.align(&mut nw_alignment, false);
    let sim = nw.similarity(&mut nw_alignment);
    let dist = nw.distance(&mut nw_alignment);
    let norm_sim = nw.normalized_similarity(&mut nw_alignment);
    let norm_dist = nw.normalized_distance(&mut nw_alignment);

    assert_eq!(aligned[0], "AAAA\nTTTT");
    assert_eq!(sim, -4_i32 * nw.scores.mismatch as i32); // mismatch score = 1
    assert_eq!(dist, (4 * nw.scores.mismatch) as i32);
    assert_eq!(norm_sim, 0.0);
    assert_eq!(norm_dist, 1.0);
}

#[test]
fn test_different_length() {
    let nw = setup_nw();
    let test_cases = vec![
        ("ACTG", "ACT", "ACTG\nACT-"), // Longer query
        ("ACT", "ACTG", "ACT-\nACTG"), // Longer subject
        ("ACGT", "AGT", "ACGT\nA-GT"), // Internal gap
    ];
    for (query, subject, expected) in test_cases {
        let mut nw_alignment = AlignmentData::new(query, subject);

        let aligned = nw.align(&mut nw_alignment, false);
        assert_eq!(aligned[0], expected);
    }
}

#[test]
fn test_normalisation() {
    let nw = setup_nw();
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
        let mut nw_alignment = AlignmentData::new(query, subject);
        nw.build(&mut nw_alignment);
        let norm_sim = nw.normalized_similarity(&mut nw_alignment);
        let norm_dist = nw.normalized_distance(&mut nw_alignment);

        assert_eq!(norm_sim, expected_sim);
        assert_eq!(norm_dist, expected_dist);
    }
}

#[test]
fn test_empty_sequences() {
    let nw = setup_nw();
    let gap_score = nw.scores.gap as i32;

    let test_cases = vec![
        ("", "ACTG", "----\nACTG", -4 * gap_score, 4 * gap_score),
        ("ACTG", "", "ACTG\n----", -4 * gap_score, 4 * gap_score),
        ("", "", "\n", 1, 0),
    ];

    for (query, subject, expected_align, expected_sim, expected_dist) in test_cases {
        let mut nw_alignment = AlignmentData::new(query, subject);

        let aligned = nw.align(&mut nw_alignment, false);
        let sim = nw.similarity(&mut nw_alignment);
        let dist = nw.distance(&mut nw_alignment);

        assert_eq!(aligned[0], expected_align);
        assert_eq!(sim, expected_sim);
        assert_eq!(dist, expected_dist);
    }
}

#[test]
fn test_single_character() {
    let nw = setup_nw();
    let scores = &nw.scores;

    let mut nw_match = AlignmentData::new("A", "A");
    nw.build(&mut nw_match);
    assert_eq!(nw.align(&mut nw_match, false)[0], "A\nA");
    assert_eq!(nw.similarity(&mut nw_match), scores.identity as i32);
    assert_eq!(nw.distance(&mut nw_match), 0);

    let mut nw_mismatch = AlignmentData::new("A", "T");
    nw.build(&mut nw_mismatch);
    assert_eq!(nw.align(&mut nw_mismatch, false)[0], "A\nT");
    assert_eq!(nw.similarity(&mut nw_mismatch), -(scores.mismatch as i32));
    assert_eq!(nw.distance(&mut nw_mismatch), scores.mismatch as i32);
}

#[test]
fn test_case_sensitivity() {
    let nw = setup_nw();
    let test_cases = vec![("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")];

    for (query, subject) in test_cases {
        let mut nw_alignment_mixed = AlignmentData::new(query, subject);

        let aligned_mixed = nw.align(&mut nw_alignment_mixed, false);

        let mut nw_alignment_upper = AlignmentData::new(
            query.to_uppercase().as_str(),
            subject.to_uppercase().as_str(),
        );
        nw.build(&mut nw_alignment_upper);
        let aligned_upper = nw.align(&mut nw_alignment_upper, false);

        assert_eq!(aligned_mixed[0], aligned_upper[0]);

        let sim_mixed = nw.similarity(&mut nw_alignment_mixed);
        let sim_upper = nw.similarity(&mut nw_alignment_upper);
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
    let custom_nw = NeedlemanWunsch::new(&custom_scores);

    let mut nw_alignment = AlignmentData::new("ACGT", "AGT");
    custom_nw.build(&mut nw_alignment);
    assert_eq!(custom_nw.align(&mut nw_alignment, false)[0], "ACGT\nA-GT");

    let query = "AC";
    let subject = "AT";
    let mut nw_alignment_matrix = AlignmentData::new(query, subject);
    custom_nw.build(&mut nw_alignment_matrix);

    let expected_score: Arr2D<i32> = Arr2D::from(&[[0, -3, -6], [-3, 1, -2], [-6, -2, -1]]);

    let score_matrix = &nw_alignment_matrix.score_matrix[0];
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
    let nw = setup_nw();

    let (query, subject) = ("ACCG", "ACG");
    let mut nw_alignment = AlignmentData::new(query, subject);
    let all_aligned = nw.align(&mut nw_alignment, true);

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ACCG\nAC-G".to_string()));
    assert!(all_aligned.contains(&"ACCG\nA-CG".to_string()));

    let (query, subject) = ("ATGTGTA", "ATA");
    let mut nw_alignment = AlignmentData::new(query, subject);
    nw.build(&mut nw_alignment);
    let all_aligned = nw.align(&mut nw_alignment, true);

    assert_eq!(all_aligned.len(), 3);
    assert!(all_aligned.contains(&"ATGTGTA\nAT----A".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA----TA".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA--T--A".to_string()));
}
