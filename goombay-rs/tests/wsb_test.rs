use goombay_rs::align::{GlobalAlignmentMatrix, WatermanSmithBeyer};
use goombay_rs::scoring::ExtendedGapScoring;

#[test]
fn test_identical_sequences() {
    let wsb = WatermanSmithBeyer::compute("ACTG", "ACTG");

    let aligned = wsb.align();
    let sim = wsb.similarity();

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * wsb.identity) as i32);
}

#[test]
fn test_completely_different() {
    let wsb = WatermanSmithBeyer::compute("AAAA", "TTTT");
    let aligned = wsb.align();
    let sim = wsb.similarity();

    assert_eq!(aligned[0], "AAAA\nTTTT");
    assert_eq!(sim, -4_i32 * wsb.mismatch as i32);
}

#[test]
fn test_different_length() {
    // With affine gap penalty, a single contiguous gap should be preferred
    // over scattered gaps
    let wsb = WatermanSmithBeyer::compute("ACGT", "AGT");
    let aligned = wsb.align();
    assert_eq!(aligned[0], "ACGT\nA-GT");
}

#[test]
fn test_empty_sequences() {
    let custom_scores = ExtendedGapScoring {
        identity: 2,
        mismatch: 1,
        gap: 3,
        extended_gap: 1,
    };
    let custom_wsb = WatermanSmithBeyer::set_scores(&custom_scores);

    // Empty query
    let wsb = custom_wsb.calculate_matrix("", "ACTG");
    let aligned = wsb.align();
    assert_eq!(aligned[0], "----\nACTG");

    // Empty subject
    let wsb = custom_wsb.calculate_matrix("ACTG", "");
    let aligned = wsb.align();
    assert_eq!(aligned[0], "ACTG\n----");

    // Both empty
    let wsb = custom_wsb.calculate_matrix("", "");
    let sim = wsb.similarity();
    assert_eq!(sim, 1);
}

#[test]
fn test_single_character() {
    let wsb_match = WatermanSmithBeyer::compute("A", "A");
    assert_eq!(wsb_match.align()[0], "A\nA");
    assert_eq!(wsb_match.similarity(), wsb_match.identity as i32);

    let wsb_mismatch = WatermanSmithBeyer::compute("A", "T");
    assert_eq!(wsb_mismatch.align()[0], "A\nT");
    assert_eq!(wsb_mismatch.similarity(), -(wsb_mismatch.mismatch as i32));
}

#[test]
fn test_case_sensitivity() {
    let test_cases = vec![("ACTG", "actg"), ("AcTg", "aCtG"), ("actg", "ACTG")];

    for (query, subject) in test_cases {
        let wsb_mixed = WatermanSmithBeyer::compute(query, subject);
        let wsb_upper = WatermanSmithBeyer::compute(
            query.to_uppercase().as_str(),
            subject.to_uppercase().as_str(),
        );

        assert_eq!(wsb_mixed.align()[0], wsb_upper.align()[0]);
        assert_eq!(wsb_mixed.similarity(), wsb_upper.similarity());
    }
}

#[test]
fn test_affine_gap_prefers_contiguous_gaps() {
    // With affine gap penalty (gap_open=3, gap_extend=1), a single gap of 2
    // costs 3+1=4, while two gaps of 1 each cost 3+3=6.
    // So the algorithm should prefer contiguous gaps over scattered gaps.
    let scores = ExtendedGapScoring {
        identity: 2,
        mismatch: 5,
        gap: 3,
        extended_gap: 1,
    };
    let wsb = WatermanSmithBeyer::set_scores(&scores);

    // "AACT" vs "AT": should get A--T aligned (one gap of 2) rather than
    // A-C-T or similar with scattered gaps
    let result = wsb.calculate_matrix("AACT", "AT");
    let aligned = result.align();
    assert_eq!(aligned[0], "AACT\nA--T");
}

#[test]
fn test_all_alignments_1() {
    let wsb = WatermanSmithBeyer::compute("ACCGT", "CT");
    let all_aligned = wsb.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ACCGT\n---CT".to_string()));
    assert!(all_aligned.contains(&"ACCGT\nC---T".to_string()));
}

#[test]
fn test_all_alignments_2() {
    let wsb = WatermanSmithBeyer::compute("ATGTGTA", "ATA");
    let all_aligned = wsb.all_alignments(true).align();

    assert_eq!(all_aligned.len(), 2);
    assert!(all_aligned.contains(&"ATGTGTA\nAT----A".to_string()));
    assert!(all_aligned.contains(&"ATGTGTA\nA----TA".to_string()));
}

#[test]
fn test_custom_scoring() {
    let scores = ExtendedGapScoring {
        identity: 1,
        mismatch: 1,
        gap: 2,
        extended_gap: 1,
    };
    let wsb = WatermanSmithBeyer::set_scores(&scores);
    let result = wsb.calculate_matrix("ACGT", "AGT");

    let aligned = result.align();
    assert_eq!(aligned[0], "ACGT\nA-GT");
}

#[test]
fn test_similarity_and_distance() {
    let wsb = WatermanSmithBeyer::compute("ACTG", "ACTG");
    assert_eq!(wsb.similarity(), (4 * wsb.identity) as i32);
    assert_eq!(wsb.distance(), 0);
    assert_eq!(wsb.normalized_similarity(), 1.0);
    assert_eq!(wsb.normalized_distance(), 0.0);
}
