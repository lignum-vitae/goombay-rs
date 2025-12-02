# Goombay-rs Examples

## Algorithms

The following algorithms are currently supported:

- Needleman-Wunsch
- Wagner-Fischer

## Scoring

These stucts can be accessed using the import `use goombay-rs::scoring::<name>`.
The purpose of these structs for the end user is to allow for custom scoring.

```rust
#[derive(Clone)]
pub struct LevenshteinScoring {
    pub substitution: usize,
    pub gap: usize,
}

#[derive(Clone)]
pub struct GeneralScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
}

#[derive(Clone)]
pub struct TransposeScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub transpose: usize,
}

#[derive(Clone)]
pub struct ExtendedGapScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub extended_gap: usize,
}
```

