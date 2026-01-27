# goombay-rs

Rust implementation of my goombay python project. Implements several sequence
alignment algorithms such as Waterman-Smith-Beyer, Gotoh, and Needleman-Wunsch
intended to calculate distance, show alignment, and display the underlying matrices.

## Installation

Add Goombay-rs as a dependency in your Cargo.toml:

```nginx
[dependencies]
goombay-rs = { git = "https://github.com/lignum-vitae/goombay-rs" }`
```

Then run:

`cargo build`

## Project layout

| Module  | Description                                                                       |
| ------- | --------------------------------------------------------------------------------- |
| align   | Alignment algorithms such as Needleman-Wunsch, Wagner-Fischer, and Smith-Waterman |
| scoring | Scoring structs for custom scoring in alignment algorithms                        |

### Running Examples

Working examples of the available algorithms as well as a full list of available algorithms can be found in the
[`examples/`](https://github.com/lignum-vitae/goombay-rs/tree/main/goombay-rs/examples) directory.

Run any example with the following command:

`cargo run --example <example_name>`

Do not include `.rs` when running examples.

## Contributing

We welcome contributions! Please read our:

- [Code of Conduct](https://github.com/lignum-vitae/goombay-rs/blob/main/docs/CODE_OF_CONDUCT.md)
- [Contribution Guidelines](https://github.com/lignum-vitae/goombay-rs/blob/main/docs/CONTRIBUTING.md)

> [!NOTE]
> Before submitting a PR, install [just](https://github.com/casey/just) and run `just check`
> to pull the latest changes from the main branch as well as to format, test, and lint your code.
> Just can be installed using `cargo install just`, curl, or your favourite package manager.

## Stability

This project is in the alpha stage. APIs may change without warning until version
1.0.0.

## License

This project is licensed under the MIT License - see the
[LICENSE](https://github.com/lignum-vitae/goombay-rs/blob/main/LICENSE) file for details.
