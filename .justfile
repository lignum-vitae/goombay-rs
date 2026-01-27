set shell := ["sh", "-c"]
set windows-shell := ["cmd.exe", "/c"]

ROOT_DIR := justfile_directory()
MANIFEST_PATH_FLAG := f"--manifest-path {{ROOT_DIR}}/goombay-rs/Cargo.toml"

recipes:
    just -l

test:
    cargo test {{MANIFEST_PATH_FLAG}}

lint:
    cargo clippy --all-targets --all-features {{MANIFEST_PATH_FLAG}}

format:
    cargo fmt {{MANIFEST_PATH_FLAG}}

latest:
    git pull origin main

check: latest format test lint
