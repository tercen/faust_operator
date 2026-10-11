# syntax=docker/dockerfile:1.7
# faust_operator — static tier (create-rust-operator §5): one musl binary on scratch.

# ---- builder ----
FROM rust:1.94-bookworm AS builder
RUN apt-get update && apt-get install -y --no-install-recommends \
        musl-tools protobuf-compiler pkg-config git ca-certificates \
 && rm -rf /var/lib/apt/lists/* && rustup target add x86_64-unknown-linux-musl
WORKDIR /build
RUN cargo install cargo-chef --locked

# Dependency layer. cargo-chef's recipe carries the package version, and a release commit always
# bumps it, so normalise the version for the recipe only, then restore the real manifests.
COPY Cargo.toml Cargo.lock ./
COPY crates/faust/Cargo.toml crates/faust/Cargo.toml
RUN cp Cargo.toml Cargo.toml.keep \
 && sed -i 's/^version = ".*"/version = "0.0.0"/' Cargo.toml \
 && mkdir -p crates/faust/src && echo '' > crates/faust/src/lib.rs \
 && cargo chef prepare --recipe-path recipe.json
RUN cargo chef cook --release --target x86_64-unknown-linux-musl --recipe-path recipe.json
RUN mv Cargo.toml.keep Cargo.toml

COPY src ./src
COPY crates ./crates
COPY operator.json ./
RUN cargo build --release --target x86_64-unknown-linux-musl --bin faust_operator \
 && ls -l target/x86_64-unknown-linux-musl/release/faust_operator
RUN mkdir -p /tmp-op && chown 1000:1000 /tmp-op

# ---- runtime ----
FROM scratch
COPY --from=builder /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/
COPY --from=builder /build/target/x86_64-unknown-linux-musl/release/faust_operator /usr/local/bin/faust_operator
COPY --from=builder --chown=1000:1000 /tmp-op /tmp
USER 1000:1000
WORKDIR /operator
ENV RUST_BACKTRACE=1 RUST_LOG=info TMPDIR=/tmp
ENTRYPOINT ["/usr/local/bin/faust_operator"]
