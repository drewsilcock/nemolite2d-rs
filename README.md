[![Build Status](https://github.com/drewsilcock/nemolite2d-rs/workflows/NEMOLite2D%20Rust/badge.svg)](https://github.com/drewsilcock/nemolite2d-rs/actions)

# NEMOLite2D Rust

This is a port of the NEMOLite2D codebase in Rust, for benchmarking both Rust and WebAssembly for scientific computing.

For more information about NEMOLite2D, see:
  - [The STFC site for NEMOLite2D](https://edata.stfc.ac.uk/handle/edata/729)
  - [The NERC site for GOcean](https://puma.nerc.ac.uk/trac/GOcean), which is part of the same family of projects and benchmarking suites.

The original code can be found at either of the two above web pages, but the [GitHub repository for PSyclone benchmarking](https://github.com/stfc/PSycloneBench) is the easiest way to quickly get the most up to date version of the original code. Note that this contains multiple versions of multiple projects. The original serial version which was used for this port can be found at [PSycloneBench/ocean/nemo/nemolite2d/original/nemolite2d.f90](https://github.com/stfc/PSycloneBench/blob/master/ocean/nemo/nemolite2d/original/nemolite2d.f90).

## Running natively

To run natively, install `rustup` from https://rustup.rs/ if you haven't already and run:

```sh
$ cargo run --release
```

## Running on Wasm

### WASI - WebAssembly System Interface

Running with a WASI-compatible runtime means that you can simply build the `.wasm` file for the project and run it exactly like the native executable using one of several available Wasm runtimes.

To build the `.wasm` output, you first need to install the WASI libraries for Rust:

```sh
$ rustup target add wasm32-wasi
$ cargo build --target=wasm32-wasi --release
```

The major WASI-compatible runtimes are wasmtime and Wasmer. To run with wasmtime:

```sh
$ cargo install wasmtime
$ wasmtime target/wasm32-wasi/release/nemolite2d-rs.wasm --dir=.
```

To run with Wasmer:

```sh
$ curl https://get.wasmer.io -sSfL | sh
$ wasmer run target/wasm32-wasi/release/nemolite2d-rs.wasm --dir .
```

### Running in a web browser via JavaScript

TODO: Add necessary wasm-bindgen attributes and create HTML/JS that runs simulation and shows performance report.
