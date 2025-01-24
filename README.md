<h1 align="center">GKR protocol</h1>

`gkr-protocol` is a Rust library that implements the GKR protocol

The GKR is an interactive proof protocol widely used in ZK/Validity proofs.

> [!Warning]
> This is a side project and hasn't been audited.
> Use at your own risk.

## Overfiew

1. Prover P and Verifier V both agree on input `x` (element of finite field F) and some fan-2 arithmetic circuit `C` over `F`, where `x` is an input, and `y` is an output. 
2. GKR is interactive proof protocol for claim `C(x) = y`
3. For instance, circuit input size is 4, circuit size is 11 gates.
4. If circuit size is much bigger than input, V saves a bunch of time on it.
5. GKR protocol has a layered structure : 
    - input layer
    - intermediate layers
    - output layer
6. For every layer gate-value function `W` is computed
7. For every layer MLE(multilinear extension) of `W` id computed 

## Protocol steps

1. P sends a claimed output `y`. Let `y_tilda` be an MLE of `y`.
2. V picks tandom `r` from field `F`
3. By Schwarz-Zipper lemma it is "safe" for V to belive `C(x) = y`
4. V checks if `s_tilda(r) = y_tilda(r)` holds, where `s_tilda` MLE computed by V, `y_tilda` is claimed output(MLE) by P. If it does hold, `a_tilda` and `y_tilda` are equal. Therefore, polynomials are equal with the high probability, according to Schwarz-Zippel lemma.
5. Apply sumcheck protocol (described below) to compute the sum for every layer
6. For instance for protocol with 3 layers, it looks like that, 
    - `w_0_tilda(r)` (output layet ) -> `w_1_tilda(r0)` (intemediate layer) -> `w_2_tilda(r1)` (input layer )
    - `r`,`r0`,`r1` is random field elements for every layer
    - `w0_tilda` is lagrange interpolation of gate-value function for layer 0 
    - `w1_tilda` is lagrange interpolation of gate-value function for layer 1
    - `w2_tilda` is lagrange interpolation of gate-value function for layer 2

## Important: 
The whole point of GKR  protocol is to avoid committing to intermediate values of a circuit. 
Assumption: Without commitments, V knows an input and goes layer by layer.

## Sumcheck protocol overview

Let `g` be n-variate polynomial over finite field `F`. Let `g` have degree 3.

Compute sum of `g(x)` over input `x = {0,1}^n`

Task: offload hard work of computing a sum to prover P.

This is public coin procool, so we can apply Fiat-Shamir to make it non-interactive.

Procotol has of `n` rounds ( number of variables in polynomial `g`).

Verifier `V` time O(n) field ops

## Sumcheck protocol steps:

1. `P` sends claimed answer C,  `S1(x)` claimed to be equal `h(x)  = sum of g(x) over input (x1...xn)`
2. V picks random `r` from finite field `F` and sends it to P
3. V checks if `S(r) = h(r)` holds. Completeness: if prover P honest, this check will pass.
4. Repeat this process `n` rounds.
    Soundness error <= n/|F|. As long as field F is big enouph, we keep this probability negligible.

**WARNING**: This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

**WARNING**: This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo` (the standard Rust build tool) to build the library:
```bash
git clone https://github.com/cyberbono3/gkr-protocol.git
cd gkr-protocol
cargo build --release
```

This library comes with some unit and integration tests. Run these tests with:
```bash
cargo test
```

Lastly, this library is instrumented with profiling infrastructure that prints detailed traces of execution time. To enable this, compile with `cargo build --features print-trace`.

## License

This library is licensed under either of the following licenses, at your discretion.

* [Apache License Version 2.0](LICENSE-APACHE)
* [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## Reference Paper
[Proofs, Arguments, and Zero-Knowledge](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf) <br/>
Justin Thaler



