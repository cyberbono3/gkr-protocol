use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

use gkr_protocol::{gkr::prover::Prover, poly::HyperCube, poly::MultiPoly};

use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;

// Example function that creates a "dummy" polynomial for testing
fn create_test_multipoly(num_vars: usize) -> SparsePolynomial<ScalarField, SparseTerm> {
    // A very minimal polynomial for demonstration
    // You can expand coefficients and terms as needed
    SparsePolynomial::from_coefficients_vec(num_vars, vec![])
}

fn bench_prover_creation(c: &mut Criterion) {
    c.bench_function("Prover::new", |b| {
        b.iter(|| {
            let poly = black_box(create_test_multipoly(2));
            let _prover = Prover::new(&poly);
        })
    });
}

fn bench_gen_uni_polynomial(c: &mut Criterion) {
    // Setup: create a Prover with a simple polynomial
    let poly = create_test_multipoly(2);
    let mut prover = Prover::new(&poly);

    c.bench_function("Prover::gen_uni_polynomial", |b| {
        b.iter_batched(
            || None, // input setup: for instance, no random 'r' for this bench
            |r| {
                // The function you want to measure
                black_box(prover.gen_uni_polynomial(r))
            },
            BatchSize::SmallInput,
        )
    });
}

// Register the benchmarks
criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = bench_prover_creation, bench_gen_uni_polynomial
);
criterion_main!(benches);
