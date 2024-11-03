use super::fiat_shamir::FiatShamir;
use crate::poly::HyperCube;

use ark_bn254::Fr as ScalarField;
use ark_ff::Field;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::{DenseMVPolynomial, Polynomial};
use rand::Rng;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<ScalarField>;

// Simulates memory of a single prover instance
#[derive(Debug, Clone)]
pub struct Prover {
    pub g: MultiPoly,
    pub fs: FiatShamir,
    pub r_vec: Vec<ScalarField>,
}

impl Prover {
    const N_ROUNDS: usize = 2;
    pub fn new(g: &MultiPoly) -> Self {
        Prover {
            g: g.clone(),
            fs: FiatShamir::new(vec![], Self::N_ROUNDS),
            r_vec: vec![],
        }
    }

    pub fn gen_uni_polynomial(&mut self, r: Option<ScalarField>) -> UniPoly {
        if let Some(r) = r {
            self.r_vec.push(r);
        }

        let k = self.g.num_vars() - self.r_vec.len();
        let mut result = UniPoly::from_coefficients_vec(vec![(0, 0u32.into())]);

        for n in 0..(2u32.pow(k as u32 - 1)) {
            let term = self.evaluate_gj(HyperCube::new(n as usize, k).0.clone());
            result = result + term;
        }

        result
    }

    // Evaluates gj over a vector permutation of points, folding all evaluated terms together into one univariate polynomial
    pub fn evaluate_gj(&self, points: Vec<ScalarField>) -> UniPoly {
        self.g.terms().iter().fold(
            UniPoly::from_coefficients_vec(vec![]),
            |sum, (coeff, term)| {
                let (coeff_eval, fixed_term) = self.evaluate_term(term, &points);
                let curr = if let Some(term) = fixed_term {
                    UniPoly::from_coefficients_vec(vec![(term.degree(), *coeff * coeff_eval)])
                } else {
                    UniPoly::from_coefficients_vec(vec![(0, *coeff * coeff_eval)])
                };
                curr + sum
            },
        )
    }

    // Evaluates a term with a fixed univar, returning (new coefficent, fixed term)
    pub fn evaluate_term(
        &self,
        term: &SparseTerm,
        point: &[ScalarField],
    ) -> (ScalarField, Option<SparseTerm>) {
        let mut fixed_term: Option<SparseTerm> = None;
        let coeff: ScalarField =
            term.iter()
                .fold(1u32.into(), |product, (var, power)| match *var {
                    j if j == self.r_vec.len() => {
                        fixed_term = Some(SparseTerm::new(vec![(j, *power)]));
                        product
                    }
                    j if j < self.r_vec.len() => self.r_vec[j].pow([*power as u64]) * product,
                    _ => point[*var - self.r_vec.len()].pow([*power as u64]) * product,
                });
        (coeff, fixed_term)
    }

    // Verify prover's claim c_1
    pub fn verify(&mut self, c_1: ScalarField) -> bool {
        // 1st round
        let mut gi = self.gen_uni_polynomial(None);
        let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(c_1, expected_c);
        let lookup_degree = max_degrees(&self.g);
        assert!(gi.degree() <= lookup_degree[0]);

        // middle rounds
        // for j in 1..self.g.num_vars() {
        for degree in lookup_degree.iter().take(self.g.num_vars()).skip(1) {
            let r = self.get_r();
            expected_c = gi.evaluate(&r.unwrap());
            gi = self.gen_uni_polynomial(r);
            let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
            assert_eq!(expected_c, new_c);
            assert!(gi.degree() <= *degree);
        }
        // final round
        let r = self.get_r();
        expected_c = gi.evaluate(&r.unwrap());
        self.r_vec.push(r.unwrap());
        let new_c = self.g.evaluate(&self.r_vec);
        assert_eq!(expected_c, new_c);
        true
    }

    // Verify prover's claim c_1
    pub fn verify_non_interactive(&mut self, c_1: ScalarField) -> bool {
        self.fs.circuit_input = vec![c_1];
        // 1st round
        let mut gi = self.gen_uni_polynomial(None);
        let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(c_1, expected_c);
        let lookup_degree = max_degrees(&self.g);
        assert!(gi.degree() <= lookup_degree[0]);

        // middle rounds
        //for j in 1..self.g.num_vars() {
        for degree in lookup_degree.iter().take(self.g.num_vars()).skip(1) {
            let r = self.fs.get_r(gi.clone());
            expected_c = gi.evaluate(&r);
            gi = self.gen_uni_polynomial(Some(r));
            let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
            assert_eq!(expected_c, new_c);
            assert!(gi.degree() <= *degree);
        }
        // final round
        let r = self.fs.get_r(gi.clone());
        expected_c = gi.evaluate(&r);
        self.r_vec.push(r);
        let new_c = self.g.evaluate(&self.r_vec);
        assert_eq!(expected_c, new_c);
        true
    }

    // Verifier procedures
    pub fn get_r(&self) -> Option<ScalarField> {
        let mut rng = rand::thread_rng();
        let r: ScalarField = rng.gen();
        Some(r)
    }
}

// A degree look up table for all variables in g
pub fn max_degrees(g: &MultiPoly) -> Vec<usize> {
    let mut lookup: Vec<usize> = vec![0; g.num_vars()];
    g.terms().iter().for_each(|(_, term)| {
        term.iter().for_each(|(var, power)| {
            if *power > lookup[*var] {
                lookup[*var] = *power
            }
        });
    });
    lookup
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr as ScalarField;
    use ark_poly::{
        polynomial::multivariate::{SparseTerm, Term},
        DenseMVPolynomial,
    };

    // Helper function to create a sample ScalarField element
    fn sample_scalar(val: u64) -> ScalarField {
        ScalarField::from(val)
    }

    // Helper function to create a multivariate polynomial with a specified number of variables
    fn sample_multivariate_poly(num_vars: usize) -> MultiPoly {
        let terms = (0..num_vars)
            .map(|i| {
                let coeff = sample_scalar((i + 1) as u64); // Coefficients 1, 2, 3, etc.
                let term = SparseTerm::new(vec![(i, 1)]); // Single term per variable
                (coeff, term)
            })
            .collect::<Vec<_>>();

        MultiPoly::from_coefficients_vec(num_vars, terms)
    }

    // TODO revise test
    // #[test]
    // fn test_gen_uni_polynomial_no_r() {
    //     // Create a sample prover with a multivariate polynomial with 2 variables
    //     let g = sample_multivariate_poly(2);
    //     let mut prover = Prover::new(&g);

    //     // Call gen_uni_polynomial without providing r
    //     let uni_poly = prover.gen_uni_polynomial(None);

    //     // Check basic properties of the result (e.g., degrees, terms)
    //     assert!(uni_poly.degree() >= 0); // Ensure polynomial is not empty
    // }

    // #[test]
    // fn test_gen_uni_polynomial_with_r() {
    //     // Create a sample prover with a multivariate polynomial with 3 variables
    //     let g = sample_multivariate_poly(3);
    //     let mut prover = Prover::new(&g);

    //     // Provide a specific r value
    //     let r = sample_scalar(5);
    //     let uni_poly = prover.gen_uni_polynomial(Some(r));

    //     // Check if r was added to r_vec
    //     assert_eq!(prover.r_vec.len(), 1);
    //     assert_eq!(prover.r_vec[0], r);

    //     // Check basic properties of the result (e.g., degrees, terms)
    //     assert!(uni_poly.degree() >= 0); // Ensure polynomial is not empty
    // }

    #[test]
    fn test_gen_uni_polynomial_multiple_r() {
        // Create a sample prover with a multivariate polynomial with 4 variables
        let g = sample_multivariate_poly(4);
        let mut prover = Prover::new(&g);

        // Provide multiple r values
        prover.gen_uni_polynomial(Some(sample_scalar(5)));
        prover.gen_uni_polynomial(Some(sample_scalar(7)));

        // Check if r_vec has both values
        assert_eq!(prover.r_vec.len(), 2);
        assert_eq!(prover.r_vec[0], sample_scalar(5));
        assert_eq!(prover.r_vec[1], sample_scalar(7));
    }

    #[test]
    pub fn test_can_verify() {
        let poly = MultiPoly::from_coefficients_vec(
            1,
            vec![(ScalarField::from(2), SparseTerm::new(vec![(0, 3)]))],
        );

        let sum = ScalarField::from(2);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);

        let poly = MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::from(2), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (ScalarField::from(5), SparseTerm::new(vec![(2, 1)])),
                (ScalarField::from(1), SparseTerm::new(vec![])),
            ],
        );
        let sum = ScalarField::from(32);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);
    }

    #[test]
    pub fn test_can_verify_non_interactive() {
        let poly = MultiPoly::from_coefficients_vec(
            1,
            vec![(ScalarField::from(2), SparseTerm::new(vec![(0, 3)]))],
        );
        let sum = ScalarField::from(2);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);

        let poly = MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::from(2), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (ScalarField::from(5), SparseTerm::new(vec![(2, 1)])),
                (ScalarField::from(1), SparseTerm::new(vec![])),
            ],
        );
        let sum = ScalarField::from(32);

        let mut prover = Prover::new(&poly);

        prover.verify_non_interactive(sum);
    }

    #[test]
    #[should_panic]
    pub fn test_fails() {
        let poly = MultiPoly::from_coefficients_vec(
            1,
            vec![(ScalarField::from(2), SparseTerm::new(vec![(0, 3)]))],
        );
        let sum = ScalarField::from(2);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);

        let poly = MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::from(2), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (ScalarField::from(5), SparseTerm::new(vec![(2, 1)])),
                (ScalarField::from(1), SparseTerm::new(vec![])),
            ],
        );
        let sum = ScalarField::from(36);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);
    }

    #[test]
    #[should_panic]
    pub fn test_fails_non_interactive() {
        let poly = MultiPoly::from_coefficients_vec(
            1,
            vec![(ScalarField::from(2), SparseTerm::new(vec![(0, 3)]))],
        );
        let sum = ScalarField::from(2);

        let mut prover = Prover::new(&poly);

        prover.verify(sum);

        let poly = MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::from(2), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (ScalarField::from(5), SparseTerm::new(vec![(2, 1)])),
                (ScalarField::from(1), SparseTerm::new(vec![])),
            ],
        );
        let sum = ScalarField::from(36);

        let mut prover = Prover::new(&poly);

        prover.verify_non_interactive(sum);
    }
}
