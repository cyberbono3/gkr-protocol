use core::str::Chars;
use std::cmp::max;
use std::ops::{Add, Mul};

use ark_ff::Zero;

use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::DenseMVPolynomial;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;

// TODO rename MLPoly to MLPoly
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MLPoly(pub MultiPoly);

impl MLPoly {
    pub fn new(multi_poly: MultiPoly) -> Self {
        Self(multi_poly)
    }

    pub fn neg_shift_poly_by_k(self, k: usize) -> Self {
        let terms = &self.0.terms;
        let current_num_vars = self.0.num_vars;
        let mut shifted_terms = Vec::with_capacity(terms.len());
        for (unit, term) in terms {
            let shifted_term = SparseTerm::new((*term).iter().map(|&c| (c.0 - k, c.1)).collect());
            shifted_terms.push((*unit, shifted_term));
        }
        MultiPoly::from_coefficients_vec(current_num_vars - k, shifted_terms).into()
    }

    pub fn evaluate_variable(&self, r: &Vec<ScalarField>) -> Self {
        if self.0.is_zero() {
            return Self(self.0.clone());
        }
        let mut new_coefficients = Vec::with_capacity(self.0.terms.len());
        let new_num_vars = self.0.num_vars;
        for (unit, term) in &self.0.terms {
            let mut new_unit = *unit;
            let mut new_term = vec![];
            for (var, power) in (*term).iter() {
                if var < &r.len() {
                    for _ in 0..*power {
                        new_unit = new_unit * r[*var];
                    }
                } else {
                    new_term.push((*var, *power));
                }
            }
            new_coefficients.push((new_unit, SparseTerm::new(new_term)));
        }
        Self::new(MultiPoly::from_coefficients_vec(
            new_num_vars,
            new_coefficients,
        ))
    }

    pub fn sum_last_k_var(self, k: usize) -> Self {
        if self.0.is_zero() {
            Self::new(self.0.clone());
        }
        let terms = &self.0.terms;
        let mut new_coefficients = Vec::with_capacity(terms.len());
        let new_num_vars = self.0.num_vars() - k;
        for (unit, term) in terms {
            let mut new_term = Vec::with_capacity(term.len());
            let mut num_reduced_terms = 0;
            for (var, power) in (*term).iter() {
                if (self.0.num_vars() - var) <= k {
                    num_reduced_terms += 1;
                } else {
                    new_term.push((*var, *power));
                }
            }
            let mut new_unit = *unit;
            for _ in 0..2_i32.pow((k - num_reduced_terms) as u32) {
                new_unit += unit;
            }
            new_coefficients.push((new_unit, SparseTerm::new(new_term)));
        }
        Self::new(MultiPoly::from_coefficients_vec(
            new_num_vars,
            new_coefficients,
        ))
    }

    pub fn shift_poly_by_k(&self, k: usize) -> Self {
        let terms = &self.0.terms;
        let current_num_vars = self.0.num_vars;
        let mut shifted_terms = vec![];
        for (unit, term) in terms {
            let shifted_term = SparseTerm::new((*term).iter().map(|c| (c.0 + k, c.1)).collect());
            shifted_terms.push((*unit, shifted_term));
        }
        Self::new(MultiPoly::from_coefficients_vec(
            current_num_vars + k,
            shifted_terms,
        ))
    }

    pub fn restrict_poly_to_line(&self, line: &[UniPoly]) -> UniPoly {
        let mut restricted_poly = UniPoly::zero();
        for (unit, term) in &self.0.terms {
            let variables: Vec<_> = (*term).to_vec();
            let mut term_poly = UniPoly::from_coefficients_slice(&[(0, *unit)]);
            for (var, power) in variables {
                let mut var_poly = line[var].clone();
                for _ in 0..(power - 1) {
                    var_poly = var_poly.mul(&var_poly)
                }
                term_poly = term_poly.mul(&var_poly);
            }
            restricted_poly = restricted_poly + term_poly;
        }
        restricted_poly
    }
}

impl From<MultiPoly> for MLPoly {
    fn from(multi: MultiPoly) -> Self {
        Self(multi)
    }
}

impl Mul for MLPoly {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        let p1_terms = self.0.terms;
        let p2_terms = other.0.terms;
        let num_vars = max(self.0.num_vars, other.0.num_vars);
        let mut mult_terms = Vec::with_capacity(p1_terms.len() * p2_terms.len());
        for (unit_1, term_1) in &p1_terms {
            for (unit_2, term_2) in &p2_terms {
                let mut mult_term: Vec<_> = (*term_1).to_vec();
                mult_term.append(&mut term_2.to_vec());
                mult_terms.push((unit_1 * unit_2, SparseTerm::new(mult_term)));
            }
        }
        MLPoly::new(MultiPoly::from_coefficients_vec(num_vars, mult_terms))
    }
}

impl Add for MLPoly {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        (self.0 + other.0).into()
    }
}
pub struct Binary<'a> {
    inputs: Vec<Chars<'a>>,
    evals: Vec<ScalarField>,
}

impl<'a> Binary<'a> {
    pub fn new(inputs: Vec<Chars<'a>>, evals: Vec<ScalarField>) -> Self {
        //let inputs: Vec<Chars<'a>> = input_strings.into_iter().map(|s| s.chars()).collect();
        Self { inputs, evals }
    }
}

impl<'a> From<Binary<'a>> for MLPoly {
    fn from(binary: Binary<'a>) -> Self {
        let mut terms: Vec<(ScalarField, SparseTerm)> = vec![];
        let num_vars = binary
            .inputs
            .iter()
            .map(|c| c.clone().count())
            .max()
            .unwrap();
        // let mut offset = 0;
        for (input, unit) in binary.inputs.iter().zip(binary.evals) {
            let mut current_term: Vec<(ScalarField, SparseTerm)> =
                Vec::with_capacity(input.clone().count());
            for (idx, char) in input.clone().enumerate() {
                // x_i
                if char == '1' {
                    if current_term.len() == 0 {
                        current_term.append(&mut vec![(unit, SparseTerm::new(vec![(idx, 1)]))])
                    } else {
                        for term in &mut current_term {
                            let mut coeffs = (*term.1.clone()).to_vec();
                            coeffs.push((idx, 1));
                            term.1 = SparseTerm::new(coeffs);
                        }
                    }
                }
                // 1 - x_i
                else if char == '0' {
                    if current_term.len() == 0 {
                        current_term.append(&mut vec![
                            (unit, SparseTerm::new(vec![])),
                            (-unit, SparseTerm::new(vec![(idx, 1)])),
                        ])
                    } else {
                        //  we check the original terms but push a new set of terms multiplied by -x_i
                        let mut new_terms = vec![];
                        for term in &current_term {
                            let mut coeffs = (*term.1.clone()).to_vec();
                            coeffs.push((idx, 1));
                            new_terms.push((-term.0, SparseTerm::new(coeffs)));
                        }
                        current_term.append(&mut new_terms);
                    }
                }
            }
            terms.append(&mut current_term)
        }

        Self::new(MultiPoly::from_coefficients_vec(num_vars, terms))
    }
}

pub struct Input {
    inputs: Vec<usize>,
    evals: Vec<ScalarField>,
    k: usize,
}

impl From<Input> for MLPoly {
    fn from(input: Input) -> Self {
        let k = input.k;
        let str_vec: Vec<String> = input
            .inputs
            .iter()
            .map(|curr| format!("{:0k$b}", curr, k = k))
            .collect();
        let chars_vec: Vec<Chars> = str_vec.iter().map(|s| s.chars()).collect();
        let binary = Binary::new(chars_vec, input.evals);
        binary.into()
    }
}
// pub fn multilinear_polynomial_from_evals(
//     inputs: Vec<usize>,
//     evals: Vec<ScalarField>,
//     k: usize,
// ) -> MultiPoly {
//     let mut binary_inputs = vec![];
//     for curr in inputs {
//         // index of current node in layer as a binary string
//         let curr_string = format!("{:0k$b}", curr, k = k);
//         binary_inputs.push(curr_string);
//     }
//     let input: Vec<Chars> = binary_inputs.iter().map(|s| s.chars()).collect();
//     polynomial_from_binary(input, evals)
// }

pub type UniPoly = UniSparsePolynomial<ScalarField>;
pub struct UVPoly(pub UniPoly);

// TODO to make a from trait
// Converts i into an index in {0,1}^k
pub fn n_to_vec(i: usize, k: usize) -> Vec<ScalarField> {
    format!("{:0k$b}", i, k = k)
        .chars()
        .map(|x| if x == '1' { 1.into() } else { 0.into() })
        .collect()
}

// TODO add Poly wrapper and update API
pub fn mult_poly(p1: &MultiPoly, p2: &MultiPoly) -> MultiPoly {
    let p1_terms = p1.terms();
    let p2_terms = p2.terms();
    let num_vars = max(p1.num_vars(), p2.num_vars());
    let mut mult_terms = vec![];
    for (unit_1, term_1) in p1_terms {
        for (unit_2, term_2) in p2_terms {
            let mut mult_term: Vec<_> = (*term_1).to_vec();
            mult_term.append(&mut term_2.to_vec());
            mult_terms.push((unit_1 * unit_2, SparseTerm::new(mult_term)));
        }
    }
    MultiPoly::from_coefficients_vec(num_vars, mult_terms)
}

pub fn neg_shift_poly_by_k(p: &MultiPoly, k: usize) -> MultiPoly {
    let terms = p.terms();
    let current_num_vars = p.num_vars();
    let mut shifted_terms = vec![];
    for (unit, term) in terms {
        let shifted_term = SparseTerm::new((*term).iter().map(|c| (c.0 - k, c.1)).collect());
        shifted_terms.push((*unit, shifted_term));
    }
    MultiPoly::from_coefficients_vec(current_num_vars - k, shifted_terms)
}

pub fn evaluate_variable(p: &MultiPoly, r: &Vec<ScalarField>) -> MultiPoly {
    if p.is_zero() {
        return p.clone();
    }
    let mut new_coefficients = vec![];
    let new_num_vars = p.num_vars();
    for (unit, term) in p.terms() {
        let mut new_unit = *unit;
        let mut new_term = vec![];
        for (var, power) in (*term).iter() {
            if var < &r.len() {
                for _ in 0..*power {
                    new_unit = new_unit * r[*var];
                }
            } else {
                new_term.push((*var, *power));
            }
        }
        new_coefficients.push((new_unit, SparseTerm::new(new_term)));
    }
    MultiPoly::from_coefficients_vec(new_num_vars, new_coefficients)
}

pub fn sum_last_k_var(p: &MultiPoly, k: usize) -> MultiPoly {
    if p.is_zero() {
        return p.clone();
    }
    let mut new_coefficients = vec![];
    let new_num_vars = p.num_vars() - k;
    for (unit, term) in p.terms() {
        let mut new_term = vec![];
        let mut num_reduced_terms = 0;
        for (var, power) in (*term).iter() {
            if (p.num_vars() - var) <= k {
                num_reduced_terms += 1;
            } else {
                new_term.push((*var, *power));
            }
        }
        let mut new_unit = *unit;
        for _ in 0..2_i32.pow((k - num_reduced_terms) as u32) {
            new_unit += unit;
        }
        new_coefficients.push((new_unit, SparseTerm::new(new_term)));
    }
    MultiPoly::from_coefficients_vec(new_num_vars, new_coefficients)
}

// TODO define from trait
pub fn polynomial_from_binary(inputs: Vec<Chars>, evals: Vec<ScalarField>) -> MultiPoly {
    let mut terms: Vec<(ScalarField, SparseTerm)> = vec![];
    let num_vars = inputs.iter().map(|c| c.clone().count()).max().unwrap();
    // let mut offset = 0;
    for (input, unit) in inputs.iter().zip(evals) {
        let mut current_term: Vec<(ScalarField, SparseTerm)> = vec![];
        for (idx, char) in input.clone().enumerate() {
            // x_i
            if char == '1' {
                if current_term.len() == 0 {
                    current_term.append(&mut vec![(unit, SparseTerm::new(vec![(idx, 1)]))])
                } else {
                    for term in &mut current_term {
                        let mut coeffs = (*term.1.clone()).to_vec();
                        coeffs.push((idx, 1));
                        term.1 = SparseTerm::new(coeffs);
                    }
                }
            }
            // 1 - x_i
            else if char == '0' {
                if current_term.len() == 0 {
                    current_term.append(&mut vec![
                        (unit, SparseTerm::new(vec![])),
                        (-unit, SparseTerm::new(vec![(idx, 1)])),
                    ])
                } else {
                    //  we check the original terms but push a new set of terms multiplied by -x_i
                    let mut new_terms = vec![];
                    for term in &current_term {
                        let mut coeffs = (*term.1.clone()).to_vec();
                        coeffs.push((idx, 1));
                        new_terms.push((-term.0, SparseTerm::new(coeffs)));
                    }
                    current_term.append(&mut new_terms);
                }
            }
        }
        terms.append(&mut current_term)
    }

    MultiPoly::from_coefficients_vec(num_vars, terms)
}

pub fn shift_poly_by_k(p: &MultiPoly, k: usize) -> MultiPoly {
    let terms = p.terms();
    let current_num_vars = p.num_vars();
    let mut shifted_terms = vec![];
    for (unit, term) in terms {
        let shifted_term = SparseTerm::new((*term).iter().map(|c| (c.0 + k, c.1)).collect());
        shifted_terms.push((*unit, shifted_term));
    }
    MultiPoly::from_coefficients_vec(current_num_vars + k, shifted_terms)
}

pub fn multilinear_polynomial_from_evals(
    inputs: Vec<usize>,
    evals: Vec<ScalarField>,
    k: usize,
) -> MultiPoly {
    let mut binary_inputs = vec![];
    for curr in inputs {
        // index of current node in layer as a binary string
        let curr_string = format!("{:0k$b}", curr, k = k);
        binary_inputs.push(curr_string);
    }
    let input: Vec<Chars> = binary_inputs.iter().map(|s| s.chars()).collect();
    polynomial_from_binary(input, evals)
}

pub fn restrict_poly_to_line(p: MultiPoly, line: &[UniPoly]) -> UniPoly {
    let mut restricted_poly = UniPoly::zero();
    for (unit, term) in p.terms() {
        let variables: Vec<_> = (*term).to_vec();
        let mut term_poly = UniPoly::from_coefficients_slice(&[(0, *unit)]);
        for (var, power) in variables {
            let mut var_poly = line[var].clone();
            for _ in 0..(power - 1) {
                var_poly = var_poly.mul(&var_poly)
            }
            term_poly = term_poly.mul(&var_poly);
        }
        restricted_poly = restricted_poly + term_poly;
    }
    restricted_poly
}

pub fn unique_univariate_line(b: &[ScalarField], c: &[ScalarField]) -> Vec<UniPoly> {
    let mut lines = vec![];
    for (b_i, c_i) in b.iter().zip(c) {
        lines.push(UniPoly::from_coefficients_slice(&[
            (0, *b_i),
            (1, c_i - b_i),
        ]));
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr as ScalarField;
    use ark_ff::One;
    use ark_poly::polynomial::multivariate::SparseTerm;

    #[test]
    fn test_mul_basic_case() {
        // Create two multivariate polynomials
        let poly1 = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        let poly2 = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        let result = poly1 * poly2;

        // Expected result should be polynomial multiplication output
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 2)])), // x^2 term
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])), // xy term
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])), // xy term
                (ScalarField::one(), SparseTerm::new(vec![(1, 2)])), // y^2 term
            ],
        ));

        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_zero() {
        // Multiply a polynomial by zero
        let poly1 = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![(ScalarField::one(), SparseTerm::new(vec![(0, 1)]))],
        ));

        let zero_poly = MLPoly::new(MultiPoly::from_coefficients_vec(2, vec![]));

        let result = poly1 * zero_poly;
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(2, vec![]));

        assert_eq!(result, expected);
    }

    #[test]
    fn test_neg_shift_poly_by_k_basic() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=1
        let shifted_poly = poly.neg_shift_poly_by_k(1);

        // Expected result after shift by k=1: just x^1 (the y term shifts out)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            1,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term after shift
            ],
        ));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_neg_shift_poly_by_k_large_shift() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=3 (too large for the number of variables)
        let shifted_poly = poly.neg_shift_poly_by_k(3);

        // Expected result should be an empty polynomial (both variables shifted out)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(0, vec![]));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_evaluate_variable_basic() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Evaluate x=2, y=3
        let evaluation =
            poly.evaluate_variable(&vec![ScalarField::from(2u32), ScalarField::from(3u32)]);

        // Expected result: 2x + 3y evaluated at x=2, y=3 -> 2 + 3 = 5
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            0,
            vec![
                (ScalarField::from(5u32), SparseTerm::new(vec![])), // Constant term 5
            ],
        ));

        assert_eq!(evaluation, expected);
    }

    #[test]
    fn test_evaluate_variable_partial() {
        // Polynomial: x^1 + z^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
            ],
        ));

        // Evaluate x=2, leaving z unevaluated
        let evaluation = poly.evaluate_variable(&vec![ScalarField::from(2u32)]);

        // Expected result: x = 2, z term unevaluated
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::from(2u32), SparseTerm::new(vec![])), // Constant term 2 from x
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term unchanged
            ],
        ));

        assert_eq!(evaluation, expected);
    }

    #[test]
    fn test_sum_last_k_var_basic() {
        // Polynomial: x^1 + y^1 + z^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
            ],
        ));

        // Sum last 1 variable (z)
        let summed_poly = poly.sum_last_k_var(1);

        // Expected result: x^1 + y^1 with z summed out
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (
                    ScalarField::one() + ScalarField::one(),
                    SparseTerm::new(vec![(0, 1)]),
                ), // x term with summed z
                (
                    ScalarField::one() + ScalarField::one(),
                    SparseTerm::new(vec![(1, 1)]),
                ), // y term with summed z
            ],
        ));

        assert_eq!(summed_poly, expected);
    }

    #[test]
    fn test_sum_last_k_var_no_var() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Expected result: same as original
        let expected = poly.clone();

        // Sum last 0 variable (no summing)
        let summed_poly = poly.sum_last_k_var(0);

        assert_eq!(summed_poly, expected);
    }

    #[test]
    fn test_sum_last_k_var_large_k() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Sum last 3 variables (larger than number of variables)
        let summed_poly = poly.sum_last_k_var(3);

        // Expected result: empty polynomial (all variables summed out)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(0, vec![]));

        assert_eq!(summed_poly, expected);
    }

    #[test]
    fn test_sum_last_k_var_partial() {
        // Polynomial: x^1 + y^1 + z^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
            ],
        ));

        // Sum last 2 variables (y and z)
        let summed_poly = poly.sum_last_k_var(2);

        // Expected result: just x^1, y and z summed out
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            1,
            vec![
                (ScalarField::from(2u32), SparseTerm::new(vec![(0, 1)])), // x term with y and z summed out
            ],
        ));

        assert_eq!(summed_poly, expected);
    }

    #[test]
    fn test_binary_to_MLPoly_single_input() {
        // Binary: input "10", eval 1
        let inputs = vec!["10".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (x_0 * (1 - x_1))
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MLPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MLPoly_multiple_inputs() {
        // Binary: inputs "10", "01", evals 1, -1
        let inputs = vec!["10".chars(), "01".chars()];
        let evals = vec![ScalarField::one(), -ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial:
        // 1 * (x_0 * (1 - x_1)) + (-1) * ((1 - x_0) * x_1)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(1, 1)])),
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MLPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MLPoly_all_zeros() {
        // Binary: input "00", eval 1
        let inputs = vec!["00".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (1 - x_0) * (1 - x_1)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(1, 1)])),
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MLPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MLPoly_all_ones() {
        // Binary: input "11", eval 1
        let inputs = vec!["11".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (x_0 * x_1)
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![(ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)]))],
        ));

        let ml_poly = MLPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_shift_poly_by_k_basic() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=1
        let shifted_poly = poly.shift_poly_by_k(1);

        // Expected result: x^1 -> z^1, y^1 -> (z+1)^1
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // shifted x term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // shifted y term
            ],
        ));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_shift_poly_by_k_large_shift() {
        // Polynomial: x^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            1,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
            ],
        ));

        // Shift by k=3
        let shifted_poly = poly.shift_poly_by_k(3);

        // Expected result: x^1 -> (z+3)^1
        let expected = MLPoly::new(MultiPoly::from_coefficients_vec(
            4,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(3, 1)])), // shifted x term
            ],
        ));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_shift_poly_by_k_zero_shift() {
        // Polynomial: x^1 + y^1
        let poly = MLPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=0 (no shift)
        let shifted_poly = poly.shift_poly_by_k(0);

        // Expected result: same as original
        let expected = poly.clone();

        assert_eq!(shifted_poly, expected);
    }

    // TODO Add tests restrict_poly_to_line
}
