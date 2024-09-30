use core::str::Chars;
use std::cmp::max;

use ark_ff::Zero;

use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::DenseMVPolynomial;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<ScalarField>;

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
