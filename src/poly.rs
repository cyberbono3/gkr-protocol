use std::cmp::max;

use ark_ff::Zero;

use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::DenseMVPolynomial;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;

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
