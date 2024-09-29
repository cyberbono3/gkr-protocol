use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;