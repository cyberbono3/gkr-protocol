use core::str::Chars;
use std::cmp::max;
use std::ops::{Add, Mul};

use ark_ff::Zero;

use ark_bn254::Fr as ScalarField;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::DenseMVPolynomial;
use thiserror::Error;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;

#[derive(Error, Debug, PartialEq)]
pub enum PolyError {
    #[error("subtract with overflow is not allowed")]
    SubtractWithOverflow,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MVPoly(pub MultiPoly);

pub struct SparseTermWrapper(pub SparseTerm, pub usize);

impl TryFrom<SparseTermWrapper> for Vec<(usize, usize)> {
    type Error = PolyError;
    fn try_from(wrapper: SparseTermWrapper) -> Result<Vec<(usize, usize)>, Self::Error> {
        let k = wrapper.1;
        let mut sparse_term: Vec<(usize, usize)> = Vec::with_capacity(wrapper.0.len());
        for c in wrapper.0.iter() {
            if let Some(subtracted) = c.0.checked_sub(k) {
                sparse_term.push((subtracted, c.1));
            } else {
                return Err(PolyError::SubtractWithOverflow);
            }
        }
        Ok(sparse_term)
    }
}

impl MVPoly {
    pub fn new(multi_poly: MultiPoly) -> Self {
        Self(multi_poly)
    }

    pub fn neg_shift_by_k(self, k: usize) -> Result<Self, PolyError> {
        let terms = &self.0.terms;
        let current_num_vars = self.0.num_vars;
        let mut shifted_terms = Vec::with_capacity(terms.len());
        for (unit, term) in terms {
            let sparse_term: Vec<(usize, usize)> = SparseTermWrapper(term.clone(), k).try_into()?;
            let shifted_term = SparseTerm::new(sparse_term);
            shifted_terms.push((*unit, shifted_term));
        }
        Ok(MultiPoly::from_coefficients_vec(current_num_vars - k, shifted_terms).into())
    }

    pub fn evaluate_variable(self, r: &[ScalarField]) -> Self {
        if self.0.is_zero() {
            return Self(self.0.clone());
        }
        let mut new_coefficients = Vec::with_capacity(self.0.terms.len());
        let new_num_vars = self.0.num_vars;
        for (unit, term) in &self.0.terms {
            let mut new_unit = *unit;
            let mut new_term = Vec::with_capacity(term.len());
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

    // TODO test it
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

    pub fn shift_by_k(&self, k: usize) -> Self {
        let terms = &self.0.terms;
        let current_num_vars = self.0.num_vars;
        let mut shifted_terms = Vec::with_capacity(terms.len());
        for (unit, term) in terms {
            let shifted_term = SparseTerm::new((*term).iter().map(|c| (c.0 + k, c.1)).collect());
            shifted_terms.push((*unit, shifted_term));
        }
        Self::new(MultiPoly::from_coefficients_vec(
            current_num_vars + k,
            shifted_terms,
        ))
    }

    pub fn restrict_to_line(&self, line: &[UniPoly]) -> UniPoly {
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

impl From<MultiPoly> for MVPoly {
    fn from(multi: MultiPoly) -> Self {
        Self(multi)
    }
}

impl Mul for MVPoly {
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
        MVPoly::new(MultiPoly::from_coefficients_vec(num_vars, mult_terms))
    }
}

impl Add for MVPoly {
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

impl<'a> From<Binary<'a>> for MVPoly {
    fn from(binary: Binary<'a>) -> Self {
        let mut terms: Vec<(ScalarField, SparseTerm)> = Vec::with_capacity(binary.inputs.len());
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

pub struct PolyInput {
    inputs: Vec<usize>,
    evals: Vec<ScalarField>,
    k: usize,
}

impl PolyInput {
    pub fn new(inputs: Vec<usize>, evals: Vec<ScalarField>, k: usize) -> Self {
        Self { inputs, evals, k }
    }
}

impl From<PolyInput> for MVPoly {
    fn from(input: PolyInput) -> Self {
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

pub type UniPoly = UniSparsePolynomial<ScalarField>;
//pub struct UVPoly(pub UniPoly);

//{0,1}^k
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HyperCube(pub Vec<ScalarField>);

impl HyperCube {
    pub fn new(n: usize, k: usize) -> Self {
        let vec: Vec<ScalarField> = format!("{:0k$b}", n, k = k)
            .chars()
            .map(|x| if x == '1' { 1.into() } else { 0.into() })
            .collect();
        Self(vec)
    }
}

impl From<Vec<ScalarField>> for HyperCube {
    fn from(vec: Vec<ScalarField>) -> Self {
        Self(vec)
    }
}

pub fn unique_univariate_line(b: &[ScalarField], c: &[ScalarField]) -> Vec<UniPoly> {
    let mut lines = Vec::with_capacity(b.iter().len());
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
    use ark_ff::{BigInteger256, One};
    use ark_poly::polynomial::multivariate::SparseTerm;

    #[test]
    fn test_hyper_cube_vec_basic() {
        // Test with i = 0 and k = 3, expecting [0, 0, 0]
        let result = HyperCube::new(0, 3);
        let scalar_field_vec = vec![0.into(), 0.into(), 0.into()];
        let expected: HyperCube = scalar_field_vec.into();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hyper_cube_single_one() {
        // Test with i = 1 and k = 3, expecting [0, 0, 1]
        let result = HyperCube::new(1, 3);
        let scalar_field_vec = vec![0.into(), 0.into(), 1.into()];
        let expected: HyperCube = scalar_field_vec.into();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hyper_cube_to_vec_mixed_bits() {
        // Test with i = 5 and k = 3, expecting [1, 0, 1]
        let result = HyperCube::new(5, 3);
        let expected: HyperCube = vec![1.into(), 0.into(), 1.into()].into();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hyper_cube_to_vec_full_bits() {
        // Test with i = 7 and k = 3, expecting [1, 1, 1]
        let result = HyperCube::new(7, 3);
        let expected: HyperCube = vec![1.into(), 1.into(), 1.into()].into();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hyper_cube_to_vec_large_k() {
        // Test with i = 2 and k = 5, expecting [0, 0, 0, 1, 0]
        let result = HyperCube::new(2, 5);
        let expected: HyperCube = vec![0.into(), 0.into(), 0.into(), 1.into(), 0.into()].into();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_unique_univariate_line_basic() {
        // Test with simple input arrays
        let b: Vec<ScalarField> = vec![1.into(), 2.into()];
        let c: Vec<ScalarField> = vec![3.into(), 4.into()];

        let result = unique_univariate_line(&b, &c);

        let expected: Vec<UniSparsePolynomial<ScalarField>> = vec![
            UniSparsePolynomial::from_coefficients_slice(&[(0, 1.into()), (1, 2.into())]), // 1 + 2*x
            UniSparsePolynomial::from_coefficients_slice(&[(0, 2.into()), (1, 2.into())]), // 2 + 2*x
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_unique_univariate_line_zero() {
        // Test with zero values in arrays
        let b: Vec<ScalarField> = vec![0.into(), 0.into()];
        let c: Vec<ScalarField> = vec![0.into(), 0.into()];

        let result = unique_univariate_line(&b, &c);

        let expected: Vec<UniSparsePolynomial<ScalarField>> = vec![
            UniSparsePolynomial::from_coefficients_slice(&[(0, 0.into()), (1, 0.into())]), // 0 + 0*x
            UniSparsePolynomial::from_coefficients_slice(&[(0, 0.into()), (1, 0.into())]), // 0 + 0*x
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_unique_univariate_line_mixed() {
        // Test with mixed positive and zero values
        let b: Vec<ScalarField> = vec![0.into(), 2.into()];
        let c: Vec<ScalarField> = vec![1.into(), 2.into()];

        let result = unique_univariate_line(&b, &c);

        let expected: Vec<UniSparsePolynomial<ScalarField>> = vec![
            UniSparsePolynomial::from_coefficients_slice(&[(0, 0.into()), (1, 1.into())]), // 0 + 1*x
            UniSparsePolynomial::from_coefficients_slice(&[(0, 2.into()), (1, 0.into())]), // 2 + 0*x
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_unique_univariate_line_different_lengths() {
        // Test when b and c have different lengths
        let b: Vec<ScalarField> = vec![1.into()];
        let c: Vec<ScalarField> = vec![2.into(), 3.into()];

        let result = unique_univariate_line(&b, &c);

        let expected: Vec<UniSparsePolynomial<ScalarField>> = vec![
            UniSparsePolynomial::from_coefficients_slice(&[(0, 1.into()), (1, 1.into())]), // 1 + 1*x
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_sparse_term_wrapper_conversion_success() {
        // Create a SparseTerm with some values
        let term = SparseTerm::new(vec![(3, 1), (8, 1)]);

        // Create a SparseTermWrapper with k = 1
        let wrapper = SparseTermWrapper(term, 1);

        // Convert the wrapper to a Vec<(usize, usize)>
        let result = Vec::<(usize, usize)>::try_from(wrapper);

        // Expected result after subtraction
        let expected = vec![(2, 1), (7, 1)];

        // Check if the conversion was successful and matches the expected result
        assert_eq!(result, Ok(expected));
    }

    #[test]
    fn test_sparse_term_wrapper_conversion_overflow_error() {
        // Create a SparseTerm where one of the values will cause an overflow
        let term = SparseTerm::new(vec![(1, 1), (10, 1)]);

        // Create a SparseTermWrapper with k = 2, which will cause an overflow on the first value
        let wrapper = SparseTermWrapper(term, 2);

        // Try converting the wrapper to a Vec<(usize, usize)>
        let result = Vec::<(usize, usize)>::try_from(wrapper);

        // Check if the result is an error and matches the SubtractWithOverflow error
        assert_eq!(result, Err(PolyError::SubtractWithOverflow));
    }

    #[test]
    fn test_sparse_term_wrapper_empty_term() {
        // Create an empty SparseTerm
        let term = SparseTerm::new(vec![]);

        // Create a SparseTermWrapper with k = 2
        let wrapper = SparseTermWrapper(term, 2);

        // Convert the wrapper to a Vec<(usize, usize)>
        let result = Vec::<(usize, usize)>::try_from(wrapper);

        // Expected result is an empty vector since the SparseTerm was empty
        let expected: Vec<(usize, usize)> = vec![];

        // Check if the conversion was successful and matches the expected result
        assert_eq!(result, Ok(expected));
    }

    #[test]
    fn test_sparse_term_wrapper_no_overflow() {
        // Create a SparseTerm where subtraction won't cause an overflow
        let term = SparseTerm::new(vec![(3, 1), (5, 1)]);

        // Create a SparseTermWrapper with k = 2
        let wrapper = SparseTermWrapper(term, 2);

        // Convert the wrapper to a Vec<(usize, usize)>
        let result = Vec::<(usize, usize)>::try_from(wrapper);

        // Expected result after subtraction
        let expected = vec![(1, 1), (3, 1)];

        // Check if the conversion was successful and matches the expected result
        assert_eq!(result, Ok(expected));
    }

    #[test]
    fn test_mul_basic_case() {
        let term_vec = vec![
            (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
            (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
        ];
        // Create two multivariate polynomials
        let poly1: MVPoly = (MultiPoly::from_coefficients_vec(2, term_vec.clone())).into();

        let poly2: MVPoly = (MultiPoly::from_coefficients_vec(2, term_vec.clone())).into();

        let result = poly1 * poly2;

        // Expected result should be polynomial multiplication output
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
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
        let poly1 = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![(ScalarField::one(), SparseTerm::new(vec![(0, 1)]))],
        ));

        let zero_poly = MVPoly::new(MultiPoly::from_coefficients_vec(2, vec![]));

        let result = poly1 * zero_poly;
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(2, vec![]));

        assert_eq!(result, expected);
    }

    // #[test]
    // fn test_neg_shift_by_k_basic() {
    //     // Polynomial: x^1 + y^1
    //     let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         2,
    //         vec![
    //             (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
    //             (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
    //         ],
    //     ));

    //     // Shift by k=1
    //     let shifted_poly = poly.neg_shift_by_k(1).unwrap();

    //     // Expected result after shift by k=1: just x^1 (the y term shifts out)
    //     let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         1,
    //         vec![
    //             (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term after shift
    //         ],
    //     ));

    //     assert_eq!(shifted_poly, expected);
    // }
    #[test]
    fn test_neg_shift_by_k_basic() {
        // Polynomial: x^2 (index 1) + y^2 (index 1)
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(1, 2)])),
                (ScalarField::one(), SparseTerm::new(vec![(1, 2)])),
            ],
        ));

        // Shift by k=1
        let shifted_poly = poly.neg_shift_by_k(1).unwrap();

        // BigInt([2, 0, 0, 0]) is represented as ScalarField::from(2)
        let coefficient = ScalarField::from(BigInteger256::new([2, 0, 0, 0]));

        // x_0^2 is represented as SparseTerm::new(vec![(0, 2)])
        let term = SparseTerm::new(vec![(0, 2)]); // x_0^2, where 0 is the variable index and 2 is the power

        // Create the multivariate polynomial using from_coefficients_vec
        let multi_poly = MultiPoly::from_coefficients_vec(1, vec![(coefficient, term)]);

        // Wrap it in MVPoly
        let expected = MVPoly::new(multi_poly);

        assert_eq!(shifted_poly, expected);
    }
    #[test]
    fn test_neg_shift_by_k_large_shift() {
        // Polynomial: x^1 + y^1
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=3 (too large for the number of variables)
        let shifted_poly = poly.neg_shift_by_k(3);

        assert!(shifted_poly.is_err());
    }

    #[test]
    fn test_evaluate_variable_basic() {
        // Polynomial: x^1 + y^1
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
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
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
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
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
            ],
        ));

        // Evaluate x=2, leaving z unevaluated
        let evaluation = poly.evaluate_variable(&vec![ScalarField::from(2u32)]);

        // Expected result: x = 2, z term unevaluated
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
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
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
            ],
        ));

        // Sum last 1 variable (z)
        let summed_poly = poly.sum_last_k_var(1);

        // Coefficients from BigInt
        let coeff2 = ScalarField::from(BigInteger256::new([2, 0, 0, 0]));
        let coeff3 = ScalarField::from(BigInteger256::new([3, 0, 0, 0]));

        // Constant term (BigInt([2, 0, 0, 0]))
        let constant_term = (coeff2, SparseTerm::new(vec![])); // Constant term has no variables

        // Linear term for x_1 (BigInt([3, 0, 0, 0]) * x_1)
        let x1_term = (coeff3, SparseTerm::new(vec![(1, 1)])); // x_1^1

        // Linear term for x_0 (BigInt([3, 0, 0, 0]) * x_0)
        let x0_term = (coeff3, SparseTerm::new(vec![(0, 1)])); // x_0^1

        // Create the multivariate polynomial using from_coefficients_vec
        let multi_poly = MultiPoly::from_coefficients_vec(2, vec![constant_term, x1_term, x0_term]);

        // Wrap it in MVPoly
        let expected = MVPoly::new(multi_poly);

        assert_eq!(summed_poly, expected);
    }

    // TODO fix it
    // #[test]
    // fn test_sum_last_k_var_no_var() {
    //     // Polynomial: x^1 + y^1
    //     let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         2,
    //         vec![
    //             (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
    //             (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
    //         ],
    //     ));
    //     println!("original poly: {:?}", poly);
    //     // BigInt([1, 0, 0, 0])  * x_1
    //     // BigInt([1, 0, 0, 0])  * x_0)

    //     // Expected result: same as original
    //     let expected = poly.clone();

    //     // Sum last 0 variable (no summing)
    //     let summed_poly = poly.sum_last_k_var(0);

    //     assert_eq!(summed_poly, expected);
    // }

    // TODO fix it
    // #[test]
    // fn test_sum_last_k_var_large_k() {
    //     // Polynomial: x^1 + y^1
    //     let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         2,
    //         vec![
    //             (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
    //             (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
    //         ],
    //     ));

    //     // Sum last 3 variables (larger than number of variables)
    //     let summed_poly = poly.sum_last_k_var(3);

    //     // Expected result: empty polynomial (all variables summed out)
    //     let expected = MVPoly::new(MultiPoly::from_coefficients_vec(0, vec![]));

    //     assert_eq!(summed_poly, expected);
    // }

    // TODO fix it
    // #[test]
    // fn test_sum_last_k_var_partial() {
    //     // Polynomial: x^1 + y^1 + z^1
    //     let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         3,
    //         vec![
    //             (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
    //             (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
    //             (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // z term
    //         ],
    //     ));

    //     // Sum last 2 variables (y and z)
    //     let summed_poly = poly.sum_last_k_var(2);

    //     // Expected result: just x^1, y and z summed out
    //     let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
    //         1,
    //         vec![
    //             (ScalarField::from(2u32), SparseTerm::new(vec![(0, 1)])), // x term with y and z summed out
    //         ],
    //     ));

    //     assert_eq!(summed_poly, expected);
    // }

    #[test]
    fn test_binary_to_MVPoly_single_input() {
        // Binary: input "10", eval 1
        let inputs = vec!["10".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (x_0 * (1 - x_1))
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MVPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MVPoly_multiple_inputs() {
        // Binary: inputs "10", "01", evals 1, -1
        let inputs = vec!["10".chars(), "01".chars()];
        let evals = vec![ScalarField::one(), -ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial:
        // 1 * (x_0 * (1 - x_1)) + (-1) * ((1 - x_0) * x_1)
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(1, 1)])),
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MVPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MVPoly_all_zeros() {
        // Binary: input "00", eval 1
        let inputs = vec!["00".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (1 - x_0) * (1 - x_1)
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![])),
                (-ScalarField::one(), SparseTerm::new(vec![(0, 1)])),
                (-ScalarField::one(), SparseTerm::new(vec![(1, 1)])),
                (ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)])),
            ],
        ));

        let ml_poly = MVPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_binary_to_MVPoly_all_ones() {
        // Binary: input "11", eval 1
        let inputs = vec!["11".chars()];
        let evals = vec![ScalarField::one()];
        let binary = Binary { inputs, evals };

        // Expected polynomial: 1 * (x_0 * x_1)
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![(ScalarField::one(), SparseTerm::new(vec![(0, 1), (1, 1)]))],
        ));

        let ml_poly = MVPoly::from(binary);

        assert_eq!(ml_poly, expected);
    }

    #[test]
    fn test_shift_by_k_basic() {
        // Polynomial: x^1 + y^1
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=1
        let shifted_poly = poly.shift_by_k(1);

        // Expected result: x^1 -> z^1, y^1 -> (z+1)^1
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            3,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // shifted x term
                (ScalarField::one(), SparseTerm::new(vec![(2, 1)])), // shifted y term
            ],
        ));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_shift_by_k_large_shift() {
        // Polynomial: x^1
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            1,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
            ],
        ));

        // Shift by k=3
        let shifted_poly = poly.shift_by_k(3);

        // Expected result: x^1 -> (z+3)^1
        let expected = MVPoly::new(MultiPoly::from_coefficients_vec(
            4,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(3, 1)])), // shifted x term
            ],
        ));

        assert_eq!(shifted_poly, expected);
    }

    #[test]
    fn test_shift_by_k_zero_shift() {
        // Polynomial: x^1 + y^1
        let poly = MVPoly::new(MultiPoly::from_coefficients_vec(
            2,
            vec![
                (ScalarField::one(), SparseTerm::new(vec![(0, 1)])), // x term
                (ScalarField::one(), SparseTerm::new(vec![(1, 1)])), // y term
            ],
        ));

        // Shift by k=0 (no shift)
        let shifted_poly = poly.shift_by_k(0);

        // Expected result: same as original
        let expected = poly.clone();

        assert_eq!(shifted_poly, expected);
    }

    // TODO Add tests restrict_to_line
}
