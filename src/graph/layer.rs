use crate::poly::*;
use ark_bn254::Fr as ScalarField;

use ark_poly::Polynomial;

// #[derive(Debug, PartialEq, Clone)]
// struct BaseLayer {
//     k: usize,
//     prev_k: usize,
//     add: MultiPoly,
//     mult: MultiPoly,
//     w_b: MultiPoly,
//     w_c: MultiPoly,
// }

// impl BaseLayer {
//     pub fn new(k: usize, prev_k: usize, add: MultiPoly,  mult: MultiPoly,  w_b: MultiPoly, w_c: MultiPoly) -> Self {
//         Self {
//             k,
//             prev_k,
//             add,
//             mult,
//             w_b,
//             w_c
//         }
//     }
// }

#[derive(Debug, PartialEq, Clone)]
pub enum Layer {
    OutputLayer {
        k: usize,
        prev_k: usize,
        add: MultiPoly,
        mult: MultiPoly,
        w_b: MultiPoly,
        w_c: MultiPoly,
        d: MultiPoly,
    },
    InterLayer {
        k: usize,
        prev_k: usize,
        add: MultiPoly,
        mult: MultiPoly,
        w_b: MultiPoly,
        w_c: MultiPoly,
    },
    InputLayer {
        k: usize,
        input_ext: MultiPoly,
    },
}

impl Layer {
    pub fn evaluation_ext(&self) -> MultiPoly {
        match self {
            Self::InputLayer { k: _, input_ext } => input_ext.clone(),
            Self::OutputLayer { d, .. } => d.clone(),
            Self::InterLayer { .. } => panic!(),
        }
    }

    pub fn w_ext_gate_eval(&self, r: &Vec<ScalarField>) -> MultiPoly {
        match self {
            Self::InputLayer { k: _, input_ext } => input_ext.clone(),
            Self::InterLayer {
                add,
                mult,
                w_b,
                w_c,
                ..
            }
            | Self::OutputLayer {
                add,
                mult,
                w_b,
                w_c,
                ..
            } => {
                let mut reduced_add_poly = mult_poly(&evaluate_variable(add, r), &(w_b + w_c));
                reduced_add_poly = neg_shift_poly_by_k(&reduced_add_poly, r.len());
                let mut reduced_mult_poly =
                    mult_poly(&evaluate_variable(mult, r), &mult_poly(&w_b, &w_c));
                reduced_mult_poly = neg_shift_poly_by_k(&reduced_mult_poly, r.len());
                reduced_add_poly + reduced_mult_poly
            }
        }
    }

    pub fn w_ext(&self) -> MultiPoly {
        match self {
            Self::InputLayer { k: _, input_ext } => input_ext.clone(),
            Self::InterLayer {
                k: _,
                prev_k,
                add,
                mult,
                w_b,
                w_c,
            }
            | Self::OutputLayer {
                k: _,
                prev_k,
                add,
                mult,
                w_b,
                w_c,
                ..
            } => {
                let f = mult_poly(add, &(w_b + w_c)) + mult_poly(mult, &mult_poly(&w_b, &w_c));
                sum_last_k_var(&f, 2 * prev_k)
            }
        }
    }

    pub fn w_ext_line_restricted_values(
        &self,
        r: &Vec<ScalarField>,
        q0: ScalarField,
        q1: ScalarField,
    ) -> ScalarField {
        match self {
            Self::InputLayer { .. } => panic!(),
            Self::InterLayer { add, mult, .. } | Self::OutputLayer { add, mult, .. } => {
                add.evaluate(r) * (q0 + q1) + mult.evaluate(r) * (q0 * q1)
            }
        }
    }

    pub fn k(&self) -> usize {
        match self {
            Self::InputLayer { k, .. }
            | Self::InterLayer { k, .. }
            | Self::OutputLayer { k, .. } => *k,
        }
    }
}
