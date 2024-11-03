use crate::poly::*;
use ark_bn254::Fr as ScalarField;

use ark_poly::Polynomial;

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
    pub fn new_output_layer(
        k: usize,
        prev_k: usize,
        add: MultiPoly,
        mult: MultiPoly,
        w_b: MultiPoly,
        w_c: MultiPoly,
        d: MultiPoly,
    ) -> Self {
        Self::OutputLayer {
            k,
            prev_k,
            add,
            mult,
            w_b,
            w_c,
            d,
        }
    }

    // Constructor for InterLayer
    pub fn new_inter_layer(
        k: usize,
        prev_k: usize,
        add: MultiPoly,
        mult: MultiPoly,
        w_b: MultiPoly,
        w_c: MultiPoly,
    ) -> Self {
        Self::InterLayer {
            k,
            prev_k,
            add,
            mult,
            w_b,
            w_c,
        }
    }

    // Constructor for InputLayer
    pub fn new_input_layer(k: usize, input_ext: MultiPoly) -> Self {
        Self::InputLayer { k, input_ext }
    }
    pub fn evaluation_ext(&self) -> MultiPoly {
        match self {
            Self::InputLayer { k: _, input_ext } => input_ext.clone(),
            Self::OutputLayer { d, .. } => d.clone(),
            Self::InterLayer { .. } => panic!(),
        }
    }

    pub fn w_ext_gate_eval(&self, r: &[ScalarField]) -> Result<MVPoly, PolyError> {
        match self {
            Self::InputLayer { k: _, input_ext } => Ok(input_ext.clone().into()),
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
                let left = MVPoly(add.clone()).evaluate_variable(r);
                let right: MVPoly = (w_b + w_c).into();
                let mut reduced_add_poly: MVPoly = left * right;
                reduced_add_poly = reduced_add_poly.neg_shift_by_k(r.len())?;
                let left = MVPoly(mult.clone()).evaluate_variable(r);
                let right = MVPoly::new(w_b.clone()) * MVPoly::new(w_c.clone());
                let mut reduced_mult_poly = left * right;
                reduced_mult_poly = reduced_mult_poly.neg_shift_by_k(r.len())?;
                Ok(reduced_add_poly + reduced_mult_poly)
            }
        }
    }

    pub fn w_ext(&self) -> MVPoly {
        match self {
            Self::InputLayer { k: _, input_ext } => input_ext.clone().into(),
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
                //let f = mult_poly(add, &(w_b + w_c)) + mult_poly(mult, &mult_poly(&w_b, &w_c));
                let add_poly = MVPoly(add.clone());
                let right_poly1: MVPoly = (w_b + w_c).into();
                let mult_poly = MVPoly(mult.clone());
                let right_poly2: MVPoly = MVPoly(w_b.clone()) * MVPoly(w_c.clone());
                let f = (add_poly + right_poly1) + (mult_poly * right_poly2);
                f.sum_last_k_var(2 * prev_k)
            }
        }
    }

    pub fn w_ext_line_restricted_values(
        &self,
        r: &[ScalarField],
        q0: ScalarField,
        q1: ScalarField,
    ) -> ScalarField {
        match self {
            Self::InputLayer { .. } => panic!(),
            Self::InterLayer { add, mult, .. } | Self::OutputLayer { add, mult, .. } => {
                add.evaluate(&r.to_vec()) * (q0 + q1) + mult.evaluate(&r.to_vec()) * (q0 * q1)
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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_poly::DenseMVPolynomial;

    fn mock_multipoly() -> MultiPoly {
        MultiPoly::from_coefficients_vec(4, vec![])
    }

    #[test]
    fn test_output_layer_initialization() {
        let add = mock_multipoly();
        let mult = mock_multipoly();
        let w_b = mock_multipoly();
        let w_c = mock_multipoly();
        let d = mock_multipoly();

        let output_layer = Layer::new_output_layer(10, 5, add, mult, w_b, w_c, d);

        if let Layer::OutputLayer { k, prev_k, .. } = output_layer {
            assert_eq!(k, 10);
            assert_eq!(prev_k, 5);
        } else {
            panic!("Layer was not correctly initialized as OutputLayer");
        }
    }

    #[test]
    fn test_inter_layer_initialization() {
        let add = mock_multipoly();
        let mult = mock_multipoly();
        let w_b = mock_multipoly();
        let w_c = mock_multipoly();

        let inter_layer = Layer::new_inter_layer(8, 4, add, mult, w_b, w_c);

        if let Layer::InterLayer { k, prev_k, .. } = inter_layer {
            assert_eq!(k, 8);
            assert_eq!(prev_k, 4);
        } else {
            panic!("Layer was not correctly initialized as InterLayer");
        }
    }

    #[test]
    fn test_input_layer_initialization() {
        let prev_k = 3;

        let input_layer = Layer::new_input_layer(6, mock_multipoly());

        if let Layer::InputLayer { k, .. } = input_layer {
            assert_eq!(k, 6);
            assert_eq!(prev_k, 3);
        } else {
            panic!("Layer was not correctly initialized as InputLayer")
        }
    }

    #[test]
    fn test_layer_equality() {
        let add = mock_multipoly();
        let mult = mock_multipoly();
        let w_b = mock_multipoly();
        let w_c = mock_multipoly();

        let layer1 =
            Layer::new_inter_layer(7, 3, add.clone(), mult.clone(), w_b.clone(), w_c.clone());

        let layer2 =
            Layer::new_inter_layer(7, 3, add.clone(), mult.clone(), w_b.clone(), w_c.clone());

        assert_eq!(layer1, layer2);
    }
}
