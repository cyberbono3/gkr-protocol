use std::ops::Deref;

use ark_bn254::Fr as ScalarField;
use light_poseidon::{Poseidon, PoseidonHasher};

use crate::poly::UniPoly;

#[derive(Hash)]
enum Input {
    First {
        circuit_input: Vec<ScalarField>,
        g: UniPoly,
    },
    Subsequent {
        circuit_input: Vec<ScalarField>,
        i: usize,
        prev_r: ScalarField,
        g: UniPoly,
    },
}

impl Input {
    pub fn to_field_vec(&self) -> Vec<ScalarField> {
        match self {
            Self::First { circuit_input, g } => {
                let mut scalar_vec = circuit_input.clone();
                for (_, coeff) in g.deref() {
                    scalar_vec.push(*coeff);
                }
                scalar_vec
            }
            Self::Subsequent {
                circuit_input,
                i,
                prev_r,
                g,
            } => {
                let mut scalar_vec = circuit_input.clone();
                scalar_vec.push(ScalarField::from(*i as u64));
                scalar_vec.push(*prev_r);
                for (_, coeff) in g.deref() {
                    scalar_vec.push(*coeff);
                }
                scalar_vec
            }
        }
    }

    pub fn new_first(circuit_input: Vec<ScalarField>, g: UniPoly) -> Self {
        Self::First { circuit_input, g }
    }

    pub fn new_subsequent(
        circuit_input: Vec<ScalarField>,
        i: usize,
        prev_r: ScalarField,
        g: UniPoly,
    ) -> Self {
        Self::Subsequent {
            circuit_input,
            i,
            prev_r,
            g,
        }
    }
}

// Simulates memory of a single prover instance
#[derive(Clone, Debug)]
pub struct FiatShamir {
    pub r_vec: Vec<ScalarField>,
    pub circuit_input: Vec<ScalarField>,
    pub n_rounds: usize,
}

impl FiatShamir {
    pub fn new(input: Vec<ScalarField>, n_rounds: usize) -> Self {
        Self {
            r_vec: vec![],
            circuit_input: input,
            n_rounds,
        }
    }
    // Use hash-chaining
    pub fn get_r(&mut self, g: UniPoly) -> ScalarField {
        let input = if self.r_vec.is_empty() {
            Input::new_first(self.circuit_input.clone(), g)
        } else {
            Input::new_subsequent(
                self.circuit_input.clone(),
                self.r_vec.len(),
                *self.r_vec.last().unwrap(),
                g,
            )
        };
        let field_vec = input.to_field_vec();
        let mut hasher = Poseidon::<ScalarField>::new_circom(field_vec.len()).unwrap();
        let r = hasher.hash(&field_vec).unwrap();
        self.r_vec.push(r);
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Zero;

    #[test]
    pub fn test_input_to_field_vec() {
        let circuit_input = vec![ScalarField::zero(), ScalarField::from(1)];
        let input = Input::new_first(circuit_input.clone(), UniPoly::zero());

        assert_eq!(input.to_field_vec(), circuit_input);

        let input = Input::new_subsequent(
            circuit_input.clone(),
            32,
            ScalarField::from(22),
            UniPoly::zero(),
        );
        assert_eq!(
            input.to_field_vec(),
            vec![
                ScalarField::zero(),
                ScalarField::from(1),
                ScalarField::from(32),
                ScalarField::from(22)
            ]
        );

        let input = Input::new_first(
            circuit_input.clone(),
            UniPoly::from_coefficients_vec(vec![
                (0, ScalarField::from(1)),
                (1, ScalarField::from(1)),
                (3, ScalarField::from(3)),
                (2, ScalarField::from(2)),
            ]),
        );
        assert_eq!(
            input.to_field_vec(),
            vec![
                ScalarField::zero(),
                ScalarField::from(1),
                ScalarField::from(1),
                ScalarField::from(1),
                ScalarField::from(2),
                ScalarField::from(3)
            ]
        );
    }
}
