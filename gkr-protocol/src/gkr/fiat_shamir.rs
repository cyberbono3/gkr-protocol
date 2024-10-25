use super::mimc::Mimc7;
use crate::poly::UniPoly;
use ark_bn254::Fr as ScalarField;
use ark_ff::Zero;
use std::ops::Deref;

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
#[derive(Debug, Clone)]
pub struct FiatShamir {
    pub r_vec: Vec<ScalarField>,
    pub circuit_input: Vec<ScalarField>,
    pub hash_func: Mimc7,
}

impl FiatShamir {
    pub fn new(input: Vec<ScalarField>, n_rounds: usize) -> Self {
        Self {
            r_vec: vec![],
            circuit_input: input,
            hash_func: Mimc7::new(n_rounds),
        }
    }
    // Use hash-chaining
    pub fn get_r(&mut self, g: UniPoly) -> ScalarField {
        let input = if self.r_vec.len() == 0 {
            Input::new_first(self.circuit_input.clone(), g)
        } else {
            Input::new_subsequent(
                self.circuit_input.clone(),
                self.r_vec.len(),
                self.r_vec.last().unwrap().clone(),
                g,
            )
        };
        let r = self
            .hash_func
            .multi_hash(input.to_field_vec(), &ScalarField::zero());
        self.r_vec.push(r);
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
