use ark_bn254::Fr as ScalarField;
use ark_ff::Zero;
use ark_poly::Polynomial;
use crate::poly::MultiPoly;

#[derive(Clone, Debug, PartialEq, Eq, Hash, Copy)]
pub struct InputValue {
    pub id: usize,
    pub value: ScalarField,
}




