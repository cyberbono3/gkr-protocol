use crate::poly::MultiPoly;


#[derive(Debug, PartialEq, Clone)]
pub struct BaseLayer {
    k: usize,
    prev_k: usize,
    add: MultiPoly,
    mult: MultiPoly,
    w_b: MultiPoly,
    w_c: MultiPoly,
}

#[derive(Debug, PartialEq, Clone)]
pub enum Layer {
    OutputLayer {
        base_layer: BaseLayer,
        d: MultiPoly,
    },
    InterLayer {
        base_layer: BaseLayer,
    },
    InputLayer {
        k: usize,
        input_ext: MultiPoly,
    },
}