use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum GraphError {
    /// must generate trace first
    #[error("graph trace has not yet been generated")]
    TraceNotGenerated,
    /// node does not exist
    #[error("a queried node does not exist")]
    NodeExistence,
    /// node does not exist in trace
    #[error("a queried trace result does not exist")]
    TraceNodeExistence,
    /// inputs are of wrong length or an id mismatches
    #[error("bad inputs")]
    BadInputs,
    /// inputs
    #[error("a node that should be an input is not an input")]
    NonInput,
    /// graph format related error
    #[error("bad graph")]
    Format,
}
