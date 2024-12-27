use crate::core::error::GraphError;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum GKRError {
    /// graph related error
    #[error("graph error: {0}")]
    GraphError(#[from] GraphError),
}
