use thiserror::Error;

use crate::graph::error::GraphError;

#[derive(Error, Debug, PartialEq)]
pub enum GKRError {
    /// graph related error
    #[error("graph error")]
    GraphError(GraphError),
}
