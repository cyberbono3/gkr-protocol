use std::cmp::max;
use std::collections::{btree_map::Entry, BTreeMap, HashMap, HashSet};

use ark_bn254::Fr as ScalarField;

use super::error::GraphError;
use super::layer::Layer;

#[derive(Clone, Debug, PartialEq, Eq, Hash, Copy)]
pub struct InputValue {
    pub id: usize,
    pub value: ScalarField,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Node {
    Add { id: usize, inputs: [Box<Node>; 2] },
    Mult { id: usize, inputs: [Box<Node>; 2] },
    Input { id: usize },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Graph {
    pub nodes: BTreeMap<usize, Vec<Node>>,
    pub last_trace: HashMap<Node, ScalarField>,
    pub layers: Vec<Layer>,
}

impl Graph {
    pub fn new() -> Self {
        Self {
            nodes: BTreeMap::new(),
            last_trace: HashMap::new(),
            layers: vec![],
        }
    }

    /// Retrieves a node, as specified by idx, from the Graph of bucketed nodes.
    pub fn get_level(&self, node: &Node) -> Result<usize, GraphError> {
        self.nodes
            .iter()
            .find_map(|(&level, nodes)| {
                if nodes.contains(node) {
                    Some(level)
                } else {
                    None
                }
            })
            .ok_or(GraphError::NodeExistence)
    }

    // Insert the node with given idx
    fn insert(&mut self, idx: usize, node: &Node) {
        self.nodes.entry(idx).or_default().push(node.clone());
    }
}

impl TryFrom<Vec<&Node>> for Graph {
    type Error = GraphError;
    fn try_from(mut nodes: Vec<&Node>) -> Result<Graph, Self::Error> {
        // remove duplicates
        let mut seen = HashSet::new();
        nodes.retain(|item| seen.insert(*item));

        let mut graph = Graph::new();

        // next step (we identify inputs)
        let mut labelled = false;
        while !labelled {
            labelled = true;
            for node in &nodes {
                if graph.get_level(node).is_err() {
                    match node {
                        Node::Input { .. } => {
                            graph.insert(0, node);
                        }
                        Node::Add { inputs, .. } | Node::Mult { inputs, .. } => {
                            let (first_input, second_input) = (&&*inputs[0], &&*inputs[1]);
                            if !nodes.contains(first_input) | !nodes.contains(second_input) {
                                return Err(GraphError::NodeExistence);
                            }
                            let first_level = graph.get_level(&inputs[0]);
                            let second_level = graph.get_level(&inputs[1]);
                            if first_level.is_err() || second_level.is_err() {
                                labelled = false;
                            } else {
                                // can safely unwrap as we checked for errors
                                let idx = max(first_level.unwrap(), second_level.unwrap()) + 1;
                                graph.insert(idx, node);
                            }
                        }
                    }
                }
            }
        }

        Ok(graph)
    }
}
