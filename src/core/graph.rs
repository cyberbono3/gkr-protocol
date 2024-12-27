use std::{
    cmp::max,
    collections::{BTreeMap, HashMap, HashSet},
};

use ark_bn254::Fr as ScalarField;
use ark_ff::Zero;

use super::{error::GraphError, layer::Layer, node::Node};
use crate::poly::*;

#[derive(Clone, Debug, PartialEq, Eq, Hash, Copy)]
pub struct InputValue {
    pub id: usize,
    pub value: ScalarField,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct Graph {
    pub nodes: BTreeMap<usize, Vec<Node>>,
    pub last_trace: HashMap<Node, ScalarField>,
    pub layers: Vec<Layer>,
}

impl Graph {
    pub fn new() -> Self {
        Self::default()
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

    /// Retrieves a node, as specified by idx, from the Graph of bucketed nodes.
    pub fn num_levels(&self) -> usize {
        self.nodes.len()
    }

    fn handle_inputs(op: &Node, trace: &mut HashMap<Node, ScalarField>) -> Result<(), GraphError> {
        match op {
            Node::Add { inputs, .. } => {
                // We expect exactly 2 inputs
                if let (Some(left_node), Some(right_node)) = (inputs.first(), inputs.get(1)) {
                    let left_val = *trace.get(left_node).ok_or(GraphError::TraceNodeExistence)?;
                    let right_val = *trace
                        .get(right_node)
                        .ok_or(GraphError::TraceNodeExistence)?;
                    trace.insert(op.clone(), left_val + right_val);
                } else {
                    return Err(GraphError::NonInput);
                }
            }
            Node::Mult { inputs, .. } => {
                if let (Some(left_node), Some(right_node)) = (inputs.first(), inputs.get(1)) {
                    let left_val = *trace.get(left_node).ok_or(GraphError::TraceNodeExistence)?;
                    let right_val = *trace
                        .get(right_node)
                        .ok_or(GraphError::TraceNodeExistence)?;
                    trace.insert(op.clone(), left_val * right_val);
                } else {
                    return Err(GraphError::NonInput);
                }
            }
            Node::Input { .. } => {
                // Should not appear in levels > 0
                return Err(GraphError::NonInput);
            }
        }
        Ok(())
    }

    /// forward pass for a graph
    pub fn forward(&mut self, mut values: Vec<InputValue>) -> Result<(), GraphError> {
        // remove duplicates
        let mut seen = HashSet::new();
        values.retain(|item| seen.insert(*item));

        let mut trace: HashMap<Node, ScalarField> = HashMap::new();
        for input_node in &self.nodes[&0] {
            match input_node {
                Node::Input { id } => {
                    let relevant_input: Vec<&InputValue> =
                        values.iter().filter(|v| v.id == *id).collect();
                    if relevant_input.is_empty() || relevant_input.len() > 1 {
                        return Err(GraphError::BadInputs);
                    } else {
                        trace.insert(input_node.clone(), relevant_input[0].value);
                    }
                }
                _ => return Err(GraphError::NonInput),
            }
        }

        for i in 1..self.num_levels() {
            for op in &self.nodes[&i] {
                Self::handle_inputs(op, &mut trace)?
            }
        }
        self.last_trace = trace;

        Ok(())
    }

    pub fn get_multivariate_extension(&mut self) -> Result<(), GraphError> {
        fn get_k(n: usize) -> usize {
            let mut k = 0;
            let mut m = n;
            while m > 1 {
                m >>= 1;
                k += 1;
            }
            if n & (n - 1) == 0 {
                k
            } else {
                k + 1
            }
        }

        if self.last_trace.is_empty() {
            return Err(GraphError::TraceNotGenerated);
        }

        let mut layers: Vec<Layer> = Vec::with_capacity(self.nodes.len());
        for (index, layer_nodes) in &self.nodes {
            let k = get_k(layer_nodes.len());
            let layer = if index > &0 {
                let mut add_ext = MultiPoly::zero();
                let mut mult_ext = MultiPoly::zero();
                for (curr, node) in layer_nodes.iter().enumerate() {
                    if let Node::Add { inputs, .. } | Node::Mult { inputs, .. } = node {
                        // index of current node in layer as a binary string
                        let curr_string = format!("{:0k$b}", curr, k = k);

                        // get index of inbound nodes to the current gate
                        let prev_nodes = &self.nodes[&(index - 1)];
                        let prev_k = get_k(prev_nodes.len());
                        let left_index = prev_nodes.iter().position(|r| *r == *inputs[0]).unwrap();
                        let right_index = prev_nodes.iter().position(|r| *r == *inputs[1]).unwrap();

                        // wiring predicates as binary string
                        let left_string = format!("{:0k$b}", left_index, k = prev_k);
                        let right_string = format!("{:0k$b}", right_index, k = prev_k);
                        // total input as current node + inbound node 1 + inbound node 2
                        let input = format!("{}{}{}", curr_string, left_string, right_string);

                        let poly: MVPoly =
                            Binary::new(vec![input.chars()], vec![ScalarField::from(1)]).into();
                        // polynomial_from_binary(vec![input.chars()], vec![ScalarField::from(1)]);

                        if let Node::Add { .. } = node {
                            add_ext = add_ext + poly.0;
                        } else if let Node::Mult { .. } = node {
                            mult_ext = mult_ext + poly.0;
                        }
                    } else {
                        return Err(GraphError::Format);
                    }
                }

                let prev_layer = &layers[*index - 1];

                let poly = prev_layer.w_ext();

                let w_b = poly.shift_by_k(max(k, 1));
                let w_c = poly.shift_by_k(prev_layer.k() + max(k, 1));
                if *index == &self.nodes.len() - 1 {
                    let poly_input = PolyInput::new(
                        (0..layer_nodes.len()).collect(),
                        layer_nodes
                            .iter()
                            .map(|n| *self.last_trace.get(n).unwrap())
                            .collect(),
                        k,
                    );
                    let output_poly: MVPoly = poly_input.into();
                    Layer::OutputLayer {
                        k: get_k(layer_nodes.len()),
                        prev_k: prev_layer.k(),
                        add: add_ext,
                        mult: mult_ext,
                        w_b: w_b.0,
                        w_c: w_c.0,
                        d: output_poly.0,
                    }
                } else {
                    Layer::InterLayer {
                        k: get_k(layer_nodes.len()),
                        prev_k: prev_layer.k(),
                        add: add_ext,
                        mult: mult_ext,
                        w_b: w_b.0,
                        w_c: w_c.0,
                    }
                }
            } else {
                let poly_input = PolyInput::new(
                    (0..layer_nodes.len()).collect(),
                    layer_nodes
                        .iter()
                        .map(|n| *self.last_trace.get(n).unwrap())
                        .collect(),
                    k,
                );
                let input_poly: MVPoly = poly_input.into();
                Layer::InputLayer {
                    k: get_k(layer_nodes.len()),
                    input_ext: input_poly.0,
                }
            };
            layers.push(layer);
        }
        self.layers = layers;
        Ok(())
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
                            let (first_input, second_input) = (&*inputs[0], &*inputs[1]);
                            if !seen.contains(first_input) | !seen.contains(second_input) {
                                return Err(GraphError::NodeExistence);
                            }

                            let first_level = graph.get_level(first_input);
                            let second_level = graph.get_level(second_input);

                            if let (Ok(first_level), Ok(second_level)) = (first_level, second_level)
                            {
                                let idx = max(first_level, second_level) + 1;
                                graph.insert(idx, node);
                            } else {
                                labelled = false;
                            }
                        }
                    }
                }
            }
        }

        Ok(graph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_poly::{
        polynomial::multivariate::{SparseTerm, Term},
        DenseMVPolynomial,
    };

    #[test]
    fn test_add_nodes_to_graph() {
        let input1 = Node::new_input(1);
        let input2 = Node::new_input(2);
        let add_node = Node::new_add(3, input1.clone(), input2.clone());

        let mut graph = Graph::new(); // Assuming you have a Graph::new method
        graph.insert(0, &input1); // Assuming add_node method
        graph.insert(1, &add_node);

        assert!(graph.nodes.contains_key(&0));
        assert!(graph.nodes.contains_key(&1));
        let vec1 = graph.nodes.get(&0).unwrap();
        assert_eq!(vec1.len(), 1);
        let vec2 = graph.nodes.get(&1).unwrap();
        assert_eq!(vec2.len(), 1);

        assert_eq!(input1, vec1.clone().pop().unwrap());
        assert_eq!(add_node, vec2.clone().pop().unwrap());
    }

    #[test]
    fn test_disconnected_nodes() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let mut graph = Graph::new();
        graph.insert(0, &input1);
        graph.insert(0, &input2);

        assert_eq!(graph.get_level(&input1).unwrap(), 0);
        assert_eq!(graph.get_level(&input2).unwrap(), 0);

        // Check that both nodes are at level 0 and are not connected
        assert_eq!(graph.nodes[&0].len(), 2, "Expected two nodes at level 0.");
    }

    #[test]
    fn test_insert_and_get_level() {
        let input_node = Node::new_input(0);
        let add_node = Node::new_add(1, input_node.clone(), input_node.clone());

        let mut graph = Graph::new();
        graph.insert(0, &input_node);
        graph.insert(1, &add_node);

        assert_eq!(graph.get_level(&input_node).unwrap(), 0);
        assert_eq!(graph.get_level(&add_node).unwrap(), 1);
    }

    #[test]
    fn test_forward_with_valid_inputs() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());

        let mut graph = Graph::try_from(vec![&input1, &input2, &add_node]).unwrap();

        let result = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(3),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(5),
            },
        ]);

        assert!(result.is_ok());
        assert_eq!(graph.last_trace[&input1], ScalarField::from(3));
        assert_eq!(graph.last_trace[&input2], ScalarField::from(5));
        assert_eq!(graph.last_trace[&add_node], ScalarField::from(8)); // 3 + 5 = 8
    }

    #[test]
    fn test_complex_graph_multiple_layers() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());
        let mult_node = Node::new_mult(3, add_node.clone(), input1.clone());

        let mut graph = Graph::try_from(vec![&input1, &input2, &add_node, &mult_node]).unwrap();

        // Perform a forward pass
        let result = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(2),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(3),
            },
        ]);

        assert!(result.is_ok());

        // Verify the results in last_trace
        assert_eq!(graph.last_trace[&input1], ScalarField::from(2));
        assert_eq!(graph.last_trace[&input2], ScalarField::from(3));
        assert_eq!(graph.last_trace[&add_node], ScalarField::from(5)); // 2 + 3
        assert_eq!(graph.last_trace[&mult_node], ScalarField::from(10)); // 5 * 2
    }

    #[test]
    fn test_forward_with_invalid_input_id() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());

        let mut graph = Graph::try_from(vec![&input1, &input2, &add_node]).unwrap();

        let result = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(3),
            },
            InputValue {
                id: 2,
                value: ScalarField::from(5),
            }, // Invalid ID
        ]);

        assert_eq!(result, Err(GraphError::BadInputs));
    }

    #[test]
    fn test_forward_with_missing_inputs() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());

        let mut graph = Graph::try_from(vec![&input1, &input2, &add_node]).unwrap();

        let result = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(3),
            },
            // Missing input for id: 1
        ]);

        assert_eq!(result, Err(GraphError::BadInputs));
    }

    #[test]
    fn test_try_from_missing_input_node() {
        // Define some nodes
        let input_node = Node::new_input(1);
        let missing_input_node = Node::new_input(2);

        // Define a node with dependencies on other nodes
        let add_node = Node::new_add(3, input_node.clone(), missing_input_node.clone());

        // Attempt to create a Graph with a missing dependency node (i.e., `missing_input_node` is not in the list)
        let nodes = vec![&input_node, &add_node];

        // The result should be an error because `missing_input_node` is not included in the vector
        let result = Graph::try_from(nodes);
        assert!(result.is_err());
        assert_eq!(result, Err(GraphError::NodeExistence));
    }

    #[test]
    fn test_graph_init_add() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let third_input = Node::new_input(2);
        let add_node = Node::new_add(0, first_input.clone(), second_input.clone());

        let res = Graph::try_from(vec![&first_input, &second_input, &add_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();

        assert_eq!(graph.get_level(&first_input).unwrap(), 0);
        assert_eq!(graph.get_level(&second_input).unwrap(), 0);
        assert_eq!(graph.get_level(&add_node).unwrap(), 1);
        assert_eq!(graph.nodes[&0].len(), 2);

        assert_eq!(
            graph.get_level(&third_input),
            Err(GraphError::NodeExistence)
        );

        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(res.is_ok());

        assert_eq!(graph.last_trace[&first_input], ScalarField::from(1));
        assert_eq!(graph.last_trace[&second_input], ScalarField::from(2));
        assert_eq!(graph.last_trace[&add_node], ScalarField::from(3));

        let res = graph.forward(vec![InputValue {
            id: 0,
            value: ScalarField::from(1),
        }]);
        assert_eq!(res, Err(GraphError::BadInputs));

        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 22,
                value: ScalarField::from(2),
            },
        ]);
        assert_eq!(res, Err(GraphError::BadInputs));
    }

    #[test]
    fn test_get_multivariate_extension() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());

        let mut graph = Graph::try_from(vec![&input1, &input2, &add_node]).unwrap();
        graph
            .forward(vec![
                InputValue {
                    id: 0,
                    value: ScalarField::from(3),
                },
                InputValue {
                    id: 1,
                    value: ScalarField::from(5),
                },
            ])
            .unwrap();

        let result = graph.get_multivariate_extension();
        assert!(result.is_ok());

        // Check that the number of layers is as expected after extension
        assert_eq!(graph.layers.len(), 2); // Assuming one input layer and one output layer
    }

    #[test]
    fn test_graph_wiring_multiple_gates() {
        let input1 = Node::new_input(0);
        let input2 = Node::new_input(1);
        let add_node = Node::new_add(2, input1.clone(), input2.clone());
        let mult_node = Node::new_mult(3, input1.clone(), add_node.clone()); // Mult using input1 and add_node

        let res = Graph::try_from(vec![&input1, &input2, &add_node, &mult_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();

        let result = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(4),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(result.is_ok());

        assert_eq!(graph.last_trace[&add_node], ScalarField::from(6)); // 4 + 2
        assert_eq!(graph.last_trace[&mult_node], ScalarField::from(24)); // 4 * 6
    }

    #[test]
    fn test_graph_init_mult() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let third_input = Node::new_input(2);
        let mult_node = Node::new_mult(0, first_input.clone(), second_input.clone());

        let res = Graph::try_from(vec![&first_input, &second_input, &mult_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();

        assert_eq!(graph.get_level(&first_input).unwrap(), 0);
        assert_eq!(graph.get_level(&second_input).unwrap(), 0);
        assert_eq!(graph.get_level(&mult_node).unwrap(), 1);
        assert_eq!(graph.nodes[&0].len(), 2);

        assert_eq!(
            graph.get_level(&third_input),
            Err(GraphError::NodeExistence)
        );

        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(res.is_ok());

        assert_eq!(graph.last_trace[&first_input], ScalarField::from(1));
        assert_eq!(graph.last_trace[&second_input], ScalarField::from(2));
        assert_eq!(graph.last_trace[&mult_node], ScalarField::from(2));

        let res = graph.forward(vec![InputValue {
            id: 0,
            value: ScalarField::from(1),
        }]);
        assert_eq!(res, Err(GraphError::BadInputs));

        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 22,
                value: ScalarField::from(2),
            },
        ]);
        assert_eq!(res, Err(GraphError::BadInputs));
    }

    #[test]
    fn test_fails_if_not_in_init() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let third_input = Node::new_input(2);
        let add_node = Node::new_add(0, first_input.clone(), third_input);

        let res = Graph::try_from(vec![&first_input, &second_input, &add_node]);
        assert!(res.is_err());
        assert_eq!(res, Err(GraphError::NodeExistence));
    }

    #[test]
    fn test_graph_wiring_add() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let add_node = Node::new_add(0, first_input.clone(), second_input.clone());

        let res = Graph::try_from(vec![&first_input, &second_input, &add_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();
        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(res.is_ok());

        let res = graph.get_multivariate_extension();
        assert!(res.is_ok());

        assert_eq!(
            graph.layers[0],
            Layer::InputLayer {
                k: 1,
                input_ext: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![])),
                        (ScalarField::from(1), SparseTerm::new(vec![(0, 1)]))
                    ],
                )
            }
        );

        let poly: MVPoly = graph.layers[0].evaluation_ext().into();

        assert_eq!(
            graph.layers[1],
            Layer::OutputLayer {
                k: 0,
                prev_k: 1,
                mult: MultiPoly::zero(),
                add: MultiPoly::from_coefficients_vec(
                    3,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![(2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(1, 1), (2, 1)])),
                        (
                            ScalarField::from(1),
                            SparseTerm::new(vec![(0, 1), (1, 1), (2, 1)])
                        )
                    ],
                ),
                w_b: poly.shift_by_k(1).0,
                w_c: poly.shift_by_k(2).0,
                d: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(-3), SparseTerm::new(vec![(0, 1)])),
                        (ScalarField::from(3), SparseTerm::new(vec![]))
                    ],
                ),
            }
        );
    }

    #[test]
    fn test_graph_wiring_mult() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let mult_node = Node::new_mult(0, first_input.clone(), second_input.clone());

        let res = Graph::try_from(vec![&first_input, &second_input, &mult_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();
        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(res.is_ok());

        let res = graph.get_multivariate_extension();
        assert!(res.is_ok());
        res.unwrap();

        assert_eq!(
            graph.layers[0],
            Layer::InputLayer {
                k: 1,
                input_ext: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![])),
                        (ScalarField::from(1), SparseTerm::new(vec![(0, 1)]))
                    ],
                )
            }
        );

        let poly: MVPoly = graph.layers[0].evaluation_ext().into();

        assert_eq!(
            graph.layers[1],
            Layer::OutputLayer {
                k: 0,
                prev_k: 1,
                add: MultiPoly::zero(),
                mult: MultiPoly::from_coefficients_vec(
                    3,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![(2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(1, 1), (2, 1)])),
                        (
                            ScalarField::from(1),
                            SparseTerm::new(vec![(0, 1), (1, 1), (2, 1)])
                        )
                    ],
                ),
                w_b: poly.shift_by_k(1).0,
                w_c: poly.shift_by_k(2).0,
                d: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(-2), SparseTerm::new(vec![(0, 1)])),
                        (ScalarField::from(2), SparseTerm::new(vec![]))
                    ],
                ),
            }
        );
    }

    #[test]
    fn test_graph_wiring_2_gate() {
        let first_input = Node::new_input(0);
        let second_input = Node::new_input(1);
        let add_node = Node::new_add(0, first_input.clone(), second_input.clone());
        let mult_node = Node::new_mult(1, first_input.clone(), second_input.clone());

        let res = Graph::try_from(vec![&first_input, &second_input, &add_node, &mult_node]);
        assert!(res.is_ok());
        let mut graph = res.unwrap();
        let res = graph.forward(vec![
            InputValue {
                id: 0,
                value: ScalarField::from(1),
            },
            InputValue {
                id: 1,
                value: ScalarField::from(2),
            },
        ]);
        assert!(res.is_ok());

        let res = graph.get_multivariate_extension();
        assert!(res.is_ok());

        assert_eq!(
            graph.layers[0],
            Layer::InputLayer {
                k: 1,
                input_ext: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![])),
                        (ScalarField::from(1), SparseTerm::new(vec![(0, 1)]))
                    ],
                )
            }
        );

        let poly: MVPoly = graph.layers[0].evaluation_ext().into();
        assert_eq!(
            graph.layers[1],
            Layer::OutputLayer {
                k: 1,
                prev_k: 1,
                mult: MultiPoly::from_coefficients_vec(
                    3,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                        (
                            ScalarField::from(-1),
                            SparseTerm::new(vec![(0, 1), (1, 1), (2, 1)])
                        )
                    ],
                ),
                add: MultiPoly::from_coefficients_vec(
                    3,
                    vec![
                        (ScalarField::from(1), SparseTerm::new(vec![(2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                        (ScalarField::from(-1), SparseTerm::new(vec![(1, 1), (2, 1)])),
                        (
                            ScalarField::from(1),
                            SparseTerm::new(vec![(0, 1), (1, 1), (2, 1)])
                        )
                    ],
                ),
                w_b: poly.shift_by_k(1).0,
                w_c: poly.shift_by_k(2).0,
                d: MultiPoly::from_coefficients_vec(
                    1,
                    vec![
                        (ScalarField::from(-1), SparseTerm::new(vec![(0, 1)])),
                        (ScalarField::from(3), SparseTerm::new(vec![]))
                    ],
                ),
            }
        );
    }
}
