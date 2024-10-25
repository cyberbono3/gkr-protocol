#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Node {
    Add { id: usize, inputs: [Box<Node>; 2] },
    Mult { id: usize, inputs: [Box<Node>; 2] },
    Input { id: usize },
}

impl Node {
    pub fn new_add(id: usize, input1: Node, input2: Node) -> Self {
        Node::Add {
            id,
            inputs: [Box::new(input1), Box::new(input2)],
        }
    }

    pub fn new_mult(id: usize, input1: Node, input2: Node) -> Self {
        Node::Mult {
            id,
            inputs: [Box::new(input1), Box::new(input2)],
        }
    }

    pub fn new_input(id: usize) -> Self {
        Node::Input { id }
    }

    pub fn id(&self) -> usize {
        match *self {
            Node::Add { id, .. } => id,
            Node::Mult { id, .. } => id,
            Node::Input { id } => id,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_input_node() {
        let input_node = Node::new_input(1);
        match input_node {
            Node::Input { id } => assert_eq!(id, 1),
            _ => panic!("Expected Node::Input"),
        }
    }

    #[test]
    fn test_new_add_node() {
        let input1 = Node::new_input(1);
        let input2 = Node::new_input(2);
        let add_node = Node::new_add(3, input1.clone(), input2.clone());

        match add_node {
            Node::Add { id, inputs } => {
                assert_eq!(id, 3);
                assert_eq!(*inputs[0], input1);
                assert_eq!(*inputs[1], input2);
            }
            _ => panic!("Expected Node::Add"),
        }
    }

    #[test]
    fn test_new_mult_node() {
        let input1 = Node::new_input(1);
        let input2 = Node::new_input(2);
        let mult_node = Node::new_mult(4, input1.clone(), input2.clone());

        match mult_node {
            Node::Mult { id, inputs } => {
                assert_eq!(id, 4);
                assert_eq!(*inputs[0], input1);
                assert_eq!(*inputs[1], input2);
            }
            _ => panic!("Expected Node::Mult"),
        }
    }
}
