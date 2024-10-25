use std::cmp::max;

use ark_bn254::Fr as ScalarField;
use ark_ff::Zero;
use ark_poly::Polynomial;
use rand::Rng;

use super::error::GKRError;
use super::fiat_shamir::FiatShamir;
use super::sumcheck::Prover as SumCheckProver;
use crate::graph::graph::{Graph, InputValue};
use crate::graph::node::Node;
use crate::poly::{unique_univariate_line, PolyError};

#[derive(Debug, Clone)]
pub struct Prover {
    pub graph: Graph,
    pub fs: FiatShamir,
    pub sumcheck_proofs: Vec<SumCheckProver>,
}

impl Prover {
    pub fn new(nodes: Vec<&Node>, mut values: Vec<InputValue>) -> Result<Self, GKRError> {
        let mut graph = Graph::try_from(nodes).map_err(|e| GKRError::GraphError(e))?;
        values.sort_by_key(|x| x.id);
        graph
            .forward(values.clone())
            .map_err(|e| GKRError::GraphError(e))?;
        graph
            .get_multivariate_extension()
            .map_err(|e| GKRError::GraphError(e))?;
        Ok(Self {
            graph,
            fs: FiatShamir::new(values.iter().map(|x| x.value).collect(), 2),
            sumcheck_proofs: vec![],
        })
    }

    pub fn verify(&self) {
        // 1st round
        let last_layer = self.graph.layers.last().unwrap();
        let mut r_i: Vec<ScalarField> = (0..max(last_layer.k(), 1)).map(|_| self.get_r()).collect();

        let mut m_i = last_layer.evaluation_ext().evaluate(&r_i);

        // recursive sumchecks
        for (prev_idx, layer) in self.graph.layers[1..].iter().enumerate().rev() {
            let f_i = layer.w_ext_gate_eval(&r_i).unwrap();
            let mut sumcheck_prover = SumCheckProver::new(&f_i.0);
            sumcheck_prover.verify(m_i);
            let prev_layer = &self.graph.layers[prev_idx];
            let b = sumcheck_prover.r_vec[0..prev_layer.k()].to_vec();
            let c = sumcheck_prover.r_vec[prev_layer.k()..].to_vec();
            assert_eq!(b.len(), c.len());

            let lines = unique_univariate_line(&b, &c);

            assert_eq!(
                b,
                lines
                    .iter()
                    .map(|l| l.evaluate(&ScalarField::zero()))
                    .collect::<Vec<ScalarField>>()
            );

            assert_eq!(
                c,
                lines
                    .iter()
                    .map(|l| l.evaluate(&ScalarField::from(1)))
                    .collect::<Vec<ScalarField>>()
            );

            let poly = prev_layer.w_ext();
            let restricted_poly = poly.restrict_poly_to_line(&lines);

            assert_eq!(
                f_i.0.evaluate(&sumcheck_prover.r_vec),
                // verifier's calc
                layer.w_ext_line_restricted_values(
                    &[r_i.as_slice(), sumcheck_prover.r_vec.as_slice()].concat(),
                    restricted_poly.evaluate(&ScalarField::zero()),
                    restricted_poly.evaluate(&ScalarField::from(1))
                )
            );
            let r_star = self.get_r();
            r_i = lines.iter().map(|l| l.evaluate(&r_star)).collect();
            m_i = restricted_poly.evaluate(&r_star);
        }

        // final round
        let input_layer = self.graph.layers.first().unwrap();
        assert_eq!(m_i, input_layer.w_ext().0.evaluate(&r_i));
    }

    pub fn verify_non_interactive(&mut self) {
        // 1st round
        let last_layer = self.graph.layers.last().unwrap();
        // let mut r_i: Vec<ScalarField> = (0..max(last_layer.k(), 1)).map(|_| self.get_r()).collect();
        let mut r_i = vec![ScalarField::zero()];

        let mut m_i = last_layer.evaluation_ext().evaluate(&r_i);

        // recursive sumchecks
        for (prev_idx, layer) in self.graph.layers[1..].iter().enumerate().rev() {
            let f_i = layer.w_ext_gate_eval(&r_i).unwrap();
            let mut sumcheck_prover = SumCheckProver::new(&f_i.0);
            sumcheck_prover.verify(m_i);
            let prev_layer = &self.graph.layers[prev_idx];
            let b = sumcheck_prover.r_vec[0..prev_layer.k()].to_vec();
            let c = sumcheck_prover.r_vec[prev_layer.k()..].to_vec();
            assert_eq!(b.len(), c.len());

            let lines = unique_univariate_line(&b, &c);

            assert_eq!(
                b,
                lines
                    .iter()
                    .map(|l| l.evaluate(&ScalarField::zero()))
                    .collect::<Vec<ScalarField>>()
            );

            assert_eq!(
                c,
                lines
                    .iter()
                    .map(|l| l.evaluate(&ScalarField::from(1)))
                    .collect::<Vec<ScalarField>>()
            );
            let poly = prev_layer.w_ext();
            let restricted_poly = poly.restrict_poly_to_line(&lines);

            assert_eq!(
                f_i.0.evaluate(&sumcheck_prover.r_vec),
                // verifier's calc
                layer.w_ext_line_restricted_values(
                    &[r_i.as_slice(), sumcheck_prover.r_vec.as_slice()].concat(),
                    restricted_poly.evaluate(&ScalarField::zero()),
                    restricted_poly.evaluate(&ScalarField::from(1))
                )
            );
            let r_star = self.fs.get_r(restricted_poly.clone());
            r_i = lines.iter().map(|l| l.evaluate(&r_star)).collect();
            m_i = restricted_poly.evaluate(&r_star);
        }

        // final round
        let input_layer = self.graph.layers.first().unwrap();
        assert_eq!(m_i, input_layer.w_ext().0.evaluate(&r_i));
    }

    // Verifier procedures
    pub fn get_r(&self) -> ScalarField {
        let mut rng = rand::thread_rng();
        rng.gen()
    }
}
