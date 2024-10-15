

Implementation of the GKR protocol as detailed in Justin Thaler's  book  [Proof, Arguments, and Zero Knowledge](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf). 

Protocol prerequisites:
1. Prover P and Verifier V both agree on input x (element of Field F) and some fan-2 arithmetic circuit C over F, where x is an input, and y is an output. 
2. GKR is interactive proof protocol for claim C(x) = y
3. For instance, circuit input size is 4, circuit size is 11 gates.
4. If circuit size is much bigger than input, V saves a bunch of time on it.
5. GKR protocol has a layered structure : 
    - input layer
    - intermediate layers
    - output layer
6. For every layer gate-value finction W is computed
7. Foer every layer MLE of W id computed 

Protocol:
1. P sends a claimed output y. Let y_tilda be an MLE of y.
2. V picks tandom r from field F
3. By Schwarz-Zipper lemma it is "safe" for V to belive C(x) = y 
4. V checks if s_tilda(r) = y_tilda(r) holds, where s_tilda MLE computed by V, y_tilda is claimed output(MLE) by P. If it does hold, a_tilda and y_tilda are equal, so it means, polynomials are equal with the high probability, according to Schwarz-Zipper lemma.



Summary: GKR  protocol is to avoid commiting to intermediate values of a circuit. No commitment or cryptography. Assumption: V knows an input 





