import numpy as np

"""
import markovchain

labels = np.array(['-', 'X'])
transitions = np.array([
    [ 85, 15],
    [ 85, 15],
])
mc = markovchain.MarkovChain(transitions)
print("".join(labels[i] for i in mc.get_samples(100)))

labels = np.array(['-','-','X','X','X','X'])
transitions = np.array([
    [9, 1, 0, 0, 0, 0],
    [0, 9, 1, 0, 0, 0],
    [0, 0, 4, 6, 0, 0],
    [0, 0, 0, 4, 6, 0],
    [0, 0, 0, 0, 5, 5],
    [5, 0, 0, 0, 0, 5]])
mc = markovchain.MarkovChain(transitions)
print("".join(labels[i] for i in mc.get_samples(100)))

"""

class MarkovChain:
    def __init__(self, transition):
        self.P = transition / transition.sum(1, keepdims=True)
        self.steady_state = self.get_steady_state(self.P)
        
    def get_steady_state(self, P):
        dim  = P.shape[0]
        Q    = (P - np.eye(dim))
        ones = np.ones(dim)
        Q    = np.c_[Q, ones]
        QTQ  = np.dot(Q, Q.T)
        bQT  = np.ones(dim)
        return np.linalg.solve(QTQ, bQT)
   
    def get_samples(self, n):
        rng    = np.random.random(n + 1)
        P_c    = self.P.cumsum(axis=1)
        state  = (self.steady_state.cumsum() > rng[0]).argmax() # start state
        states = []
        for i in rng[1:]:
            state = (P_c[state] > i).argmax() # next state
            states += [state]
        return np.array(states)


