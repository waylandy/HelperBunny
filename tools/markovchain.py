import numpy as np

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
   
    def simulate(self, n):
        rng    = np.random.random(n + 1)
        P_c    = self.P.cumsum(axis=1)
        state  = (self.steady_state.cumsum() > rng[0]).argmax()
        states = []
        for i in rng[1:]:
            state = (P_c[state] > i).argmax()
            states += [state]
        return np.array(states)

