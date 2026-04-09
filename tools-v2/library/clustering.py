import sys
import numpy as np

import torch
import torch.nn.functional as F

"""
import numpy as np
import matplotlib.pyplot as plt

import clustering

x = [
    'AAAAAAAAAA', 'AGGAAAGAAA', 'AGGAAAAGAA', 'AGAAAAAAAA', 'AAAGAAAAAA', 
    'AGGFFAGAAA', 'AGGFFAGABF', 'AAAAAGAAAA', 'AAAAAAAABA', 'AGGFFAGAAA',
    'AGGFFAGAAA', 'AGGFFAGAFF', 'AGGFFFGAFF',
    ]
clustering.get_clusters_greedy(x, 2, metric='hamming')

x = np.random.random((2400, 2))
labels = clustering.get_clusters_greedy(x, 0.2, metric='euclidean', device='cpu')
indices = np.arange(labels.shape[0])
samples = [np.where(labels == i)[0].min() for i in set(labels)]

plt.figure(figsize=(7, 7))
plt.scatter(*x.T, s=3, c="lightgrey")
plt.scatter(*x[samples].T, s=12, c="red")
plt.show()

"""

def get_clusters_greedy(X, cutoff, metric='euclidean', device='cpu'):

    match metric:
        case "hamming":
            if isinstance(X[0], str):
                assert {str} == set(map(type,X))
                X = np.array([list(i) for i in X]).view(np.int32)
            dist = lambda x1, x2: torch.cdist(x1, x2, p=0.0)
        case "cityblock":
            dist = lambda x1, x2: torch.cdist(x1, x2, p=1.0)
        case "euclidean":
            dist = lambda x1, x2: torch.cdist(x1, x2, p=2.0)
        case "cosine":
            dist = lambda x1, x2: (1 - F.cosine_similarity(x1.unsqueeze(1), x2, dim=-1)).abs()
    
    x = torch.tensor(X.copy(), dtype=torch.float32, device=device)
    indices = torch.arange(x.shape[0], device=device)
    clusters = torch.full([x.shape[0]], -1, device=device)
    total, n_iter, n_cls = indices.shape[0], 0, 0

    while True:
        i = indices.argmin()
        mask = dist(x[i].unsqueeze(0), x)[0] <= cutoff
        clusters[indices[mask]] = indices[i]
        x = x[~mask]
        indices = indices[~mask]
        n_iter += 1
        n_cls += mask.sum()
        sys.stderr.write(f'{n_cls} {total} ({n_iter})\r')
        if mask.all():
            sys.stderr.write(f'\n')
            return clusters.cpu().numpy()
