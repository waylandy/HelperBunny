import sys
import numpy as np

import torch
import torch.nn.functional as F

def cluster_greedy(X, cutoff, metric='euclidean', device='cpu'):
    metrics = {
        'hamming'   : lambda x1, x2: torch.cdist(x1, x2, p=0.0),
        'cityblock' : lambda x1, x2: torch.cdist(x1, x2, p=1.0),
        'euclidean' : lambda x1, x2: torch.cdist(x1, x2, p=2.0),
        'cosine'    : lambda x1, x2: (1 - F.cosine_similarity(x1.unsqueeze(1), x2, dim=-1)).abs(),
    }
    assert metric in metrics
    
    x = torch.tensor(X.copy(), dtype=torch.float32, device=device)
    indices = torch.arange(x.shape[0], device=device)
    clusters = torch.full([x.shape[0]], -1, device=device)
    total, n_iter, n_cls = indices.shape[0], 0, 0

    while True:
        i = indices.argmin()
        mask = metrics[metric](x[i].unsqueeze(0), x)[0] <= cutoff
        clusters[indices[mask]] = indices[i]
        x = x[~mask]
        indices = indices[~mask]
        n_iter += 1
        n_cls += mask.sum()
        sys.stderr.write(f'{n_cls} / {total} ({n_iter})\r')
        if mask.all():
            sys.stderr.write(f'\n')
            return clusters.cpu().numpy()


