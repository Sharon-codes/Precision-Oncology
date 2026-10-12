"""Synthetic counterexample: high within-sequence cosine is not information loss.

No protein embeddings, trained neural models, or biological labels are used.
The shared offset is orthogonal to a predeclared classification direction.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
rows = []
for seed in range(30):
    rng = np.random.default_rng(seed)
    n, length, dim = 1000, 20, 32
    y = np.tile([0, 1], n // 2)
    x = rng.normal(0, 1, (n, length, dim))
    x[:, :, 0] += 0.6 * (2*y - 1)[:, None]
    original_centered = x - x.mean(axis=1, keepdims=True)
    score = x[:, :, 0].mean(axis=1)
    for offset in [0, 2, 4, 8, 16, 32, 64]:
        shifted = x.copy()
        shifted[:, :, -1] += offset
        normed = shifted / np.linalg.norm(shifted, axis=-1, keepdims=True)
        # Mean off-diagonal cosine without materializing N x L x L matrix.
        cos = ((normed.sum(axis=1)**2).sum(axis=1)-length)/(length*(length-1))
        centered = shifted-shifted.mean(axis=1, keepdims=True)
        err = np.max(np.abs(centered-original_centered))
        assert err < 1e-12
        assert np.array_equal(score, shifted[:, :, 0].mean(axis=1))
        rows.append(dict(seed=seed, offset=offset, mean_token_cosine=float(cos.mean()),
                         auc=float(roc_auc_score(y, score)), centered_max_abs_error=float(err)))
out=ROOT/'results';out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out/'geometry_control.csv',index=False)
meta={'purpose':'synthetic mathematical counterexample, not biological experiment',
      'n_samples_per_seed':1000,'tokens':20,'dimensions':32,'n_seeds':30,
      'class_signal':0.6,'noise_sd':1,'score':'mean over tokens of coordinate 0; fixed a priori',
      'offset_coordinate':31,'claim':'cosine can approach one while fixed classification scores and token-centered representations remain unchanged'}
(out/'geometry_control.json').write_text(json.dumps(meta,indent=2))
print(pd.DataFrame(rows).groupby('offset')[['mean_token_cosine','auc','centered_max_abs_error']].mean().to_string())
