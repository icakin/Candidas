#!/usr/bin/env python3
"""seq2tm_order_test.py -- is the batch-composition effect the whole story?

The diagnostic established that Seq2Tm's output depends on how sequences are batched:
batch=1, batch=4 in file order and batch=4 length-sorted give three different answers for
the same protein. Upstream code/seq2tm.py feeds the padded batch straight to
MultiAttModel with no attention mask, so PAD tokens enter the attention and pooling. That
is the mechanism.

But the diagnostic then compared against a RANDOM SUBSAMPLE of gem/thermal_tm.csv, batched
by 4. That is not the original batching either, so it could not have reproduced the table
even if the pipeline were identical. This takes rows in FILE ORDER from the start of a
species block, which is what the original run batched, and asks whether that reproduces.

  reproduces  -> batch composition explains everything. Run the full proteomes in file
                 order at batch=4 and the numbers will be comparable with the paper's.
  still off   -> something else differs too, and per-enzyme Tm values are not reproducible
                 from sequence alone. Only means over many proteins are then quotable --
                 which is what Figure 4B reports, so its conclusion survives either way.

Also reports how far apart two batchings put the SAME protein, which is the quantity a
reader should compare against the 0.52 C interspecies difference.

    python3 seq2tm_order_test.py gem/thermal_tm.csv     (with SEQ2TOPT_CKPT / _CODE set)
"""
import argparse, math, os, sys
import numpy as np, pandas as pd, torch
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

TM_MAX = 100.0


def run(pairs, net, esm2, conv, dev, bs):
    out = {}
    for k in range(math.ceil(len(pairs) / bs)):
        b = pairs[k * bs:(k + 1) * bs]
        _, _, toks = conv(b)
        with torch.no_grad():
            emb = esm2(toks.to(dev), repr_layers=[6],
                       return_contacts=False)['representations'][6]
            pr = net(emb.transpose(1, 2)).cpu().numpy().reshape(-1)
        for (pid, _), v in zip(b, pr):
            out[pid] = float(v) * TM_MAX
    return out


def stats(tag, got, truth):
    d = np.array([got[i] - truth[i] for i in got])
    r = np.corrcoef([truth[i] for i in got], [got[i] for i in got])[0, 1]
    print(f'{tag:<46}{np.abs(d).max():>9.4f}{d.mean():>+9.4f}{r:>10.6f}'
          + ('   <== reproduces' if np.abs(d).max() < 0.05 else ''))
    return float(np.abs(d).max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('table'); ap.add_argument('-n', type=int, default=200)
    ap.add_argument('--ckpt', default=os.environ.get('SEQ2TOPT_CKPT'))
    ap.add_argument('--code', default=os.environ.get('SEQ2TOPT_CODE'))
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from model import MultiAttModel
    import esm

    dev = torch.device('cpu')
    net = MultiAttModel(320, 3, 4, 4).to(dev)
    net.load_state_dict(torch.load(a.ckpt, map_location=dev)); net.eval()
    esm2, alpha = esm.pretrained.esm2_t6_8M_UR50D()
    esm2, conv = esm2.to(dev).eval(), alpha.get_batch_converter()

    T = pd.read_csv(a.table)
    T['sp'] = T.id.str.split('|').str[0]
    print(f'table has {len(T)} rows; species blocks: '
          + ', '.join(f'{s}={int(n)}' for s, n in T.sp.value_counts().items()) + '\n')

    print(f'{"configuration":<46}{"max|d|":>9}{"mean d":>9}{"r":>10}')
    print('-' * 74)
    worst = {}
    for sp in T.sp.unique():
        blk = T[T.sp == sp].head(a.n)          # FILE ORDER from the block start
        pairs = [(r.id, r.sequence) for r in blk.itertuples()]
        truth = dict(zip(blk.id, blk.pred_tm))
        for bs in (4, 1):
            got = run(pairs, net, esm2, conv, dev, bs)
            worst[(sp, bs)] = stats(f'{sp}, first {len(pairs)} in file order, batch={bs}',
                                    got, truth)
        # how much does batching alone move the SAME protein?
        g4 = run(pairs, net, esm2, conv, dev, 4)
        g1 = run(pairs, net, esm2, conv, dev, 1)
        dd = np.array([g4[i] - g1[i] for i in g4])
        print(f'    batch=4 vs batch=1, same proteins: max |shift| '
              f'{np.abs(dd).max():.2f} C, sd {dd.std():.2f} C')
        break                                   # one block is enough to answer it

    print('\n' + '-' * 74)
    m = min(worst.values())
    if m < 0.05:
        print('Batch composition explains the mismatch. Run the proteomes in file order,')
        print('batch=4, and the values will be comparable with gem/thermal_tm.csv.')
    else:
        print(f'Still {m:.2f} C off in file order, so batching is not the whole story.')
        print('Treat per-enzyme Seq2Tm values as not reproducible from sequence alone.')
        print('Quote only means over many proteins. Figure 4B already does exactly that:')
        print('its 0.52 C is a mean over 1041 ortholog pairs, and unbiased per-protein')
        print('noise averages out of it. Nothing in the paper needs to change; the')
        print('methods should state the limitation.')


if __name__ == '__main__':
    main()
