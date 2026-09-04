#!/usr/bin/env python3
"""seq2tm_diagnose.py -- which knob makes Seq2Tm predictions differ from gem/thermal_tm.csv?

A first attempt to re-predict enzymes already in that table agreed only to r = 0.93, with
individual proteins off by up to 3.8 C and a mean difference of -0.08 C: unbiased and
noisy rather than shifted, so not a different model or checkpoint. Three candidate causes:

  batching   upstream code/seq2tm.py uses batch_size=4 in FILE ORDER and feeds the padded
             batch straight to MultiAttModel with no attention mask, so PAD tokens enter
             the attention and pooling. A protein's prediction then depends on what else
             is in its batch. Length-sorting to reduce padding, which is what the first
             attempt did, changes the batches and therefore the numbers.
  device     Apple MPS is not bit-identical to CPU.
  truncation the first attempt cut sequences at 1022 residues (the ESM2 positional limit);
             upstream does not truncate.

This runs the same 200 enzymes through several configurations and reports which one
reproduces the table. Whichever does is the one to use for the full proteomes.

If NO configuration reproduces it, that is itself the finding: the predictor's output is
not reproducible from its inputs alone, and per-enzyme Tm values cannot be quoted -- only
means over many proteins, where unbiased noise averages out.

    python3 seq2tm_diagnose.py gem/thermal_tm.csv --ckpt ... --code ...   (or the env vars)
"""
import argparse, os, sys, math
import numpy as np, pandas as pd, torch

TM_MAX = 100.0


def build(ckpt, code, dev):
    sys.path.insert(0, code)
    from model import MultiAttModel
    import esm
    net = MultiAttModel(320, 3, 4, 4).to(dev)
    net.load_state_dict(torch.load(ckpt, map_location=dev)); net.eval()
    esm2, alpha = esm.pretrained.esm2_t6_8M_UR50D()
    return net, esm2.to(dev).eval(), alpha.get_batch_converter()


def run(pairs, net, esm2, conv, dev, bs=4, sort=False, trunc=None):
    """Return {id: pred_tm}. bs is a fixed batch size, as upstream."""
    p = list(pairs)
    if trunc: p = [(i, s[:trunc]) for i, s in p]
    if sort:  p = sorted(p, key=lambda r: len(r[1]))
    out = {}
    for k in range(math.ceil(len(p) / bs)):
        batch = p[k * bs:(k + 1) * bs]
        _, _, toks = conv(batch)
        toks = toks.to(dev)
        with torch.no_grad():
            emb = esm2(toks, repr_layers=[6], return_contacts=False)['representations'][6]
            pred = net(emb.transpose(1, 2)).cpu().numpy().reshape(-1)
        for (pid, _), v in zip(batch, pred):
            out[pid] = float(v) * TM_MAX
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('table', help='gem/thermal_tm.csv (needs id, sequence, pred_tm)')
    ap.add_argument('--ckpt', default=os.environ.get('SEQ2TOPT_CKPT'))
    ap.add_argument('--code', default=os.environ.get('SEQ2TOPT_CODE'))
    ap.add_argument('-n', type=int, default=200)
    a = ap.parse_args()

    ref = pd.read_csv(a.table).sample(a.n, random_state=0)
    pairs = [(r.id, r.sequence) for r in ref.itertuples()]
    truth = dict(zip(ref.id, ref.pred_tm))
    longest = max(len(s) for _, s in pairs)
    print(f'{len(pairs)} enzymes, longest {longest} aa\n')

    devs = [torch.device('cpu')]
    if getattr(torch.backends, 'mps', None) and torch.backends.mps.is_available():
        devs.append(torch.device('mps'))
    if torch.cuda.is_available():
        devs.append(torch.device('cuda'))

    print(f'{"configuration":<44}{"max|d|":>9}{"mean d":>9}{"r":>10}')
    print('-' * 72)
    best = None
    for dev in devs:
        net, esm2, conv = build(a.ckpt, a.code, dev)
        cfgs = [
            (f'{dev.type}, batch=4, file order, no trunc  (UPSTREAM)', dict(bs=4)),
            (f'{dev.type}, batch=1, file order, no trunc', dict(bs=1)),
            (f'{dev.type}, batch=4, length-sorted, no trunc', dict(bs=4, sort=True)),
            (f'{dev.type}, batch=4, file order, trunc 1022', dict(bs=4, trunc=1022)),
        ]
        for name, kw in cfgs:
            got = run(pairs, net, esm2, conv, dev, **kw)
            d = np.array([got[i] - truth[i] for i in truth])
            r = np.corrcoef([truth[i] for i in truth], [got[i] for i in truth])[0, 1]
            flag = '   <== reproduces' if np.abs(d).max() < 0.05 else ''
            print(f'{name:<44}{np.abs(d).max():>9.4f}{d.mean():>+9.4f}{r:>10.6f}{flag}')
            if best is None or np.abs(d).max() < best[1]:
                best = (name, float(np.abs(d).max()))

    print('\n' + '-' * 72)
    if best[1] < 0.05:
        print(f'USE: {best[0]}')
    else:
        print(f'NOTHING REPRODUCES THE TABLE. Closest: {best[0]} at {best[1]:.3f} C.')
        print('Then per-enzyme Tm values are not reproducible from sequence alone, and')
        print('only means over many proteins should be quoted. Report this rather than')
        print('picking whichever configuration happens to look closest.')


if __name__ == '__main__':
    main()
