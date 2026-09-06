#!/usr/bin/env python3
"""run_seq2tm.py -- Seq2Tm melting-temperature prediction over a whole proteome FASTA.

WHY THIS EXISTS
gem/thermal_tm.csv covers only the 3,042 enzymes that appear in the metabolic models.
The proteostasis hypothesis says the cell fails because its LEAST stable proteins unfold
and overwhelm quality control. Inside the metabolic set there is no such tail: not one
enzyme in any of the four proteomes is predicted to melt below 45 C, and the median
melting points differ by 0.52 C. If a vulnerable tail exists it is outside metabolism --
ribosomes, chaperones, mitochondrial complexes, membrane and cell-wall proteins -- none
of which were ever predicted. This runs the same predictor over the COMPLETE proteomes so
that question can be answered rather than assumed.

SAME MODEL AS THE PAPER
MultiAttModel(emb_dim=320, window=3, n_head=4, n_RD=4) on ESM2-t6-8M layer-6
embeddings, checkpoint model_tm_window.3_r2.0.76.pth from Seq2Topt release v1.0.0
(Qiu et al. 2025, Brief Bioinform 26(2)). Output columns id,pred_tm match
gem/thermal_tm.csv so the two are directly comparable.

DIFFERENCES FROM UPSTREAM code/seq2tm.py, all mechanical
  * reads FASTA rather than CSV, and takes the checkpoint path as an option instead of
    hard-coding '../../large_model_pth/...';
  * batches to a token budget with sequences length-sorted, so a batch of short proteins
    is not padded out to the longest one in the file;
  * appends to the output as it goes, so an interrupted run resumes instead of restarting;
  * uses Apple MPS or CUDA when present.
The arithmetic is unchanged: the head's output is multiplied by 100 exactly as upstream.

RUN
    python3 run_seq2tm.py in.faa out.csv --ckpt /path/model_tm_window.3_r2.0.76.pth \\
                                         --code /path/Seq2Topt/code
or set SEQ2TOPT_CKPT and SEQ2TOPT_CODE and omit them.

CHECK IT FIRST
    python3 run_seq2tm.py --selfcheck gem/thermal_tm.csv --ckpt ... --code ...
re-predicts 200 enzymes already in that table and reports the agreement. If the numbers
do not match, this pipeline is not the one that produced the paper's values and nothing
computed from its output is comparable with them. Do that before the full run.
"""
import argparse, os, sys, time
import numpy as np, pandas as pd, torch
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

TM_MAX = 100.0       # the head predicts a fraction of this; upstream convention
MAXLEN = 1022        # ESM2 positional limit, less BOS/EOS
MINLEN = 10


def device():
    if torch.cuda.is_available(): return torch.device('cuda')
    if getattr(torch.backends, 'mps', None) and torch.backends.mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


def load(ckpt, code):
    sys.path.insert(0, code)
    from model import MultiAttModel
    import esm
    dev = device()
    net = MultiAttModel(320, 3, 4, 4).to(dev)
    net.load_state_dict(torch.load(ckpt, map_location=dev)); net.eval()
    esm2, alpha = esm.pretrained.esm2_t6_8M_UR50D()
    esm2 = esm2.to(dev).eval()
    return dev, net, esm2, alpha.get_batch_converter()


def predict(pairs, dev, net, esm2, conv, budget=8000, log=None,
            fixed_bs=None, sort=True):
    """pairs: [(id, sequence)] -> yields (id, pred_tm).

    BATCHING IS NOT A FREE CHOICE HERE. MultiAttModel receives the padded batch with no
    attention mask, so PAD tokens enter the attention and pooling and a protein's
    prediction depends on which proteins share its batch. Measured on 200 C. auris
    enzymes, batch=4 against batch=1 moves the same protein by up to 5.36 C, sd 1.27 C.

      fixed_bs=1, sort=False   no padding at all. The clean configuration, and the one
                               to use for any new result.
      fixed_bs=4, sort=False   reproduces gem/thermal_tm.csv exactly (verified: max
                               difference 0.0000 C, r = 1.000000 over the first 200
                               C. auris enzymes in file order). Use only to reproduce.
      sort=True                length-sorted token budget: fastest, and NOT comparable
                               with either of the above.
    """
    if sort:
        pairs = sorted(pairs, key=lambda r: len(r[1]))
    i, n, t0 = 0, 0, time.time()
    while i < len(pairs):
        L = len(pairs[i][1]) + 2
        bs = fixed_bs if fixed_bs else max(1, min(64, budget // L))
        batch = pairs[i:i + bs]; i += bs
        _, _, toks = conv(batch)
        toks = toks.to(dev)
        with torch.no_grad():
            emb = esm2(toks, repr_layers=[6], return_contacts=False)['representations'][6]
            pred = net(emb.transpose(1, 2)).cpu().numpy().reshape(-1)
        for (pid, _), v in zip(batch, pred):
            yield pid, float(v) * TM_MAX
        n += len(batch)
        if log and n % 500 < len(batch):
            el = time.time() - t0
            print(f'  {n}/{len(pairs)}  {el/60:.1f} min  '
                  f'eta {(el/max(n,1))*(len(pairs)-n)/60:.1f} min', flush=True)


def read_fasta(p):
    out, name, buf = [], None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if name: out.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: out.append((name, ''.join(buf)))
    # '*' marks a stop codon in some proteome files and is not in the ESM alphabet
    out = [(i, s.replace('*', '').upper()[:MAXLEN]) for i, s in out]
    return [(i, s) for i, s in out if len(s) >= MINLEN]


def selfcheck(path, dev, net, esm2, conv, n=200):
    ref = pd.read_csv(path)
    if 'sequence' not in ref.columns:
        sys.exit('selfcheck needs the sequence column of gem/thermal_tm.csv')
    ref = ref.sample(min(n, len(ref)), random_state=0)
    got = dict(predict([(r.id, r.sequence[:MAXLEN]) for r in ref.itertuples()],
                       dev, net, esm2, conv, fixed_bs=4, sort=False))
    a = ref.pred_tm.values
    b = np.array([got[i] for i in ref.id])
    d = b - a
    print(f'n = {len(a)}')
    print(f'max |difference| = {np.abs(d).max():.4f} C')
    print(f'mean difference   = {d.mean():+.4f} C')
    print(f'correlation       = {np.corrcoef(a, b)[0,1]:.6f}')
    ok = np.abs(d).max() < 0.05
    print('MATCHES the paper\'s table' if ok else
          'DOES NOT MATCH -- do not use this pipeline for comparisons with thermal_tm.csv')
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fasta', nargs='?'); ap.add_argument('out', nargs='?')
    ap.add_argument('--seqs-from', metavar='CSV',
                    help='take id,sequence from this CSV (e.g. gem/thermal_tm.csv) '
                         'instead of a FASTA, preserving its row order')
    ap.add_argument('--ckpt', default=os.environ.get('SEQ2TOPT_CKPT'))
    ap.add_argument('--code', default=os.environ.get('SEQ2TOPT_CODE'))
    ap.add_argument('--selfcheck', metavar='THERMAL_TM_CSV')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--batch-tokens', type=int, default=8000)
    ap.add_argument('--batch-size', type=int, default=1,
                    help='fixed batch size in FILE ORDER. Default 1 = no padding, the '
                         'clean setting. Use 4 only to reproduce gem/thermal_tm.csv. '
                         'Pass 0 for the fast length-sorted mode, which is not '
                         'comparable with either.')
    a = ap.parse_args()
    if not a.ckpt or not a.code:
        sys.exit('need --ckpt and --code (or SEQ2TOPT_CKPT / SEQ2TOPT_CODE)')

    dev, net, esm2, conv = load(a.ckpt, a.code)
    print(f'device {dev}', flush=True)

    if a.selfcheck:
        sys.exit(0 if selfcheck(a.selfcheck, dev, net, esm2, conv) else 1)
    if a.seqs_from:
        # with --seqs-from there is only ONE positional, the output, but argparse fills
        # the first optional positional (fasta) with it. Accept either slot.
        if not a.out and a.fasta:
            a.out, a.fasta = a.fasta, None
        if not a.out: sys.exit('need an output path')
        t = pd.read_csv(a.seqs_from)
        recs = [(str(r.id), str(r.sequence)[:MAXLEN]) for r in t.itertuples()]
    elif a.fasta and a.out:
        recs = read_fasta(a.fasta)
    else:
        sys.exit('need <in.faa> <out.csv>, or --seqs-from CSV <out.csv>')
    if a.limit: recs = recs[:a.limit]
    if os.path.exists(a.out):
        done = set(pd.read_csv(a.out).id.astype(str))
        recs = [r for r in recs if r[0] not in done]
        print(f'resuming: {len(done)} done, {len(recs)} to go', flush=True)
    else:
        with open(a.out, 'w') as fh: fh.write('id,pred_tm\n')
    if not recs:
        print('nothing to do'); return

    n = 0
    for pid, tm in predict(recs, dev, net, esm2, conv, a.batch_tokens, log=True,
                           fixed_bs=(a.batch_size or None), sort=(a.batch_size == 0)):
        with open(a.out, 'a') as fh: fh.write(f'{pid},{tm:.6f}\n')
        n += 1
    print(f'done {n} -> {a.out}', flush=True)


if __name__ == '__main__':
    main()
