"""End-to-end test of bin/liftmap against a brute-force oracle.  Standard library only.

Random alignments are generated as sets of aligned base pairs, then written two
independent ways -- as PAF (CIGARs with =/X/M/I/D and zero-length ops, both strands, one
copy gzipped) and as UCSC chain (built from the base pairs, not from the PAF).  Checks:

  * import of either format gives exactly the canonical run set: the aligned base pairs as
    maximal runs, whatever the record and block boundaries were
  * to-paf / to-chain export re-imports to the same runs (every --max-gap tried)
  * single-base lifts from axis a and from axis b match the oracle, strand included
  * axis-a overlap is refused by default and accepted with --allow-overlap
  * verify passes, and robust's corruption checks pass on a generated file

usage: python3 tests/test_tools.py [bin_dir]
"""
import gzip, os, random, subprocess, sys, tempfile

BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'bin')
LM = os.path.join(BIN, 'liftmap')
fails = 0

def check(what, ok):
    global fails
    print(('  ok   ' if ok else '  FAIL ') + what)
    if not ok:
        fails += 1

def run(*args, ok=True):
    p = subprocess.run([LM] + [str(a) for a in args], capture_output=True, text=True)
    if ok and p.returncode != 0:
        sys.exit(f'liftmap {" ".join(map(str, args))} failed:\n{p.stderr}')
    return p

# ---------------------------------------------------------------- generation

def gen_records(rnd, qseqs, tseqs, n, disjoint_q):
    """Each record: (qname, tname, strand, pairs) with pairs = [(qpos, tpos)] in CIGAR
    order (target ascending).  With disjoint_q no query base is used twice."""
    used = {q: [] for q in qseqs}
    recs = []
    tries = 0
    while len(recs) < n and tries < 100 * n:
        tries += 1
        q = rnd.choice(list(qseqs)); t = rnd.choice(list(tseqs))
        span = rnd.randint(50, 1500)
        if qseqs[q] <= span + 2 or tseqs[t] <= span + 2:
            continue
        q0 = rnd.randrange(qseqs[q] - span); t0 = rnd.randrange(tseqs[t] - span)
        if disjoint_q and any(not (q0 + span <= a or b <= q0) for a, b in used[q]):
            continue
        used[q].append((q0, q0 + span))
        strand = rnd.choice('+-')
        ops = []      # (op, len) in target-ascending order
        qi = ti = 0
        while qi < span and ti < span:
            L = rnd.randint(1, 60)
            L = min(L, span - qi, span - ti)
            for _ in range(rnd.randint(1, 3)):              # split M into =/X/M pieces
                if L <= 0:
                    break
                k = rnd.randint(1, L) if rnd.random() < 0.5 else L
                ops.append((rnd.choice('M=X'), k)); qi += k; ti += k; L -= k
            if rnd.random() < 0.1:
                ops.append((rnd.choice('ID'), 0))           # zero-length op
            g = rnd.random()
            if g < 0.35 and qi < span - 1:
                k = rnd.randint(1, min(30, span - 1 - qi)); ops.append(('I', k)); qi += k
            elif g < 0.7 and ti < span - 1:
                k = rnd.randint(1, min(30, span - 1 - ti)); ops.append(('D', k)); ti += k
        while ops and ops[-1][0] in 'ID':
            op, k = ops.pop()
            if op == 'I': qi -= k
            else: ti -= k
        qlen_used, tlen_used = qi, ti
        pairs = []
        qp = q0 + qlen_used if strand == '-' else q0
        tp = t0
        for op, k in ops:
            if op in 'M=X':
                for i in range(k):
                    if strand == '+': pairs.append((qp + i, tp + i))
                    else:             pairs.append((qp - 1 - i, tp + i))
                qp += k if strand == '+' else -k
                tp += k
            elif op == 'I':
                qp += k if strand == '+' else -k
            else:
                tp += k
        recs.append(dict(q=q, t=t, strand=strand, ops=ops, q0=q0, q1=q0 + qlen_used,
                         t0=t0, t1=t0 + tlen_used, pairs=pairs))
    return recs

def write_paf(recs, qseqs, tseqs, path, gz=False):
    op = gzip.open if gz else open
    with op(path, 'wt') as f:
        for r in recs:
            cg = ''.join(f'{k}{o}' for o, k in r['ops'])
            n = sum(k for o, k in r['ops'] if o in 'M=X')
            f.write(f"{r['q']}\t{qseqs[r['q']]}\t{r['q0']}\t{r['q1']}\t{r['strand']}\t{r['t']}\t"
                    f"{tseqs[r['t']]}\t{r['t0']}\t{r['t1']}\t{n}\t{n}\t60\ttp:A:P\tcg:Z:{cg}\n")

def write_chain(recs, qseqs, tseqs, path):
    """Chain with t = PAF query (so axis a matches), built from the base pairs."""
    with open(path, 'w') as f:
        for cid, r in enumerate(recs, 1):
            tl, ql = qseqs[r['q']], tseqs[r['t']]
            pairs = sorted(r['pairs'])                       # t (= query) ascending
            rev = r['strand'] == '-'
            conv = [(a, (ql - 1 - b) if rev else b) for a, b in pairs]   # q on its strand
            blocks = []
            for a, b in conv:
                if blocks and blocks[-1][0] + blocks[-1][2] == a and blocks[-1][1] + blocks[-1][2] == b:
                    blocks[-1][2] += 1
                else:
                    blocks.append([a, b, 1])
            ts, te = blocks[0][0], blocks[-1][0] + blocks[-1][2]
            qs, qe = blocks[0][1], blocks[-1][1] + blocks[-1][2]
            f.write(f"chain 1000 {r['q']} {tl} + {ts} {te} {r['t']} {ql} {r['strand']} {qs} {qe} {cid}\n")
            for i, (a, b, L) in enumerate(blocks):
                if i + 1 < len(blocks):
                    f.write(f"{L} {blocks[i+1][0] - a - L} {blocks[i+1][1] - b - L}\n")
                else:
                    f.write(f"{L}\n\n")

def canonical_runs(recs):
    """The aligned base-pair set as maximal runs: (a_name, a, b_name, b, len, strand)."""
    pairs = set()
    for r in recs:
        s = 1 if r['strand'] == '-' else 0
        for a, b in r['pairs']:
            pairs.add((r['q'], r['t'], s, a, b))
    runs = []
    by = {}
    for q, t, s, a, b in pairs:
        d = (a + b) if s else (b - a)
        by.setdefault((q, t, s, d), []).append(a)
    for (q, t, s, d), As in by.items():
        As.sort()
        start = prev = As[0]
        for a in As[1:] + [None]:
            if a is not None and a == prev + 1:
                prev = a; continue
            L = prev - start + 1
            b = (d - start - L + 1) if s else (d + start)
            runs.append((q, start, t, b, L, '-' if s else '+'))
            if a is not None:
                start = prev = a
    return sorted(runs), pairs

def dump(lmap):
    out = run('dump', lmap).stdout
    return sorted((f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5])
                  for f in (l.split('\t') for l in out.splitlines()))

def lift_check(label, lmap, pairs, seqs, side, rnd, n=3000):
    """Lift n random single bases from `side` ('a' or 'b') and compare with the oracle."""
    names = list(seqs)
    pts = [(rnd.choice(names), 0) for _ in range(n)]
    pts = [(s, rnd.randrange(seqs[s])) for s, _ in pts]
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
        for i, (s, p) in enumerate(pts):
            f.write(f'{s}\t{p}\t{p+1}\tr{i}\t0\t+\n')
        bed = f.name
    got = sorted(tuple(l.split('\t')[i] for i in (3, 0, 1, 5))
                 for l in run('lift', lmap, bed, '--from', side).stdout.splitlines())
    idx = {}
    for q, t, s, a, b in pairs:
        if side == 'a': idx.setdefault((q, a), []).append((t, b, s))
        else:           idx.setdefault((t, b), []).append((q, a, s))
    want = sorted((f'r{i}', o, str(p2), '-' if s else '+')
                  for i, (sq, p) in enumerate(pts) for o, p2, s in idx.get((sq, p), []))
    os.unlink(bed)
    check(f'{label}: {len(want)} single-base lifts from axis {side} match the oracle', got == want)

# ---------------------------------------------------------------- tests

def main():
    rnd = random.Random(20260922)
    wd = tempfile.mkdtemp(prefix='liftmap_test_')
    qseqs = {f'gA.chr{i}': rnd.randint(5000, 30000) for i in range(1, 5)}
    tseqs = {f'gB.scaf{i}': rnd.randint(5000, 30000) for i in range(1, 6)}

    print('disjoint query (axis a is a function):')
    recs = generate = gen_records(rnd, qseqs, tseqs, 60, disjoint_q=True)
    want, pairs = canonical_runs(recs)
    paf, pafgz, chain = (os.path.join(wd, x) for x in ('a.paf', 'a.paf.gz', 'a.chain'))
    write_paf(recs, qseqs, tseqs, paf); write_paf(recs, qseqs, tseqs, pafgz, gz=True)
    write_chain(recs, qseqs, tseqs, chain)
    L1, L2, L3 = (os.path.join(wd, x) for x in ('p.lmap', 'pz.lmap', 'c.lmap'))
    run('from-paf', paf, L1); run('from-paf', pafgz, L2); run('from-chain', chain, L3)
    check(f'from-paf gives the canonical {len(want)} runs', dump(L1) == want)
    check('from-paf on gzip input is identical', dump(L2) == want)
    check('from-chain (independently encoded) gives the same runs', dump(L3) == want)
    check('verify passes', run('verify', L1).stdout.strip() == 'OK')
    check('order a, no overlap', 'order\ta' in run('info', L1).stdout)
    for gap in (0, 7, 10000):
        for fmt in ('paf', 'chain'):
            x = os.path.join(wd, f'rt.{gap}.{fmt}'); y = x + '.lmap'
            run(f'to-{fmt}', L1, x, '--max-gap', gap); run(f'from-{fmt}', x, y)
            check(f'to-{fmt} --max-gap {gap} re-imports to the same runs', dump(y) == want)
    lift_check('a->b', L1, pairs, qseqs, 'a', rnd)
    lift_check('b->a', L1, pairs, tseqs, 'b', rnd)
    S = os.path.join(wd, 's.lmap')
    run('from-paf', paf, S, '--swap', '--allow-overlap')
    swapped = sorted((t, b, q, a, L, s) for q, a, t, b, L, s in want)
    check('--swap exchanges the axes', dump(S) == sorted(swapped) or dump(S) == swapped)

    print('overlapping query (axis a may overlap):')
    recs = gen_records(rnd, qseqs, tseqs, 60, disjoint_q=False)
    recs += recs[:5]                                          # duplicate records
    want, pairs = canonical_runs(recs)
    paf = os.path.join(wd, 'o.paf'); write_paf(recs, qseqs, tseqs, paf)
    O = os.path.join(wd, 'o.lmap')
    p = run('from-paf', paf, O, ok=False)
    check('refused without --allow-overlap', p.returncode != 0 and 'overlap' in p.stderr)
    run('from-paf', paf, O, '--allow-overlap')
    check(f'with --allow-overlap gives the canonical {len(want)} runs (duplicates collapse)', dump(O) == want)
    check('stored in order b with the overlap flag',
          'order\tb' in run('info', O).stdout and 'axis_a_overlaps\tyes' in run('info', O).stdout)
    for fmt in ('paf', 'chain'):
        x = os.path.join(wd, f'ort.{fmt}'); y = x + '.lmap'
        run(f'to-{fmt}', O, x); run(f'from-{fmt}', x, y, '--allow-overlap')
        check(f'to-{fmt} re-imports to the same runs', dump(y) == want)
    lift_check('a->b', O, pairs, qseqs, 'a', rnd)
    lift_check('b->a', O, pairs, tseqs, 'b', rnd)

    print('malformed input:')
    bad = os.path.join(wd, 'bad.paf')
    with open(bad, 'w') as f:
        f.write('q\t100\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:5M\n')
    check('a CIGAR that does not span the record is refused', run('from-paf', bad, O, ok=False).returncode != 0)
    with open(bad, 'w') as f:
        f.write('q\t100\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60\n')
    check('a PAF record without a CIGAR is refused', run('from-paf', bad, O, ok=False).returncode != 0)
    with open(bad, 'w') as f:
        f.write('q\t100\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n'
                'q\t90\t20\t30\t+\tt\t100\t20\t30\t10\t10\t60\tcg:Z:10M\n')
    check('a sequence given two lengths is refused', run('from-paf', bad, O, ok=False).returncode != 0)

    robust = os.path.join(BIN, 'robust')
    if os.path.exists(robust):
        r = subprocess.run([robust, L1], capture_output=True, text=True)
        check('robust (corruption detection) passes on a generated file', r.returncode == 0)

    print(('ALL PASS' if not fails else 'FAILURES') + f' ({fails} failures)')
    sys.exit(1 if fails else 0)

if __name__ == '__main__':
    main()
