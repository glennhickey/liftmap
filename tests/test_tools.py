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

def coarsen_oracle(runs, key, gap, lens):
    """The chaining rule, run by run: per (key member, other member), in key order, extend
    the open chain when the strand matches and both gaps are in [0, gap]; a chain's length
    is its key span; clip at the other member's end; then canonical form of the result."""
    K = 1 if key == 'b' else 0
    per = {}
    for q, a, t, b, L, s in runs:
        km, om, k, o = (t, q, b, a) if K else (q, t, a, b)
        per.setdefault((km, om), []).append((k, o, L, s))
    chains = []
    for (km, om), rs in per.items():
        rs.sort()
        cur = None
        for k, o, L, s in rs:
            if cur and cur[3] == s:
                kg = k - (cur[0] + cur[2])
                og = (cur[1] - (o + L)) if s == '-' else (o - (cur[1] + cur[2]))
                if 0 <= kg <= gap and 0 <= og <= gap:
                    cur = [cur[0], o if s == '-' else cur[1], k + L - cur[0], s]
                    continue
            if cur: chains.append((km, om, *cur))
            cur = [k, o, L, s]
        if cur: chains.append((km, om, *cur))
    pairs = set()
    for km, om, k, o, L, s in chains:
        e = o + L - lens[om]
        if e > 0:
            if e >= L: continue
            if s == '-': k += e
            L -= e
        for i in range(L):
            kk, oo = k + i, (o + L - 1 - i) if s == '-' else o + i
            a, b = (oo, kk) if K else (kk, oo)
            q, t = (om, km) if K else (km, om)
            pairs.add((q, t, 1 if s == '-' else 0, a, b))
    rec = [dict(q=q, t=t, strand='-' if s else '+', pairs=[(a, b)]) for q, t, s, a, b in pairs]
    return canonical_runs(rec)[0]

def lift_merge_check(label, lmap, pairs, seqs, side, rnd, gap, n=300):
    """Lift windows with --max-gap and compare with the merge rule applied to the oracle's
    pieces: per target sequence and strand, in source order, a piece extends an open
    interval when the source gap and the target gap are both in [0, gap]."""
    wins = []
    for i in range(n):
        s = rnd.choice(list(seqs)); w = rnd.choice([10, 200, 3000])
        lo = rnd.randrange(max(1, seqs[s] - w)); wins.append((s, lo, min(seqs[s], lo + w)))
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
        for i, (s, lo, hi) in enumerate(wins):
            f.write(f'{s}\t{lo}\t{hi}\tw{i}\t0\t+\n')
        bed = f.name
    got = sorted(tuple(l.split('\t')[i] for i in (3, 0, 1, 2, 5))
                 for l in run('lift', lmap, bed, '--from', side, '--max-gap', gap).stdout.splitlines())
    os.unlink(bed)
    by_src = {}
    for q, t, st, a, b in pairs:
        sq, sp, tq, tp = (q, a, t, b) if side == 'a' else (t, b, q, a)
        by_src.setdefault(sq, []).append((sp, tq, tp, st))
    want = []
    for i, (s, lo, hi) in enumerate(wins):
        bases = sorted(x for x in by_src.get(s, []) if lo <= x[0] < hi)
        # pieces: maximal runs of consecutive source bases on one target diagonal
        pieces = []
        for sp, tq, tp, st in bases:
            d = (sp + tp) if st else (tp - sp)
            for pc in pieces:
                if pc[1] == tq and pc[3] == st and pc[4] == d and pc[0] + pc[2] == sp:
                    pc[2] += 1; break
            else:
                pieces.append([sp, tq, 1, st, d])
        pieces = [(sp, tq, (d - sp - L + 1) if st else (d + sp), L, st) for sp, tq, L, st, d in pieces]
        pieces.sort(key=lambda p: (p[0], p[2]))
        res = []
        for sp, tq, tp, L, st in pieces:
            for r in res:
                if r['tq'] != tq or r['st'] != st or not r['open']: continue
                sg = sp - r['se']
                tg = (r['ts'] - (tp + L)) if st else (tp - r['te'])
                if 0 <= sg <= gap and 0 <= tg <= gap:
                    r['se'] = sp + L
                    if st: r['ts'] = tp
                    else: r['te'] = tp + L
                    break
            else:
                res.append(dict(tq=tq, st=st, ts=tp, te=tp + L, se=sp + L, open=True))
            for r in res:
                if sp - r['se'] > gap: r['open'] = False
        for r in res:
            want.append((f'w{i}', r['tq'], str(r['ts']), str(r['te']), '-' if r['st'] else '+'))
    check(f'{label}: {len(want)} intervals merged at --max-gap {gap} match the oracle', got == sorted(want))

def best_check(lmap, seqs, side, rnd, n=3000):
    """--best (the paralogy filter) on data with paralogs: single bases lift to at most
    one target, a subset of the unfiltered result, non-empty whenever that is; --best 1
    keeps everything; windows under --best 0 lift to a subset of the unfiltered pieces."""
    def rows(bed, *opt):
        return [tuple(l.split('\t')[i] for i in (3, 0, 1, 2, 5))
                for l in run('lift', lmap, bed, '--from', side, *opt).stdout.splitlines()]
    def bed_of(pts, name):
        f = tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False)
        for i, (sq, lo, hi) in enumerate(pts):
            f.write(f'{sq}\t{lo}\t{hi}\t{name}{i}\t0\t+\n')
        f.close(); return f.name
    names = list(seqs)
    pts = []
    for _ in range(n):
        sq = rnd.choice(names); p = rnd.randrange(seqs[sq]); pts.append((sq, p, p + 1))
    bed = bed_of(pts, 'p')
    allr, strict, loose = rows(bed), rows(bed, '--best', 0), rows(bed, '--best', 1)
    per_all, per_strict = {}, {}
    for r in allr: per_all.setdefault(r[0], []).append(r)
    for r in strict: per_strict.setdefault(r[0], []).append(r)
    multi = sum(1 for v in per_all.values() if len(v) > 1)
    check(f'--best 0: {multi} bases with several targets each keep exactly one',
          all(len(per_strict.get(k, [])) == 1 for k in per_all) and set(strict) <= set(allr))
    check('--best 1 keeps every target', sorted(loose) == sorted(allr))
    os.unlink(bed)
    wins = []
    for _ in range(300):
        sq = rnd.choice(names); lo = rnd.randrange(max(1, seqs[sq] - 500)); wins.append((sq, lo, min(seqs[sq], lo + 500)))
    bed = bed_of(wins, 'w')
    wa, ws = rows(bed), rows(bed, '--best', 0)
    check(f'--best 0 on windows keeps a subset of the pieces ({len(ws)} of {len(wa)})', set(ws) <= set(wa) and len(ws) < len(wa))
    os.unlink(bed)

def http_check(wd, files, qseqs, tseqs, rnd):
    """Serve wd over HTTP with Range support (standard library only) and require that dump,
    lift and info through http:// URLs equal the local file's.  Skipped when bin/liftmap
    was built without HTTP (make HTTP=1)."""
    import http.server, threading
    counts = {'requests': 0}
    class H(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k): super().__init__(*a, directory=wd, **k)
        def log_message(self, *a): pass
        def do_GET(self):
            rng = self.headers.get('Range')
            path = self.translate_path(self.path)
            if not rng or not os.path.isfile(path):
                return super().do_GET()
            size = os.path.getsize(path)
            lo, hi = rng.split('=')[1].split('-')
            lo = int(lo); hi = min(int(hi) if hi else size - 1, size - 1)
            counts['requests'] += 1
            with open(path, 'rb') as f:
                f.seek(lo); body = f.read(hi - lo + 1)
            self.send_response(206)
            self.send_header('Content-Range', f'bytes {lo}-{hi}/{size}')
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Accept-Ranges', 'bytes')
            self.end_headers()
            self.wfile.write(body)
    srv = http.server.ThreadingHTTPServer(('127.0.0.1', 0), H)
    t = threading.Thread(target=srv.serve_forever, daemon=True); t.start()
    base = f'http://127.0.0.1:{srv.server_address[1]}/'
    try:
        probe = run('info', base + os.path.basename(files[0]), ok=False)
        if probe.returncode != 0:
            print('  skip (bin/liftmap built without HTTP; make HTTP=1)')
            return
        for fn in files:
            url = base + os.path.basename(fn)
            check(f'dump over HTTP equals the local file ({os.path.basename(fn)})', dump(url) == dump(fn))
            check(f'info over HTTP equals the local file ({os.path.basename(fn)})',
                  run('info', url).stdout == run('info', fn).stdout)
            with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
                for i in range(200):
                    s = rnd.choice(list(qseqs)); p = rnd.randrange(qseqs[s]); f.write(f'{s}\t{p}\t{p+300}\tr{i}\t0\t+\n')
                bed = f.name
            before = counts['requests']
            a = run('lift', url, bed, '--max-gap', 20).stdout
            n_req = counts['requests'] - before
            check(f'lift over HTTP equals local ({n_req} range requests for 200 intervals)',
                  a == run('lift', fn, bed, '--max-gap', 20).stdout)
            os.unlink(bed)
        check('a missing URL is refused', run('info', base + 'nope.lmap', ok=False).returncode != 0)
    finally:
        srv.shutdown()

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
    info = run('info', L1).stdout
    check('axis names are recorded', 'meta.axis.a\tquery' in info and 'meta.axis.b\ttarget' in info)
    G = os.path.join(wd, 'g.lmap')
    run('from-paf', paf, G, '--group-sep', '.')
    ginfo = [l.split('\t') for l in run('info', G).stdout.splitlines() if l.startswith('group.')]
    want_g = sorted([('group.a', 'gA', f'{len(qseqs)} sequences', f'{sum(qseqs.values())} bp'),
                     ('group.b', 'gB', f'{len(tseqs)} sequences', f'{sum(tseqs.values())} bp')])
    check('--group-sep groups each axis by genome, with totals', sorted(map(tuple, ginfo)) == want_g)
    check('grouping does not change the runs', dump(G) == want)
    pn = os.path.join(wd, 'pansn.paf')
    with open(pn, 'w') as f:
        for q in ('HG002#1#chr1', 'HG002#2#chr1', 'HG003#1#chr2', 'nosep'):
            f.write(f'{q}\t100\t0\t10\t+\tGCA_000001635.9.chr1\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n')
    P1, P2 = os.path.join(wd, 'pn1.lmap'), os.path.join(wd, 'pn2.lmap')
    run('from-paf', pn, P1, '--group-sep', '#', '--group-fields', 2, '--allow-overlap')
    gp = sorted(l.split('\t')[1] for l in run('info', P1).stdout.splitlines() if l.startswith('group.a'))
    check("PanSN '#' 2: HG002#1, HG002#2, HG003#1; a name without separators joins none",
          gp == ['HG002#1', 'HG002#2', 'HG003#1'])
    run('from-paf', pn, P2, '--group-sep', '.', '--group-fields', 2, '--swap', '--allow-overlap')
    gp = [l.split('\t')[1] for l in run('info', P2).stdout.splitlines() if l.startswith('group.a')]
    check("dotted accession '.' 2: GCA_000001635.9", gp == ['GCA_000001635.9'])
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

    lift_merge_check('merged lift a->b', L1, pairs, qseqs, 'a', rnd, 30)
    lift_merge_check('merged lift b->a', L1, pairs, tseqs, 'b', rnd, 30)

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
    best_check(O, tseqs, 'b', rnd)

    print('spilling (the builder over its memory budget):')
    recs = gen_records(rnd, qseqs, tseqs, 400, disjoint_q=False)
    recs += recs[:40]
    want, pairs = canonical_runs(recs)
    paf = os.path.join(wd, 'sp.paf'); write_paf(recs, qseqs, tseqs, paf)
    SP, NS = os.path.join(wd, 'sp.lmap'), os.path.join(wd, 'ns.lmap')
    p = run('from-paf', paf, SP, '--allow-overlap', '--mem', 1, '--tmp-dir', wd)
    run('from-paf', paf, NS, '--allow-overlap')
    check('--mem 1 spills', '(spilled)' in p.stderr)
    check(f'spilled import gives the canonical {len(want)} runs', dump(SP) == want)
    check('and the same runs as the in-memory import', dump(SP) == dump(NS))
    check('spill files are gone afterwards', not [x for x in os.listdir(wd) if x.startswith('lmap.spill')])
    lift_check('spilled a->b', SP, pairs, qseqs, 'a', rnd)

    print('coarsening:')
    recs = gen_records(rnd, qseqs, tseqs, 60, disjoint_q=True)
    base, pairs = canonical_runs(recs)
    paf = os.path.join(wd, 'co.paf'); write_paf(recs, qseqs, tseqs, paf)
    C0 = os.path.join(wd, 'co.lmap'); run('from-paf', paf, C0)
    lens = {**qseqs, **tseqs}
    for key in ('a', 'b'):
        for gap in (0, 50):
            out = os.path.join(wd, f'co.{key}.{gap}.lmap')
            run('coarsen', C0, out, '--max-gap', gap, '--key', key)
            want = coarsen_oracle(base, key, gap, lens)
            check(f'coarsen --key {key} --max-gap {gap}: {len(want)} runs match the oracle', dump(out) == want)
            info = run('info', out).stdout
            check(f'  max_gap {gap} recorded, sequences and metadata carried over',
                  f'meta.max_gap\t{gap}' in info and 'meta.axis.a\tquery' in info)
    check('--max-gap 0 changes nothing', dump(os.path.join(wd, 'co.b.0.lmap')) == base)

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

    print('over HTTP (range requests):')
    http_check(wd, [L1, O], qseqs, tseqs, rnd)

    robust = os.path.join(BIN, 'robust')
    if os.path.exists(robust):
        r = subprocess.run([robust, L1], capture_output=True, text=True)
        check('robust (corruption detection) passes on a generated file', r.returncode == 0)

    print(('ALL PASS' if not fails else 'FAILURES') + f' ({fails} failures)')
    sys.exit(1 if fails else 0)

if __name__ == '__main__':
    main()
