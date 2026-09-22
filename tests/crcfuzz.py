"""Mutate an .lmap, then repair every CRC so the mutant is structurally valid but
semantically hostile, and feed it to a reader built with -fsanitize=address,undefined.

Without the CRC repair this finds nothing: the integrity checks reject ~100% of random
mutations before the parsers ever run. With it, ~45% of mutants open and exercise the
directory and chunk decode paths.

usage: crcfuzz.py file.lmap iters [seed] [all|tail]
expects ./pocrun: tests/fuzz_reader.c built with -fsanitize=address,undefined.
"""
import struct, zlib, random, subprocess, sys, os

def parse(buf):
    n=len(buf)
    foff,flen,fcrc = struct.unpack_from('<QQI', buf, n-32)
    secs=[]; p=foff; (nsec,)=struct.unpack_from('<I',buf,p); p+=4
    for _ in range(nsec):
        idlen=buf[p]; p+=1
        sid=buf[p:p+idlen].decode('latin1'); p+=idlen
        so,sl,sraw=struct.unpack_from('<QQQ',buf,p); p+=24
        codec=buf[p]; flags=buf[p+1]; p+=2
        crcpos=p; p+=4
        secs.append((sid,so,sl,crcpos))
    return foff,flen,secs

def repair(buf):
    b=bytearray(buf); n=len(b)
    try: foff,flen,secs=parse(bytes(b))
    except Exception: return None
    # chunk CRCs first (they live in dir.a, whose section CRC is repaired below):
    # 72-byte entries, off u64 at 40, clen u32 at 48, crc u32 at 64
    sd={sid:(so,sl) for sid,so,sl,crcpos in secs}
    if 'runs' in sd and 'dir.a' in sd:
        rso,rsl=sd['runs']; dso,dsl=sd['dir.a']
        if dso+dsl>n: return None
        for e in range(dso, dso+dsl-71, 72):
            off,=struct.unpack_from('<Q',b,e+40); clen,=struct.unpack_from('<I',b,e+48)
            if rso+off+clen<=n:
                struct.pack_into('<I',b,e+64, zlib.crc32(bytes(b[rso+off:rso+off+clen])) & 0xffffffff)
    for sid,so,sl,crcpos in secs:
        if so+sl>n: return None
        struct.pack_into('<I',b,crcpos, zlib.crc32(bytes(b[so:so+sl])) & 0xffffffff)
    struct.pack_into('<I',b,n-32+16, zlib.crc32(bytes(b[foff:foff+flen])) & 0xffffffff)
    struct.pack_into('<I',b,n-32+20, zlib.crc32(bytes(b[0:64])) & 0xffffffff)
    return bytes(b)

src=open(sys.argv[1],'rb').read()
iters=int(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
region=sys.argv[4] if len(sys.argv)>4 else 'all'
rnd=random.Random(seed)
foff,flen,secs=parse(src)
lo,hi = (0,len(src)) if region=='all' else (foff-1500000, len(src))
lo=max(lo,0)
opened=0; crashes=0
for it in range(iters):
    b=bytearray(src)
    for _ in range(rnd.randint(1,10)):
        pos=rnd.randrange(lo,hi); b[pos]=rnd.randrange(256)
    fixed=repair(bytes(b))
    if fixed is None: continue
    open('/tmp/m.lmap','wb').write(fixed)
    r=subprocess.run(['./pocrun','/tmp/m.lmap'],capture_output=True,text=True,
                     env={**os.environ,'ASAN_OPTIONS':'detect_leaks=1'})
    out=(r.stdout or '')+(r.stderr or '')
    if 'opened' in out: opened+=1
    if r.returncode not in (0,1) or 'ERROR' in out or 'runtime error' in out or 'LeakSanitizer' in out:
        crashes+=1
        print(f"  !! iter {it} rc={r.returncode}")
        print('   '+'\n   '.join(out.strip().splitlines()[:12]))
        open(f'/tmp/crash_{seed}_{it}.lmap','wb').write(fixed)
        if crashes>=3: break
print(f"seed={seed} region={region} iters={iters}: {opened} opened, {crashes} crashes/leaks")
