import re,sys
def routines(path):
    txt=open(path).read().split('\n')
    out={}; cur=None; buf=[]
    for ln in txt:
        m=re.match(r'\s*(?:pure |elemental )?(?:subroutine|function)\s+(\w+)',ln)
        e=re.match(r'\s*end (?:subroutine|function)\s+(\w+)',ln)
        if m and cur is None:
            cur=m.group(1); buf=[ln]
        elif cur is not None:
            buf.append(ln)
            if e and e.group(1)==cur:
                out[cur]='\n'.join(buf); cur=None
    return out
files=sys.argv[1:]
rs=[routines(f) for f in files]
common=set(rs[0])
for r in rs[1:]: common &= set(r)
same=sorted(k for k in common if all(rs[0][k]==r[k] for r in rs[1:]))
diff=sorted(k for k in common if not all(rs[0][k]==r[k] for r in rs[1:]))
print("IDENTICAL (%d):"%len(same)); [print("  ",k) for k in same]
print("DIFFERENT (%d):"%len(diff)); [print("  ",k) for k in diff]
for i,f in enumerate(files):
    only=sorted(set(rs[i])-common)
    print("only in",f,":",only)
