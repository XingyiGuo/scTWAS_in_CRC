#!/usr/bin/env python3
"""
24 - A streaming reader for R's serialization format, written because there is no R.

The CosMx per-cell expression of this study exists only inside
`cosmx_exp2_seuratObject.rds` (2.1 GB gzip). No R interpreter is available on either
machine and no Python package that reads RDS can be installed (the package index is
blocked by the egress policy), so the file is parsed directly.

RDS is R's own serialization format: a header, then a recursive stream of
(flags, payload) items in XDR (big-endian) encoding. This module implements enough of
it to walk a Seurat object and pull single genes out of the dgCMatrix count matrices
without ever holding the matrix in memory:

  pass 1 (`walk`)     - map the object: every slot, its class, and for each large
                        atomic vector its length and byte offset in the decompressed
                        stream. Gene names and cell names are captured; the count data
                        is skipped.
  pass 2 (`extract`)  - re-stream, and for the chosen matrix read the row-index vector
                        `i` keeping only the positions belonging to the requested genes,
                        read `p` in full (one int per cell), then read the value vector
                        `x` keeping only those positions. Memory stays proportional to
                        the counts of the requested genes, not to the matrix.

Usage:
    python3 24_rds_reader.py map   <file.rds> <out.json>
    python3 24_rds_reader.py genes <file.rds> <out_prefix> GENE1 GENE2 ...
"""
import gzip, struct, sys, json, os, time
import numpy as np

# ------------------------------------------------------------------ SEXP type codes
NILSXP, SYMSXP, LISTSXP, CLOSXP, ENVSXP, PROMSXP, LANGSXP = 0, 1, 2, 3, 4, 5, 6
SPECIALSXP, BUILTINSXP, CHARSXP, LGLSXP = 7, 8, 9, 10
INTSXP, REALSXP, CPLXSXP, STRSXP, DOTSXP = 13, 14, 15, 16, 17
VECSXP, EXPRSXP, BCODESXP, EXTPTRSXP, WEAKREFSXP, RAWSXP, S4SXP = 19, 20, 21, 22, 23, 24, 25
ALTREP_SXP, ATTRLANGSXP, ATTRLISTSXP = 238, 240, 239
BASEENV_SXP, EMPTYENV_SXP, BCREPREF, BCREPDEF = 241, 242, 243, 244
GENERICREFSXP, CLASSREFSXP, PERSISTSXP, PACKAGESXP = 245, 246, 247, 248
NAMESPACESXP, BASENAMESPACE_SXP, MISSINGARG_SXP = 249, 250, 251
UNBOUNDVALUE_SXP, GLOBALENV_SXP, NILVALUE_SXP, REFSXP = 252, 253, 254, 255

T0 = time.time()
BIG = 100_000          # vectors at least this long are mapped, not materialised


class Stream:
    """Position-tracking reader over the decompressed byte stream.

    The stream is walked item by item, so a plain gzip read per four-byte integer costs
    more than the decompression itself. Everything therefore goes through one large
    in-process buffer and `struct.unpack_from`, which is what makes a full pass over the
    2.1 GB object finish in a couple of minutes rather than tens of them.
    """

    CHUNK = 1 << 23        # 8 MB

    def __init__(self, path):
        self.f = gzip.open(path, "rb")
        self.buf = b""
        self.off = 0
        self.base = 0      # stream position of buf[0]

    @property
    def pos(self):
        return self.base + self.off

    def _fill(self, need):
        if self.off:
            self.buf = self.buf[self.off:]
            self.base += self.off
            self.off = 0
        while len(self.buf) < need:
            b = self.f.read(max(self.CHUNK, need - len(self.buf)))
            if not b:
                raise EOFError(f"eof at {self.pos}")
            self.buf += b

    def read(self, n):
        if self.off + n > len(self.buf):
            self._fill(n)
        out = self.buf[self.off:self.off + n]
        self.off += n
        return out

    def skip(self, n):
        avail = len(self.buf) - self.off
        if n <= avail:
            self.off += n
            return
        n -= avail
        self.base += len(self.buf)
        self.buf, self.off = b"", 0
        while n:
            k = min(n, self.CHUNK)
            b = self.f.read(k)
            if not b:
                raise EOFError
            n -= len(b)
            self.base += len(b)

    def seek_to(self, target):
        if target < self.pos:
            raise ValueError("stream cannot go backwards")
        self.skip(target - self.pos)

    def i32(self):
        if self.off + 4 > len(self.buf):
            self._fill(4)
        v = struct.unpack_from(">i", self.buf, self.off)[0]
        self.off += 4
        return v

    def skip_strings(self, n):
        """Walk n CHARSXPs without decoding them; returns nothing."""
        buf = self.buf
        off = self.off
        upk = struct.unpack_from
        for _ in range(n):
            if off + 8 > len(buf):
                self.off = off
                self._fill(8)
                buf, off = self.buf, self.off
            ln = upk(">i", buf, off + 4)[0]
            off += 8
            if ln > 0:
                if off + ln <= len(buf):
                    off += ln
                else:
                    self.off = off
                    self.skip(ln)
                    buf, off = self.buf, self.off
        self.off = off

    def close(self):
        self.f.close()


class Big:
    """A large atomic vector left on disk: where it is and how long it is."""
    __slots__ = ("type", "n", "at")

    def __init__(self, type_, n, at):
        self.type, self.n, self.at = type_, n, at

    def __repr__(self):
        return f"<Big {self.type} n={self.n} at={self.at}>"


class StopWalk(Exception):
    """Raised to abandon the walk once every matrix we need has been mapped."""


class RDS:
    def __init__(self, path, big=BIG, stop_after=0, verbose=False):
        self.s = Stream(path)
        self.refs = []
        self.big = big
        self.found = []          # dgCMatrix-shaped objects, in the order they close
        self.stop_after = stop_after
        self.verbose = verbose
        self.nitems = 0
        # progress goes to a file, not to stderr: a walk that outlives its time slice is
        # killed with its output buffers unflushed, and the log is how we see how far it
        # got and where the time went.
        self.log = open(os.environ.get("RDS_LOG", "/dev/null"), "w", buffering=1)
        self._read_header()

    def say(self, msg):
        self.log.write(msg + "\n")

    def _note_matrix(self, d):
        """Record a sparse matrix and, if we have enough of them, stop the walk."""
        keys = set(k for k in d if not k.startswith("__"))
        if not {"i", "p", "x", "Dim"} <= keys:
            return
        dim = d.get("Dim")
        dn = d.get("Dimnames")
        rn = dn[0] if isinstance(dn, list) and len(dn) == 2 else None
        rec = {"at_i": d["i"].at if isinstance(d["i"], Big) else None,
               "n_i": int(d["i"].n) if isinstance(d["i"], Big) else int(np.size(d["i"])),
               "at_p": d["p"].at if isinstance(d["p"], Big) else None,
               "n_p": int(d["p"].n) if isinstance(d["p"], Big) else int(np.size(d["p"])),
               "at_x": d["x"].at if isinstance(d["x"], Big) else None,
               "n_x": int(d["x"].n) if isinstance(d["x"], Big) else int(np.size(d["x"])),
               "type_x": int(d["x"].type) if isinstance(d["x"], Big) else REALSXP,
               "dim": [int(v) for v in np.asarray(dim).ravel()] if dim is not None else None,
               "rownames": rn if isinstance(rn, list) else None,
               "at_colnames": dn[1].at if (isinstance(dn, list) and len(dn) == 2
                                           and isinstance(dn[1], Big)) else None,
               "n_colnames": int(dn[1].n) if (isinstance(dn, list) and len(dn) == 2
                                              and isinstance(dn[1], Big)) else None,
               "i_inline": d["i"] if not isinstance(d["i"], Big) else None,
               "p_inline": d["p"] if not isinstance(d["p"], Big) else None,
               "x_inline": d["x"] if not isinstance(d["x"], Big) else None,
               }
        self.found.append(rec)
        if self.verbose:
            self.say(f"  matrix {len(self.found)}: dim={rec['dim']} nnz={rec['n_x']:,} "
                     f"at byte {self.s.pos:,}")
        if self.stop_after and len(self.found) >= self.stop_after:
            raise StopWalk

    def _read_header(self):
        magic = self.s.read(2)
        if magic != b"X\n":
            raise ValueError(f"not an XDR-serialized RDS (magic {magic!r})")
        self.version = self.s.i32()
        self.writer = self.s.i32()
        self.minreader = self.s.i32()
        if self.version >= 3:
            n = self.s.i32()
            self.encoding = self.s.read(n).decode("latin-1")
        else:
            self.encoding = "native"

    # -------------------------------------------------------------- primitives
    def _charsxp(self, flags):
        n = self.s.i32()
        if n == -1:
            return None
        b = self.s.read(n)
        enc = "utf-8" if (flags >> 12) & 8 else "latin-1"
        try:
            return b.decode(enc)
        except UnicodeDecodeError:
            return b.decode("latin-1")

    def _vec_numeric(self, t, keep):
        n = self.s.i32()
        width = 8 if t in (REALSXP,) else (16 if t == CPLXSXP else (1 if t == RAWSXP else 4))
        at = self.s.pos
        if n >= self.big and not keep:
            self.s.skip(n * width)
            return Big(t, n, at)
        raw = self.s.read(n * width)
        if t == REALSXP:
            return np.frombuffer(raw, dtype=">f8").astype(np.float64)
        if t == CPLXSXP:
            return np.frombuffer(raw, dtype=">c16").astype(np.complex128)
        if t == RAWSXP:
            return np.frombuffer(raw, dtype=np.uint8)
        return np.frombuffer(raw, dtype=">i4").astype(np.int32)

    def _pairlist(self, flags):
        out = {}
        order = []
        while True:
            t = flags & 255
            if t == NILVALUE_SXP:
                break
            if flags & (1 << 9):
                self.item()                       # attributes of the pairlist cell
            tag = self.item() if flags & (1 << 10) else None
            val = self.item()
            name = tag if isinstance(tag, str) else f"<{len(order)}>"
            out[name] = val
            order.append(name)
            flags = self.s.i32()
        out["__order__"] = order
        return out

    # -------------------------------------------------------------- main dispatch
    def item(self, keep=False):
        if self.verbose:
            self.nitems += 1
            if not self.nitems % 1_000_000:
                self.say(f"  ... {self.nitems/1e6:.0f}M items, byte {self.s.pos:,}, "
                         f"{time.time()-T0:.0f}s")
        flags = self.s.i32()
        return self._item(flags, keep)

    def _item(self, flags, keep=False):
        t = flags & 255
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))

        if t == NILVALUE_SXP or t == NILSXP:
            return None
        if t == REFSXP:
            idx = flags >> 8
            if idx == 0:
                idx = self.s.i32()
            return self.refs[idx - 1]
        if t in (GLOBALENV_SXP, BASEENV_SXP, EMPTYENV_SXP, UNBOUNDVALUE_SXP,
                 MISSINGARG_SXP, BASENAMESPACE_SXP):
            return f"<env:{t}>"
        if t == SYMSXP:
            name = self.item()
            self.refs.append(name)
            return name
        if t in (PACKAGESXP, NAMESPACESXP):
            v = self.item()
            self.refs.append(v)
            return v
        if t == CHARSXP:
            return self._charsxp(flags)
        if t in (LISTSXP, LANGSXP, DOTSXP, ATTRLANGSXP, ATTRLISTSXP):
            return self._pairlist(flags)
        if t == CLOSXP:
            attr = self.item() if has_attr else None
            self.item()          # environment
            self.item()          # formals
            self.item()          # body
            return "<closure>"
        if t == ENVSXP:
            self.s.i32()                       # locked
            self.refs.append({"__env__": True})
            slot = self.refs[-1]
            self.item()                        # enclos
            frame = self.item()                # frame
            self.item()                        # hashtab
            self.item()                        # attributes
            slot["frame"] = frame
            return slot
        if t == STRSXP:
            n = self.s.i32()
            if n >= self.big and not keep:
                at = self.s.pos
                self.s.skip_strings(n)
                out = Big(STRSXP, n, at)
            else:
                out = [self.item() for _ in range(n)]
            if has_attr:
                self.item()
            return out
        if t in (LGLSXP, INTSXP, REALSXP, CPLXSXP, RAWSXP):
            v = self._vec_numeric(t, keep)
            if has_attr:
                self.item()
            return v
        if t in (VECSXP, EXPRSXP):
            n = self.s.i32()
            vals = [self.item() for _ in range(n)]
            attr = self.item() if has_attr else None
            names = None
            if isinstance(attr, dict):
                nm = attr.get("names")
                if isinstance(nm, list):
                    names = nm
            if names and len(names) == n:
                d = {k: v for k, v in zip(names, vals)}
                d["__order__"] = names
                if isinstance(attr, dict) and "class" in attr:
                    d["__class__"] = attr["class"]
                return d
            return vals
        if t == S4SXP:
            attr = self.item() if has_attr else {}
            if isinstance(attr, dict):
                attr["__s4__"] = True
                self._note_matrix(attr)
            return attr
        if t == ALTREP_SXP:
            info = self.item()
            state = self.item()
            self.item()                        # attributes
            return self._altrep(info, state)
        if t in (EXTPTRSXP, WEAKREFSXP):
            self.item()
            self.item()
            return "<extptr>"
        if t == BCODESXP:
            raise NotImplementedError("byte-compiled code in stream")
        raise NotImplementedError(f"SEXP type {t} at byte {self.s.pos}")

    def _altrep(self, info, state):
        cls = None
        if isinstance(info, dict):
            order = info.get("__order__", [])
            for k in order:
                v = info[k]
                if isinstance(v, str):
                    cls = v
                    break
        if cls in ("compact_intseq", "compact_realseq") and isinstance(state, np.ndarray):
            n, start, step = int(state[0]), state[1], state[2]
            return np.arange(start, start + n * step, step)
        if cls == "wrap_string" or cls == "deferred_string":
            return state
        return state


# ------------------------------------------------------------------ object mapping
def describe(o, path="", out=None, depth=0, maxdepth=12):
    if out is None:
        out = []
    if depth > maxdepth:
        return out
    if isinstance(o, Big):
        out.append({"path": path, "kind": "big", "type": int(o.type), "n": int(o.n), "at": int(o.at)})
    elif isinstance(o, np.ndarray):
        out.append({"path": path, "kind": "vec", "n": int(o.size),
                    "head": [str(x) for x in o[:4]]})
    elif isinstance(o, dict):
        keys = [k for k in o.get("__order__", list(o.keys())) if not k.startswith("__")]
        cls = o.get("class") or o.get("__class__")
        out.append({"path": path, "kind": "s4" if o.get("__s4__") else "list",
                    "class": cls if isinstance(cls, (str, list)) else None,
                    "keys": keys[:40]})
        for k in keys:
            describe(o[k], f"{path}/{k}", out, depth + 1, maxdepth)
    elif isinstance(o, list):
        out.append({"path": path, "kind": "vlist", "n": len(o),
                    "head": [str(x)[:30] for x in o[:3]]})
        for idx, v in enumerate(o[:12]):
            if isinstance(v, (dict, Big)):
                describe(v, f"{path}[{idx}]", out, depth + 1, maxdepth)
    elif isinstance(o, str):
        out.append({"path": path, "kind": "str", "value": o[:60]})
    return out


def find_matrices(tree):
    """Every dgCMatrix-shaped slot in the mapped object."""
    by_path = {e["path"]: e for e in tree}
    hits = []
    for e in tree:
        if e["kind"] not in ("s4", "list"):
            continue
        keys = set(e.get("keys") or [])
        if {"i", "p", "Dim", "x"} <= keys:
            hits.append(e["path"])
    return hits


def main():
    cmd = sys.argv[1]
    path = sys.argv[2]

    if cmd == "map":
        out = sys.argv[3]
        stop_after = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        r = RDS(path, stop_after=stop_after, verbose=True)
        try:
            r.item()
            complete = True
        except StopWalk:
            complete = False
        mats = []
        for k, m in enumerate(r.found):
            mats.append({kk: vv for kk, vv in m.items() if not kk.endswith("_inline")})
        json.dump({"version": r.version, "complete_walk": complete,
                   "n_matrices": len(mats), "matrices": mats},
                  open(out, "w"), indent=1, default=str)
        for k, m in enumerate(mats):
            rn = m["rownames"]
            print(f"matrix {k}: dim={m['dim']} nnz={m['n_x']:,} "
                  f"rownames={'yes' if rn else 'no'}"
                  + (f" e.g. {rn[:3]}" if rn else ""))
        print(f"walk complete: {complete}; wrote {out}")
        r.s.close()
        return

    if cmd == "genes":
        mapfile = sys.argv[3]
        prefix = sys.argv[4]
        genefile = sys.argv[5]
        wanted = [ln.strip() for ln in open(genefile) if ln.strip()]
        extract(path, mapfile, prefix, wanted)
        return

    raise SystemExit("usage: map <rds> <out.json> [stop_after] | genes <rds> <map.json> <prefix> <genelist.txt>")


def _read_strings(s, at, n):
    """Read a block of n CHARSXPs starting at a known offset."""
    s.seek_to(at)
    out = []
    for _ in range(n):
        s.i32()                       # per-string flags (encoding)
        ln = s.i32()
        out.append("" if ln <= 0 else s.read(ln).decode("utf-8", "replace"))
    return out


def extract(path, mapfile, prefix, wanted, want_cells=True):
    """Pull the requested genes out of the largest matrix in the mapped object.

    One forward pass. The dgCMatrix slots are serialised in the order i, p, Dimnames, x,
    so the whole extraction is a single sweep with no seeking backwards:
      i  -> keep only the positions whose row is a requested gene
      p  -> column pointers, one int per cell, small enough to hold
      Dimnames[[2]] -> cell barcodes, so columns can be joined to the metadata
      x  -> values, read only at the positions kept from i
    """
    mp = json.load(open(mapfile))
    m = max(mp["matrices"], key=lambda d: d["n_x"])
    genes = m["rownames"]
    gidx = {g: k for k, g in enumerate(genes)}
    rows = sorted(gidx[g] for g in wanted if g in gidx)
    missing = sorted(set(wanted) - set(genes))
    print(f"panel {len(genes)} genes, {m['dim'][1]:,} cells, nnz {m['n_x']:,}")
    print(f"requested {len(wanted)}, on panel {len(rows)}, absent {len(missing)}")
    if missing:
        print("  absent:", ", ".join(missing[:20]) + (" ..." if len(missing) > 20 else ""))
    if not rows:
        raise SystemExit("none of the requested genes are on the panel")

    mask = np.zeros(len(genes), dtype=bool)
    mask[np.array(rows, dtype=np.int64)] = True

    s = Stream(path)
    CH = 1 << 22

    # ---- i
    s.seek_to(m["at_i"])
    n = m["n_i"]
    hit_pos, hit_row = [], []
    done = 0
    while done < n:
        k = min(CH, n - done)
        a = np.frombuffer(s.read(k * 4), dtype=">i4").astype(np.int32)
        sel = mask[a]
        if sel.any():
            w = np.flatnonzero(sel)
            hit_pos.append(w.astype(np.int64) + done)
            hit_row.append(a[w])
        done += k
        if not (done // CH) % 20:
            print(f"  i: {100*done/n:.0f}%", flush=True)
    hp = np.concatenate(hit_pos) if hit_pos else np.zeros(0, np.int64)
    hr = np.concatenate(hit_row) if hit_row else np.zeros(0, np.int32)
    del hit_pos, hit_row
    print(f"  entries for the requested genes: {len(hp):,}")

    # ---- p
    s.seek_to(m["at_p"])
    p = np.frombuffer(s.read(m["n_p"] * 4), dtype=">i4").astype(np.int64)

    # ---- cell barcodes
    cells = None
    if want_cells and m.get("at_colnames"):
        cells = _read_strings(s, m["at_colnames"], m["n_colnames"])
        print(f"  cell barcodes: {len(cells):,} e.g. {cells[:2]}")

    # ---- x
    s.seek_to(m["at_x"])
    nx, width = m["n_x"], (8 if m["type_x"] == REALSXP else 4)
    dt = ">f8" if width == 8 else ">i4"
    vals = np.zeros(len(hp), dtype=np.float32)
    done, cur = 0, 0
    while done < nx and cur < len(hp):
        k = min(CH, nx - done)
        raw = s.read(k * width)
        end = int(np.searchsorted(hp, done + k))
        if end > cur:
            loc = (hp[cur:end] - done).astype(np.int64)
            vals[cur:end] = np.frombuffer(raw, dtype=dt)[loc].astype(np.float32)
            cur = end
        done += k
    s.close()

    col = (np.searchsorted(p, hp, side="right") - 1).astype(np.int32)
    np.savez_compressed(prefix + "_counts.npz",
                        row=hr, col=col, val=vals,
                        gene_of_row=np.array(genes, dtype=object)[np.array(rows)],
                        rows=np.array(rows, dtype=np.int32),
                        ncell=np.int64(m["dim"][1]))
    if cells is not None:
        with gzip.open(prefix + "_cells.txt.gz", "wt") as fh:
            fh.write("\n".join(cells))
    print(f"wrote {prefix}_counts.npz  ({len(hp):,} entries)")
    print("  value sample:", vals[:8])
    print("  all integral:", bool(np.all(vals == np.round(vals))))


if __name__ == "__main__":
    main()
