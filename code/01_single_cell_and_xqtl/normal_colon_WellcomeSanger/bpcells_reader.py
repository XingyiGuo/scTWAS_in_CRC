"""Minimal reader for a BPCells on-disk 'packed-float-matrix-v2' directory.

Layout (all files start with an 8-byte type tag, e.g. b'UINT32v1'):
  shape             uint32[2]            (n_rows, n_cols)
  idxptr            uint64[n_cols+1]     CSC column pointers into the value/index streams
  val               float32[nnz]         nonzero values, column-major
  index_starts      uint32[n_blocks]     reference row index per 128-nonzero block
  index_idx         uint32[n_blocks+1]   word offsets into index_data (per 128-block)
  index_idx_offsets uint64[...]          chunk offsets to add to index_idx
  index_data        uint32[...]          BP128-packed zigzag deltas of row indices
  row_names/col_names  newline-separated strings

Row indices are stored as 128-value blocks, bit-packed with the simdcomp/BP128
lane layout (4 interleaved 32-bit lanes), holding zigzag-encoded first-order
deltas; the block's first value is index_starts[b] + unzigzag(delta[0]).
"""

import numpy as np
import os

BLOCK = 128


class PackedMatrixDir:
    def __init__(self, path):
        self.path = path
        with open(os.path.join(path, "version")) as fh:
            self.version = fh.read().strip()
        with open(os.path.join(path, "storage_order")) as fh:
            self.storage_order = fh.read().strip()
        sh = self._read("shape", np.uint32)
        self.n_rows, self.n_cols = int(sh[0]), int(sh[1])
        self.idxptr = self._read("idxptr", np.uint64)
        self.index_starts = self._read("index_starts", np.uint32)
        self.index_idx = self._read("index_idx", np.uint32).astype(np.int64)
        offs = self._read("index_idx_offsets", np.uint64).astype(np.int64)
        # index_idx is stored in chunks; each chunk c gets offs[c] added.
        # offs has one entry per chunk boundary; chunk size is inferred from it.
        if len(offs) > 1:
            n_per = int(offs[1]) if offs[1] > 0 else len(self.index_idx)
        else:
            n_per = len(self.index_idx)
        if n_per < len(self.index_idx):
            add = np.zeros(len(self.index_idx), dtype=np.int64)
            for c in range(len(offs)):
                lo = c * n_per
                hi = min((c + 1) * n_per, len(self.index_idx))
                if lo < hi:
                    add[lo:hi] = offs[c]
            self.index_idx = self.index_idx + add
        self.n_blocks = len(self.index_starts)
        self.nnz = int(self.idxptr[-1])
        self.val = np.memmap(os.path.join(path, "val"), dtype=np.float32,
                             mode="r", offset=8, shape=(self.nnz,))
        nwords = (os.path.getsize(os.path.join(path, "index_data")) - 8) // 4
        self.index_data = np.memmap(os.path.join(path, "index_data"),
                                    dtype=np.uint32, mode="r", offset=8,
                                    shape=(int(nwords),))
        self.row_names = self._names("row_names")
        self.col_names = self._names("col_names")

    def _read(self, name, dtype):
        return np.fromfile(os.path.join(self.path, name), dtype=dtype, offset=8)

    def _names(self, name):
        with open(os.path.join(self.path, name), "rb") as fh:
            txt = fh.read().decode("utf-8")
        return np.array([s for s in txt.split("\n") if s != ""], dtype=object)

    # ---- BP128 block decoding -------------------------------------------
    def _decode_blocks(self, block_ids):
        """Return (len(block_ids), 128) int64 array of row indices."""
        block_ids = np.asarray(block_ids, dtype=np.int64)
        out = np.empty((len(block_ids), BLOCK), dtype=np.int64)
        starts = self.index_idx[block_ids]
        widths = ((self.index_idx[block_ids + 1] - starts) // 4).astype(np.int64)
        for b in np.unique(widths):
            sel = np.nonzero(widths == b)[0]
            b = int(b)
            nb = len(sel)
            # words for these blocks: (nb, b, 4) with W[k, l] = word 4k + l
            idx = (starts[sel][:, None] + np.arange(4 * b)[None, :])
            W = np.asarray(self.index_data[idx.ravel()]).reshape(nb, b, 4).astype(np.uint64)
            mask = np.uint64((1 << b) - 1)
            j = np.arange(32, dtype=np.int64)
            start_bit = j * b
            w0 = start_bit // 32
            shift = (start_bit % 32).astype(np.uint64)
            lo = (W[:, w0, :] >> shift[None, :, None]) & mask
            cross = (start_bit % 32 + b) > 32
            if cross.any():
                wc = np.nonzero(cross)[0]
                hi = (W[:, w0[wc] + 1, :] << (np.uint64(32) - shift[wc])[None, :, None]) & mask
                lo[:, wc, :] |= hi
            # value order: i = 4*j + l  ->  transpose (j, l) to flat
            u = lo.transpose(0, 1, 2).reshape(nb, 32, 4).reshape(nb, BLOCK)
            u = lo.reshape(nb, 32, 4).reshape(nb, BLOCK)
            d = (u >> np.uint64(1)).astype(np.int64) ^ -(u & np.uint64(1)).astype(np.int64)
            out[sel] = self.index_starts[block_ids[sel]][:, None].astype(np.int64) + np.cumsum(d, axis=1)
        return out

    def row_indices(self, lo, hi):
        """Row indices for nonzeros [lo, hi) of the global CSC stream."""
        b0, b1 = lo // BLOCK, (hi - 1) // BLOCK + 1
        dec = self._decode_blocks(np.arange(b0, b1))
        return dec.ravel()[lo - b0 * BLOCK: hi - b0 * BLOCK]

    def get_columns(self, cols, row_subset=None):
        """Extract columns as CSC (indices, values, colptr).

        row_subset: optional boolean mask over rows; kept rows are remapped to
        their position in np.nonzero(row_subset)[0].
        """
        cols = np.asarray(cols, dtype=np.int64)
        if row_subset is not None:
            remap = np.full(self.n_rows, -1, dtype=np.int64)
            remap[np.nonzero(row_subset)[0]] = np.arange(int(row_subset.sum()))
        idx_out, val_out, counts = [], [], np.empty(len(cols), dtype=np.int64)
        for k, c in enumerate(cols):
            lo, hi = int(self.idxptr[c]), int(self.idxptr[c + 1])
            ri = self.row_indices(lo, hi)
            vv = np.asarray(self.val[lo:hi])
            if row_subset is not None:
                keep = remap[ri] >= 0
                ri, vv = remap[ri[keep]], vv[keep]
            idx_out.append(ri.astype(np.int32))
            val_out.append(vv)
            counts[k] = len(ri)
        colptr = np.concatenate([[0], np.cumsum(counts)])
        return np.concatenate(idx_out), np.concatenate(val_out), colptr
