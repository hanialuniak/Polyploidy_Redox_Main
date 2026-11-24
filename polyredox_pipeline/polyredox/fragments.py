"""Fragment overlap statistics from 10x-style fragment files.

For each fragment, x is the number of same-nucleus fragments covering its left end.
Fragments from one genome copy do not overlap, so x cannot exceed the local copy
number apart from mapping artefacts. This follows Takeuchi & Kato (2024, scPloidy).
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

HUMAN_AUTOSOMES = r"^(chr)?([1-9]|1[0-9]|2[0-2])$"
SHIFT = np.uint64(32)


def count_overlaps(codes, starts, ends, n_barcodes, k=6):
    """Return an n_barcodes x k matrix of f_x counts; column k pools x >= k.

    Coordinates are 0-based half-open, as in Cell Ranger fragment files.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    starts = np.asarray(starts, dtype=np.uint64)
    ends = np.asarray(ends, dtype=np.uint64)
    if codes.size == 0:
        return np.zeros((n_barcodes, k), dtype=np.int64)
    start_key = (codes << SHIFT) | starts
    order = np.lexsort((ends, start_key))
    start_key, ends, codes = start_key[order], ends[order], codes[order]
    unique = np.r_[True, (start_key[1:] != start_key[:-1]) | (ends[1:] != ends[:-1])]
    start_key, ends, codes = start_key[unique], ends[unique], codes[unique]
    end_key = np.sort((codes << SHIFT) | ends)
    depth = np.searchsorted(start_key, start_key, side="right") - np.searchsorted(
        end_key, start_key, side="right"
    )
    depth = np.clip(depth, 1, k).astype(np.int64)
    combined = codes.astype(np.int64) * k + depth - 1
    return np.bincount(combined, minlength=n_barcodes * k).reshape(n_barcodes, k)


def load_exclude_regions(path):
    """Read a BED file into {chrom: (merged_starts, merged_ends)}."""
    if not path:
        return {}
    bed = pd.read_csv(
        path, sep="\t", comment="#", header=None, usecols=[0, 1, 2],
        names=["chrom", "start", "end"], dtype={"chrom": str},
    )
    regions = {}
    for chrom, part in bed.groupby("chrom"):
        part = part.sort_values("start")
        starts = part["start"].to_numpy()
        ends = np.maximum.accumulate(part["end"].to_numpy())
        new = np.r_[True, starts[1:] > ends[:-1]]
        last = np.r_[new[1:], True]
        regions[chrom] = (starts[new], ends[last])
    return regions


def in_regions(positions, region):
    region_starts, region_ends = region
    index = np.searchsorted(region_starts, positions, side="right") - 1
    hit = index >= 0
    hit[hit] = positions[hit] < region_ends[index[hit]]
    return hit


def iter_chromosomes(path, barcode_index, autosome_regex, chunksize=5_000_000):
    """Yield (chrom, codes, starts, ends), keeping memory to about one chromosome.

    The file must be coordinate-sorted, so that each chromosome forms one block.
    """
    is_autosome = re.compile(autosome_regex)
    reader = pd.read_csv(
        path, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "barcode"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "barcode": str},
        chunksize=chunksize,
    )
    current, buffer, closed = None, [], set()

    def merged():
        parts = list(zip(*buffer)) if buffer else ([], [], [])
        joined = [np.concatenate(p) if len(p) else np.array([], dtype=np.int64) for p in parts]
        return current, joined[0], joined[1], joined[2]

    for chunk in reader:
        chroms = chunk["chrom"].to_numpy(dtype=object)
        cuts = np.flatnonzero(chroms[1:] != chroms[:-1]) + 1
        edges = [0, *cuts.tolist(), len(chunk)]
        for left, right in zip(edges[:-1], edges[1:]):
            chrom = chroms[left]
            if chrom != current:
                if current is not None:
                    if is_autosome.match(current):
                        yield merged()
                    closed.add(current)
                if chrom in closed:
                    raise ValueError(f"{path}: {chrom} appears in two blocks; sort the file.")
                current, buffer = chrom, []
            if not is_autosome.match(chrom):
                continue
            part = chunk.iloc[left:right]
            codes = barcode_index.get_indexer(part["barcode"].to_numpy(dtype=object))
            keep = codes >= 0
            buffer.append(
                (codes[keep], part["start"].to_numpy()[keep], part["end"].to_numpy()[keep])
            )
    if current is not None and is_autosome.match(current):
        yield merged()


def library_overlaps(path, barcodes, autosome_regex=HUMAN_AUTOSOMES, exclude=None, k=6, chunksize=5_000_000):
    """Return f_x counts for the barcodes in one fragment file."""
    barcode_index = pd.Index(barcodes)
    if not barcode_index.is_unique:
        raise ValueError(f"{path}: duplicate barcodes were requested")
    exclude = exclude or {}
    total = np.zeros((len(barcode_index), k), dtype=np.int64)
    for chrom, codes, starts, ends in iter_chromosomes(path, barcode_index, autosome_regex, chunksize):
        if chrom in exclude and codes.size:
            drop = in_regions(starts, exclude[chrom]) | in_regions(ends - 1, exclude[chrom])
            codes, starts, ends = codes[~drop], starts[~drop], ends[~drop]
        total += count_overlaps(codes, starts, ends, len(barcode_index), k)
    return total
