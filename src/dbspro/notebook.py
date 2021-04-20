"""
Utility functions for post-processing of data generated using the DBS-Pro pipeline
"""

###########
# FILTERS #
###########
import operator
import time
import pandas as pd
import numpy as np
from itertools import combinations
from umi_tools import UMIClusterer
from typing import Union

from tqdm.notebook import tqdm
from dbspro.utils import jaccard_index


def _filter_wrapper(func):
    """Wrapper to time and monitor Barcode count"""
    def myfilter(*args, **kwargs):
        count_before = len(args[0]['Barcode'].unique())
        t = time.time()
        df = func(*args, **kwargs)
        count_after = len(df['Barcode'].unique())
        count_diff = count_before-count_after
        percent_removed = round(100 * count_diff / count_before, 2)
        print(f"Barcodes = {count_after:,} (-{count_diff:,}, -{percent_removed}%, runtime:{time.time() - t} s)")
        return df
    return myfilter


@_filter_wrapper
def filter_rc(df: pd.DataFrame, threshold: int, opr=operator.gt) -> pd.DataFrame:
    """Filter read count per molecule"""
    print(f"Filtering molecules per readcount {opr.__name__} {threshold}")
    return df[opr(df["ReadCount"], threshold)]


pd.DataFrame.filter_rc = filter_rc


@_filter_wrapper
def filter_rc_sum(df: pd.DataFrame, threshold: int, opr=operator.gt) -> pd.DataFrame:
    """Filter total read count per Barcode"""
    print(f"Filtering barcodes by total read count {opr.__name__} {threshold}")
    return df[opr(df.groupby("Barcode")["ReadCount"].transform('sum'), threshold)]


pd.DataFrame.filter_rc_sum = filter_rc_sum


@_filter_wrapper
def filter_uc(df: pd.DataFrame, threshold: int, opr=operator.gt) -> pd.DataFrame:
    """Filter total UMI count per barcode"""
    print(f"Filtering barcodes by total UMI count {opr.__name__} {threshold}")
    return df[opr(df.groupby("Barcode")["UMI"].transform('count'), threshold)]


pd.DataFrame.filter_uc = filter_uc


@_filter_wrapper
def filter_quantile(df: pd.DataFrame, quantile: float, opr=operator.gt) -> pd.DataFrame:
    """Filter total UMI count per barcode"""
    print(f"Filtering barcodes individually by read count {opr.__name__} quantile {quantile}")

    def above_quantile_threshold(group, quantile):
        return opr(group, np.quantile(group, quantile))

    return df[df.groupby("Barcode")["ReadCount"].apply(above_quantile_threshold, quantile)]


pd.DataFrame.filter_quantile = filter_quantile


@_filter_wrapper
def filter_ratio(df: pd.DataFrame, threshold: int, opr=operator.gt) -> pd.DataFrame:
    """Filter barcodes by ratio of reads to UMIs"""
    print(f"Filtering barcodes by Reads/UMI ratio {opr.__name__} {threshold}")
    return df[opr(df.groupby("Barcode")["ReadCount"].transform('sum') / df.groupby("Barcode")["UMI"].transform('count'), threshold)]


pd.DataFrame.filter_ratio = filter_ratio


@_filter_wrapper
def filter_dups(df: pd.DataFrame, threshold: Union[float, int] = 2, min_len: int = 3) -> pd.DataFrame:
    """Remove barcodes that share a portion of their UMI-Targets combos based on the given threshold. If float
    then jaccard_index is used. If int then the there must be at least this many combos in common."""
    def nr_shared(set1, set2):
        return set1 & set2

    if isinstance(threshold, int):
        print(f"Filter barcodes who share >{threshold} UMI + Target combos ")
        compare_function = nr_shared
    else:
        print(f"Filter barcodes whose UMI + Target combos have a jaccard index >{threshold}")
        compare_function = jaccard_index

    d = df.copy()
    d["Target-UMI"] = d[["Target", "UMI"]].agg("-".join, axis=1)
    bcs = d.groupby("Barcode")["Target-UMI"].apply(set)
    bcs = dict(bcs[bcs.transform(len) >= min_len])

    total = (len(bcs) * (len(bcs) - 1)) / 2
    barcodes_to_remove = set()
    for bc1, bc2 in tqdm(combinations(bcs, 2), desc="Parsing pairs", total=total, unit_scale=True):
        if bc1 in barcodes_to_remove and bc2 in barcodes_to_remove:
            continue

        s1 = bcs[bc1]
        s2 = bcs[bc2]
        if compare_function(s1, s2) > threshold:
            barcodes_to_remove.update({bc1, bc2})

    return df[~df["Barcode"].isin(barcodes_to_remove)]


pd.DataFrame.filter_dups = filter_dups


@_filter_wrapper
def filter_targets(df: pd.DataFrame, threshold, opr=operator.gt) -> pd.DataFrame:
    print(f"Filtering barcodes by nr Targets {opr.__name__} {threshold}")
    return df[opr(df.groupby("Barcode")["Target"].transform('nunique'), threshold)]


pd.DataFrame.filter_targets = filter_targets


@_filter_wrapper
def filter_connected(df: pd.DataFrame, dist: int = 2) -> pd.DataFrame:
    """Remove barcodes that are within hamming dist of eachother leving isolated sequences
    Inspired by the Abseq analysis in https://www.nature.com/articles/srep44447
    """
    clusterer = UMIClusterer(cluster_method="cluster")
    # Encode each DBS for UMITools and perpare counts
    dbs_counts = {bytes(dbs, encoding='utf-8'): sum(group["ReadCount"]) for dbs, group in df.groupby("Barcode") \
                  if len(dbs) == 20}

    print(f"Pre cluster DBSs: {len(dbs_counts)}")

    # Cluster DBS
    clustered_dbss = clusterer(dbs_counts, threshold=dist)

    print(f"Clustered DBSs: {len(clustered_dbss)}")

    # Loop over clusters and save only isolated DBSs.
    dbs_kept = set()
    for cluster in clustered_dbss:
        seqs = [seq.decode("utf-8") for seq in cluster]

        if len(seqs) == 1:
            dbs_kept.add(seqs[0])
    return df[df["Barcode"].isin(dbs_kept)]


pd.DataFrame.filter_connected = filter_connected


@_filter_wrapper
def filter_shared_targets(df: pd.DataFrame, threshold, opr=operator.gt) -> pd.DataFrame:
    return df[opr(df.groupby(["UMI", "Target"])["Barcode"].transform('count'), threshold)]


pd.DataFrame.filter_shared_targets = filter_shared_targets


@_filter_wrapper
def filter_shared_umis(df: pd.DataFrame, threshold, opr=operator.gt) -> pd.DataFrame:
    return df[opr(df.groupby(["UMI"])["Barcode"].transform('count'), threshold)]


pd.DataFrame.filter_shared_umis = filter_shared_umis
