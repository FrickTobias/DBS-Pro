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
def filter_dups(df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    print(f"Filtering DBSs than share {threshold} or more UMI + Target combos")
    bcs = df.groupby(["Target", "UMI"])["Barcode"].apply(set)
    umis = list(bcs[bcs.transform(len) >= threshold])

    # Expected rate calculations
    # mu = 1 / len(umis)
    # pr_overlap = 1 - sum(poisson.pmf(k, mu) for k in range(threshold))
    # nr_expect_overlaps = pr_overlap * len(umis) * (len(umis) - 1) / 2
    # print(f"Expected overlaps {nr_expect_overlaps} (rate: {pr_overlap})")

    connected_barcodes = set()
    for c1, c2 in combinations(umis, 2):
        comb = c1 & c2
        if len(comb) >= threshold:
            connected_barcodes.update(comb)

    return df[~df["Barcode"].isin(connected_barcodes)].copy()


pd.DataFrame.filter_dups = filter_dups


@_filter_wrapper
def filter_dups_jaccard(df: pd.DataFrame, threshold: float = 0.8, min_len: int = 3) -> pd.DataFrame:
    print(f"Filter barcodes with at least {min_len} UMI + Target combos whose jaccard index > {threshold}")
    bcs = df.groupby(["Target", "UMI"])["Barcode"].apply(set)
    umis = list(bcs[bcs.transform(len) >= min_len])

    barcodes_to_remove = set()
    for c1, c2 in combinations(umis, 2):
        comb = c1 & c2
        if len(comb) > 0 and len(comb) / len(c1 | c2) > threshold:
            barcodes_to_remove.update(comb)

    return df[~df["Barcode"].isin(barcodes_to_remove)]


pd.DataFrame.filter_dups_jaccard = filter_dups_jaccard


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
