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
from umi_tools import UMIClusterer
from typing import Union
from scipy.stats.mstats import gmean
from scipy.sparse import csr_matrix


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
    return df[opr(df.groupby("Barcode")["ReadCount"].transform('sum') /
                  df.groupby("Barcode")["UMI"].transform('count'), threshold)]


pd.DataFrame.filter_ratio = filter_ratio


@_filter_wrapper
def filter_target_count(df: pd.DataFrame, threshold: int, opr=operator.gt) -> pd.DataFrame:
    """Filter UMI count per target"""
    print(f"Filtering for targets with UMI count {opr.__name__} {threshold}")
    return df[opr(df.groupby(["Barcode", "Target"])["UMI"].transform('count'), threshold)]


pd.DataFrame.filter_target_count = filter_target_count


@_filter_wrapper
def filter_dups(df: pd.DataFrame, threshold: Union[float, int] = 2, min_len: int = 3) -> pd.DataFrame:
    """Remove barcodes that share a portion of their UMI-Targets combos based on the given threshold. If float
    then jaccard_index is used. If int then the there must be at least this many combos in common."""
    def nr_shared(set1, set2):
        return set1 & set2

    d = df.copy()
    d["Target-UMI"] = d[["Target", "UMI"]].agg("-".join, axis=1)
    bcs = d.groupby("Barcode")["Target-UMI"].apply(set)
    bcs = dict(bcs[bcs.transform(len) >= min_len])

    indptr = [0]
    indices = []
    data = []
    umis = {}
    barcodes = []

    for bc, umi_set in bcs.items():
        for umi in umi_set:
            index = umis.setdefault(umi, len(umis))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))
        barcodes.append(bc)
    barcodes = np.array(barcodes)

    arr = csr_matrix((data, indices, indptr), dtype=int)
    overlapps = (arr * arr.transpose())
    bcs_rows, bcs_cols = overlapps.nonzero()
    in_lower = bcs_rows < bcs_cols
    bcs_rows = bcs_rows[in_lower]
    bcs_cols = bcs_cols[in_lower]
    nr_overlapps = np.array(overlapps[bcs_rows, bcs_cols]).flatten()

    barcodes_to_remove = set()
    if isinstance(threshold, int):
        print(f"Filter barcodes who share >{threshold} UMI + Target combos ")
        dups = nr_overlapps > threshold
        dup_cols = bcs_cols[dups]
        dup_rows = bcs_rows[dups]
        barcodes_to_remove |= set(barcodes[dup_cols])
        barcodes_to_remove |= set(barcodes[dup_rows])
    else:
        print(f"Filter barcodes whose UMI + Target combos have a jaccard index >{threshold}")
        bcs_col_counts = np.array(overlapps[bcs_cols, bcs_cols]).flatten()
        bcs_rows_counts = np.array(overlapps[bcs_rows, bcs_rows]).flatten()
        totals = bcs_col_counts + bcs_rows_counts - nr_overlapps
        jaccard_values = nr_overlapps / totals
        dups = jaccard_values > threshold

        dup_cols = bcs_cols[dups]
        dup_rows = bcs_rows[dups]
        barcodes_to_remove |= set(barcodes[dup_cols])
        barcodes_to_remove |= set(barcodes[dup_rows])

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
    dbs_counts = {bytes(dbs, encoding='utf-8'): sum(group["ReadCount"]) for dbs, group in df.groupby("Barcode")
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

###################
# TRANSFORMATIONS #
###################


def to_matrix(df: pd.DataFrame, norm: bool = False, on_cells: bool = False, qc: bool = False) -> pd.DataFrame:
    """Convert from long (raw) fromat to wide format. Each row is a barcode with the corresponding UMI count for
    each target.
    """
    samples = sorted(list(set(df["Sample"])))
    targets = sorted(list(set(df["Target"])))
    matrix = []
    for sample in samples:
        df_sample = df[df["Sample"] == sample].copy()
        matrix_sample = df_sample.groupby(["Barcode", "Target"], as_index=False)["UMI"].count().set_index("Barcode")\
            .pivot(columns="Target", values="UMI").fillna(0)

        matrix_sample["Sample"] = sample
        matrix_sample["Sample"] = matrix_sample["Sample"].astype("category")

        if norm:
            matrix_sample = clr_normalize(matrix_sample, on_cells=on_cells)

        if qc:
            matrix_sample._qc(inplace=True)
            map_rc = {bc: group["ReadCount"].sum() for bc, group in df_sample.groupby("Barcode")}
            matrix_sample["total_reads"] = matrix_sample.index.to_series().map(map_rc)

        matrix.append(matrix_sample)

    matrix = pd.concat(matrix)
    # Target_columns come first and replace NaN
    matrix = matrix[[*targets, *[c for c in matrix if c not in targets]]]
    matrix.loc[:, targets] = matrix.loc[:, targets].fillna(0)
    return matrix


pd.DataFrame.to_matrix = to_matrix


def to_readmatrix(df: pd.DataFrame, norm: bool = False, on_cells: bool = False, qc: bool = False) -> pd.DataFrame:
    """Convert from long (raw) fromat to wide format. Each row is a barcode with the corresponding read count for
    each target.
    """
    samples = sorted(list(set(df["Sample"])))
    targets = sorted(list(set(df["Target"])))
    matrix = []
    for sample in samples:
        df_sample = df[df["Sample"] == sample].copy()
        matrix_sample = df_sample.groupby(["Barcode", "Target"], as_index=False)["ReadCount"].sum()\
            .set_index("Barcode").pivot(columns="Target", values="ReadCount").fillna(0)

        matrix_sample["Sample"] = sample
        matrix_sample["Sample"] = matrix_sample["Sample"].astype("category")

        if norm:
            matrix_sample = clr_normalize(matrix_sample, on_cells=on_cells)

        if qc:
            matrix_sample._qc(inplace=True)
            map_rc = {bc: group["ReadCount"].sum() for bc, group in df_sample.groupby("Barcode")}
            matrix_sample["total_reads"] = matrix_sample.index.to_series().map(map_rc)

        matrix.append(matrix_sample)

    matrix = pd.concat(matrix)
    # Target_columns come first and replace NaN
    matrix = matrix[[*targets, *[c for c in matrix if c not in targets]]]
    matrix.loc[:, targets] = matrix.loc[:, targets].fillna(0)
    return matrix


pd.DataFrame.to_readmatrix = to_readmatrix


def _qc(df, targets=None, inplace=False):
    if inplace:
        matrix = df
    else:
        matrix = df.copy()

    if targets is None:
        targets = matrix.select_dtypes(include=np.number).columns.tolist()

    matrix["total_count"] = matrix.loc[:, targets].sum(axis=1)
    matrix["total_count"] = matrix["total_count"].astype(np.int64)
    matrix["nr_targets"] = matrix.loc[:, targets].gt(0).sum(axis=1)
    return matrix


pd.DataFrame._qc = _qc

##################
# NORMALIZATIONS #
##################


def clr_normalize(matrix: pd.DataFrame, columns=None, on_cells: bool = False) -> pd.DataFrame:
    """CLR normalisation"""
    matrix_copy = matrix.copy()

    if not columns:
        columns = matrix_copy.select_dtypes(include=np.number).columns.tolist()

    sub = matrix_copy.loc[:, columns]
    if on_cells:
        # Based on seurat CLR function
        # https://github.com/satijalab/seurat/blob/9843b843ed0c3429d86203011fda00badeb29c2e/R/preprocessing.R#L2192
        sub = sub.apply(lambda x: np.log1p(x / gmean(x[x > 0])))
    else:
        sub = sub.apply(lambda x: np.log1p(x / gmean(x[x > 0])), axis=1)
    matrix_copy.loc[:, columns] = sub
    return matrix_copy


pd.DataFrame.clr_normalize = clr_normalize
