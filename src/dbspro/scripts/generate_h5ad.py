"""
Generate a h5ad file from the output of dbspro.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

DTYPES = {
    "Barcode": "object",
    "Target": "object",
    "UMI": "object",
    "ReadCount": int,
    "Sample": "category"
}


def main(input_data, output_data):
    data = pd.read_csv(input_data, sep="\t", dtype=DTYPES)

    counts = data.groupby(["Barcode", "Target"], as_index=False)["UMI"].count().set_index("Barcode")\
        .pivot(columns="Target", values="UMI").fillna(0)

    X = counts.values
    var = pd.DataFrame(index=counts.columns.values)
    obs = pd.DataFrame(index=counts.index)

    # Add metadata
    # Get the total number of reads for each barcode
    map_rc = data.groupby("Barcode")["ReadCount"].sum().to_dict()
    obs["total_reads"] = obs.index.to_series().map(map_rc)
    del map_rc

    # Get the total number of UMIs for each barcode
    obs["total_umis"] = X.sum(axis=1)

    # Get the ratio of reads to UMIs for each barcode
    obs["reads_per_umi"] = obs["total_reads"] / obs["total_umis"]

    # Get the number of non-zero count targets for each barcode
    obs["nr_targets"] = (X > 0).sum(axis=1)

    # Get the set of unique UMIs for each target and each barcode
    map_umis = data.groupby(["Barcode", "Target"], as_index=False)["UMI"]\
        .apply(set)\
        .set_index("Barcode")\
        .pivot(columns="Target", values="UMI")\
        .map(lambda x: {} if x is np.nan else x)\
        .to_dict(orient="index")

    map_umis = {k: ";".join(f"{t}:{','.join(v)}" for t, v in d.items() if v) for k, d in map_umis.items()}
    obs["target_umis"] = obs.index.to_series().map(map_umis)
    del map_umis

    # Get the sample
    obs["sample"] = data["Sample"].iloc[0]

    adata = sc.AnnData(X=csr_matrix(X), obs=obs, var=var)
    adata.var_names = counts.columns.values
    adata.obs_names = counts.index.values

    adata.write(output_data, compression="gzip")


if __name__ == "__main__":
    input_data = snakemake.input.data  # noqa: F821
    output_data = snakemake.output.data  # noqa: F821

    main(input_data, output_data)
