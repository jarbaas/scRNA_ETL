"""Export AnnData (.h5ad) expression matrices and cell metadata to Apache Parquet."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def extract_nonzero_counts(adata: ad.AnnData, layer: str | None = None) -> pd.DataFrame:
    matrix = adata.layers[layer] if layer and layer in adata.layers else adata.X
    coo = sp.coo_matrix(matrix)
    return pd.DataFrame(
        {
            "cell_id": adata.obs_names.to_numpy()[coo.row].astype(str),
            "gene_symbol": adata.var_names.to_numpy()[coo.col].astype(str),
            "raw_count": coo.data.astype("float32"),
        }
    )


def load_pbmc3k() -> tuple[pd.DataFrame, pd.DataFrame]:
    adata_raw = sc.datasets.pbmc3k()
    adata_raw.var_names_make_unique()
    adata_proc = sc.datasets.pbmc3k_processed()

    obs_df = pd.DataFrame({"cell_id": adata_raw.obs_names.astype(str)}).set_index("cell_id")
    obs_df["sample_id"] = "pbmc3k_donor1"

    annotations = pd.DataFrame(
        {
            "cell_type": adata_proc.obs["louvain"].astype(str).values,
            "umap_1": adata_proc.obsm["X_umap"][:, 0].astype("float32"),
            "umap_2": adata_proc.obsm["X_umap"][:, 1].astype("float32"),
        },
        index=adata_proc.obs_names.astype(str),
    )

    obs_df = obs_df.join(annotations, how="left")
    obs_df["cell_type"] = obs_df["cell_type"].fillna("LowQC_Unassigned")
    obs_df["umap_1"] = obs_df["umap_1"].fillna(0.0).astype("float32")
    obs_df["umap_2"] = obs_df["umap_2"].fillna(0.0).astype("float32")

    return obs_df.reset_index(), extract_nonzero_counts(adata_raw)


def load_h5ad(
    h5ad_path: Path,
    sample_col: str,
    cell_type_col: str,
    counts_layer: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names_make_unique()

    umap = adata.obsm["X_umap"] if "X_umap" in adata.obsm else None
    obs_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "sample_id": (
                adata.obs[sample_col].astype(str).values
                if sample_col in adata.obs
                else h5ad_path.stem
            ),
            "cell_type": (
                adata.obs[cell_type_col].astype(str).values
                if cell_type_col in adata.obs
                else "Unassigned"
            ),
            "umap_1": umap[:, 0].astype("float32") if umap is not None else 0.0,
            "umap_2": umap[:, 1].astype("float32") if umap is not None else 0.0,
        }
    )

    return obs_df, extract_nonzero_counts(adata, layer=counts_layer)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("./parquet_staging"))
    parser.add_argument("--sample-col", type=str, default="sample_id")
    parser.add_argument("--cell-type-col", type=str, default="cell_type")
    parser.add_argument("--counts-layer", type=str, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.input_h5ad is None:
        logger.info("Loading Scanpy PBMC3k dataset")
        obs_df, counts_df = load_pbmc3k()
    else:
        logger.info("Reading %s", args.input_h5ad)
        obs_df, counts_df = load_h5ad(
            args.input_h5ad,
            sample_col=args.sample_col,
            cell_type_col=args.cell_type_col,
            counts_layer=args.counts_layer,
        )

    meta_path = args.output_dir / "cell_metadata.parquet"
    counts_path = args.output_dir / "raw_counts_long.parquet"

    obs_df.to_parquet(meta_path, index=False, compression="snappy")
    counts_df.to_parquet(counts_path, index=False, compression="snappy")

    logger.info(
        "Exported %s cells and %s non-zero count records to %s",
        f"{len(obs_df):,}",
        f"{len(counts_df):,}",
        args.output_dir,
    )


if __name__ == "__main__":
    main()
