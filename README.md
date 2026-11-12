# scRNA_ETL & PostgreSQL Warehouse

ETL pipeline for processing sparse matrix single-cell RNA-sequencing datasets using Apache Parquet, PySpark on Databricks Serverless, Unity Catalog Delta Lake, and PostgreSQL.

## Architecture

1. **Extraction (`scripts/export_h5ad_to_parquet.py`):** Converts `.h5ad` sparse count matrices and cell metadata into Snappy-compressed Parquet files storing non-zero coordinate records.
2. **Distributed Transformation (`notebooks/01_scrna_pyspark_etl.py`):** Computes per-cell QC metrics (UMI depth, detected features, mitochondrial read fraction) in PySpark, filters low-quality droplets, performs $\ln(\text{CP10K} + 1)$ normalization, and writes intermediate tables to Unity Catalog Delta Lake.
3. **Relational Serving (`sql/01_schema_and_indexes.sql`):** Broadcast-filters top expressed features and batch-loads normalized tables into PostgreSQL using Databricks Secret Scopes and `cascadeTruncate` idempotency.
4. **Analytical Views (`sql/02_views_and_queries.sql`):** Implements cell-weighted one-vs-rest $\log_2$ fold-change marker ranking via SQL window functions alongside B-tree indexed single-gene lookups.

## Usage

### 1. Generate Staging Parquet Files

```bash
pip install -r requirements.txt

# Default PBMC3k dataset (2,700 cells; ~2.28M non-zero records)
python scripts/export_h5ad_to_parquet.py --output-dir ./parquet_staging

# Custom .h5ad input
python scripts/export_h5ad_to_parquet.py \
    --input-h5ad path/to/dataset.h5ad \
    --sample-col sample_id \
    --cell-type-col cell_type \
    --output-dir ./parquet_staging
```

Upload `cell_metadata.parquet` and `raw_counts_long.parquet` to `/Volumes/workspace/default/scrna_raw`.

### 2. Initialize PostgreSQL Schema

Run `sql/01_schema_and_indexes.sql` followed by `sql/02_views_and_queries.sql` in PostgreSQL to create the relational schema, `idx_fact_gene_symbol` index, and `vw_cluster_markers` view.

### 3. Configure Secrets & Run Databricks Workflow

Store PostgreSQL credentials (`host` and `password`) in a Databricks Secret Scope named `neon_postgres`, import `notebooks/01_scrna_pyspark_etl.py` into your Databricks Workspace, and schedule/trigger the job defined in `workflows/databricks_job.yaml`.

## Example Output (`vw_cluster_markers`, PBMC3k)

| cell_type | marker_rank | gene_symbol | mean_in | mean_out | pct_in | pct_out | log2_fc |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: |
| B cells | 1 | CD79A | 2.79 | 0.07 | 93.6 | 4.2 | 4.09 |
| B cells | 2 | TCL1A | 1.76 | 0.04 | 62.3 | 2.2 | 3.74 |
| B cells | 3 | MS4A1 | 2.15 | 0.08 | 86.0 | 5.3 | 3.61 |
| CD14+ Monocytes | 1 | S100A8 | 3.89 | 0.22 | 97.5 | 12.1 | 3.63 |
| CD14+ Monocytes | 2 | LGALS2 | 2.52 | 0.11 | 89.8 | 6.2 | 3.61 |
| CD14+ Monocytes | 3 | CD14 | 1.51 | 0.04 | 66.3 | 2.9 | 3.48 |
| CD4 T cells | 1 | MAL | 0.50 | 0.04 | 27.2 | 2.8 | 2.06 |
| CD4 T cells | 2 | LEF1 | 0.58 | 0.08 | 31.5 | 5.1 | 1.89 |
| CD4 T cells | 3 | CCR7 | 0.69 | 0.12 | 36.5 | 6.8 | 1.83 |
| CD8 T cells | 1 | GZMK | 1.53 | 0.11 | 55.7 | 5.3 | 2.95 |
| CD8 T cells | 2 | CCL5 | 3.45 | 0.58 | 93.0 | 23.3 | 2.38 |
| CD8 T cells | 3 | CST7 | 1.92 | 0.29 | 75.9 | 12.0 | 2.38 |
| Dendritic cells | 1 | FCER1A | 2.15 | 0.02 | 83.8 | 0.9 | 4.27 |
| Dendritic cells | 2 | CLEC10A | 1.36 | 0.02 | 70.3 | 1.4 | 3.57 |
| Dendritic cells | 3 | CD1C | 0.89 | 0.03 | 51.4 | 2.0 | 2.88 |
| FCGR3A+ Monocytes | 1 | MS4A7 | 1.76 | 0.14 | 80.0 | 7.7 | 2.98 |
| FCGR3A+ Monocytes | 2 | HES4 | 0.96 | 0.04 | 58.7 | 2.3 | 2.96 |
| FCGR3A+ Monocytes | 3 | RP11-290F20.3 | 1.62 | 0.13 | 81.3 | 7.4 | 2.90 |
| Megakaryocytes | 1 | PF4 | 4.80 | 0.02 | 100.0 | 1.0 | 5.39 |
| Megakaryocytes | 2 | PPBP | 5.76 | 0.04 | 100.0 | 2.3 | 5.39 |
| Megakaryocytes | 3 | GNG11 | 3.83 | 0.01 | 86.7 | 1.0 | 5.10 |
| NK cells | 1 | GZMB | 3.60 | 0.15 | 96.1 | 6.8 | 3.91 |
| NK cells | 2 | FGFBP2 | 2.81 | 0.13 | 89.6 | 6.1 | 3.68 |
| NK cells | 3 | SPON2 | 1.89 | 0.07 | 73.4 | 3.9 | 3.58 |
