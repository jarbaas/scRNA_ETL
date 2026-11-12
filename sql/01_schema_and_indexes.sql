CREATE TABLE IF NOT EXISTS dim_cells (
    cell_id VARCHAR(64) PRIMARY KEY,
    sample_id VARCHAR(64),
    cell_type VARCHAR(100),
    umap_1 REAL,
    umap_2 REAL,
    total_counts REAL,
    n_genes INTEGER,
    pct_counts_mt REAL
);

CREATE TABLE IF NOT EXISTS fact_expression (
    cell_id VARCHAR(64) NOT NULL,
    gene_symbol VARCHAR(64) NOT NULL,
    log1p_norm REAL NOT NULL,
    CONSTRAINT fk_fact_cell FOREIGN KEY (cell_id) REFERENCES dim_cells(cell_id)
);

CREATE TABLE IF NOT EXISTS agg_cluster_expression (
    cell_type VARCHAR(100) NOT NULL,
    gene_symbol VARCHAR(64) NOT NULL,
    mean_log1p_expr REAL NOT NULL,
    pct_expressing REAL NOT NULL,
    PRIMARY KEY (cell_type, gene_symbol)
);

CREATE INDEX IF NOT EXISTS idx_fact_gene_symbol ON fact_expression(gene_symbol);
