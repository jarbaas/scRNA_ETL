CREATE OR REPLACE VIEW vw_cluster_markers AS
WITH cluster_sizes AS (
    SELECT
        cell_type,
        COUNT(*)::real AS n_cells_cluster,
        SUM(COUNT(*)) OVER ()::real AS n_cells_total
    FROM dim_cells
    GROUP BY cell_type
),
gene_vs_rest AS (
    SELECT
        a.cell_type,
        a.gene_symbol,
        a.mean_log1p_expr AS mean_in,
        a.pct_expressing AS pct_in,
        GREATEST(
            (
                SUM(a.mean_log1p_expr * c.n_cells_cluster) OVER (PARTITION BY a.gene_symbol)
                - (a.mean_log1p_expr * c.n_cells_cluster)
            ) / (c.n_cells_total - c.n_cells_cluster),
            0.0
        ) AS mean_out,
        GREATEST(
            (
                SUM((a.pct_expressing / 100.0) * c.n_cells_cluster) OVER (PARTITION BY a.gene_symbol)
                - ((a.pct_expressing / 100.0) * c.n_cells_cluster)
            ) / (c.n_cells_total - c.n_cells_cluster) * 100.0,
            0.0
        ) AS pct_out
    FROM agg_cluster_expression a
    JOIN cluster_sizes c ON a.cell_type = c.cell_type
    WHERE a.gene_symbol !~ '^(RP[SL][0-9]|RPLP|MT-|MALAT1)'
),
scored_markers AS (
    SELECT
        *,
        LOG(2.0, ((mean_in + 0.1) / (mean_out + 0.1))::numeric) AS log2_fc
    FROM gene_vs_rest
    WHERE pct_in >= 25.0
),
ranked_markers AS (
    SELECT
        *,
        ROW_NUMBER() OVER (PARTITION BY cell_type ORDER BY log2_fc DESC) AS marker_rank
    FROM scored_markers
)
SELECT
    cell_type,
    marker_rank,
    gene_symbol,
    ROUND(mean_in::numeric, 2) AS mean_in,
    ROUND(mean_out::numeric, 2) AS mean_out,
    ROUND(pct_in::numeric, 1) AS pct_in,
    ROUND(pct_out::numeric, 1) AS pct_out,
    ROUND(log2_fc, 2) AS log2_fc
FROM ranked_markers;

-- Top 3 differential markers per cell type
SELECT *
FROM vw_cluster_markers
WHERE marker_rank <= 3
ORDER BY cell_type, marker_rank;

-- Single-gene expression joined across all QC-passing cells
SELECT
    c.cell_id,
    c.cell_type,
    c.umap_1,
    c.umap_2,
    COALESCE(f.log1p_norm, 0.0) AS expression
FROM dim_cells c
LEFT JOIN fact_expression f
    ON c.cell_id = f.cell_id
    AND f.gene_symbol = 'MS4A1';

-- Index scan verification
EXPLAIN ANALYZE
SELECT cell_id, log1p_norm
FROM fact_expression
WHERE gene_symbol = 'MS4A1';
