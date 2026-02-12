import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
import seaborn as sns
import os

def subsample_dataset(dataset_path: str, random_seed=71, fraction=0.05, min_per_cluster=34,
                      filter_by=None, add_info_to_stratify_by=None, class_col_name="cluster_id"):
    """
    Subsample an AnnData dataset in backed mode (don't modify backed object).
    Stratifies sampling by grouping cells and ensuring each group is represented.

    Parameters
    ----------
    dataset_path : str
        Path to the AnnData dataset in h5ad format (backed mode).
    random_seed : int, optional
        Random seed for reproducibility. Default is 71.
    fraction : float, optional
        Fraction of cells to sample from each stratification group. Default is 0.05 (5%).
    min_per_cluster : int, optional
        Minimum number of cells to sample from each stratification group. Default is 34.
    filter_by : dict, optional
        Dictionary of column names and regex patterns to filter cells before sampling. Default is None.
    add_info_to_stratify_by : tuple, optional
        Tuple of (column name, mapping dict) to add additional info for stratification. Default is None.
    class_col_name : str, optional
        Column name in obs that contains cluster labels (or the mapped stratification column after cluster_map).
        Default is "cluster_id".

    Returns
    -------
    AnnData
        Subsampled AnnData object.
    """
    dirname = os.path.dirname(dataset_path)
    basename = os.path.basename(dataset_path)
    subset_path = os.path.join(dirname, f"{basename.split(r'.h5')[0]}_subset{random_seed}.h5ad")

    if os.path.exists(subset_path):
        print("Loading existing subset...")
        human_data = sc.read(subset_path)
        print(f"Sampled dataset: {human_data.shape[0]} × {human_data.shape[1]}")
        return human_data

    print("Creating new subset...")
    adata_backed = sc.read_h5ad(dataset_path, backed="r")
    n_obs = adata_backed.n_obs

    # Build combined boolean mask (safe for multiple filters)
    if filter_by:
        mask = np.ones(n_obs, dtype=bool)
        for key, value in filter_by.items():
            if key not in adata_backed.obs.columns:
                raise NameError(f"Column '{key}' not found in obs")
            mask &= adata_backed.obs[key].str.contains(value, regex=True, na=False).values
    else:
        mask = np.ones(n_obs, dtype=bool)

    valid_positions = np.where(mask)[0]  # integer positions in the backed AnnData
    n_cells = len(valid_positions)
    target_size = int(n_cells * fraction)
    print(f"Original: {n_obs} cells, filtered: {n_cells} cells, target {target_size}")

    # Get cluster labels by positional indexing to avoid label/position ambiguity
    original_clusters = adata_backed.obs[class_col_name].iloc[valid_positions].reset_index(drop=True)

    # Map stratify values but fall back to original cluster if mapping misses
    if add_info_to_stratify_by:
        new_col, mapping_dict = add_info_to_stratify_by
        stratify_clusters = original_clusters.map(mapping_dict).fillna('Other').astype(str)
        print(f"Stratifying by mapped column with {stratify_clusters.nunique()} unique values")
    else:
        stratify_clusters = original_clusters.astype(str)

    np.random.seed(random_seed)
    sampled_indices = []
    summary = []

    # Use np.where on the values to guarantee integer positional arrays
    unique_groups = stratify_clusters.unique()
    for clust in unique_groups:
        pos_in_filtered = np.where(stratify_clusters.values == clust)[0]  # integer positions 0..len(valid_positions)-1
        if pos_in_filtered.size == 0:
            continue
        actual_indices = valid_positions[pos_in_filtered]  # convert to dataset row indices
        n_total = len(actual_indices)

        n_sample = int(n_total * fraction) + min_per_cluster
        n_sample = min(n_sample, n_total)

        chosen = np.random.choice(actual_indices, size=n_sample, replace=False)
        chosen_indexes = [int(x) for x in chosen]
        sampled_indices.extend(chosen_indexes)
        summary.append([clust, n_total, n_sample])

    sampled_indices = sorted(sampled_indices)
    if len(sampled_indices) == 0:
        raise ValueError("No cells sampled (check filters / fraction).")

    # Load only sampled cells into memory
    human_data = adata_backed[sampled_indices, :].to_memory()
    print(f"Final sampled dataset: {human_data.shape[0]} × {human_data.shape[1]}")

    df_summary = pd.DataFrame(summary, columns=["cluster", "total_cells", "sampled_cells"])
    df_summary["fraction_sampled"] = df_summary["sampled_cells"] / df_summary["total_cells"]
    print(df_summary.sort_values("total_cells").head(5))
    print("...")
    print(df_summary.sort_values("total_cells").tail(5))

    human_data.write_h5ad(subset_path)
    print(f"Sampled dataset saved to {subset_path}")
    return human_data

def plot_qc_violin_with_thres(adata, qc_cuts=None, region_col = 'region'):
    if qc_cuts is None:
        qc_cuts = {"n_genes_by_counts":3, 
        "min_total_counts":200, 
        "max_total_counts":22000,
        "pct_counts_mt": 3, 
        'pct_counts_ribo': 7}
    tissue = adata.obs[region_col].iloc[0]
    print(f"Plotting {tissue}")
    plot_keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt", 'pct_counts_ribo']
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for i, key in enumerate(plot_keys):
        # Create violin plot with seaborn
        data_to_plot = adata.obs[[key]].copy()
        sns.violinplot(data=data_to_plot, y=key, ax=axes[i], color='skyblue')
        
        # Add jitter points on top
        sns.stripplot(data=data_to_plot, y=key, ax=axes[i], color='black', alpha=0.4, size=1, jitter=1 )
        
        # Add threshold lines
        if key == "total_counts":
            axes[i].axhline(y=qc_cuts['min_total_counts'], color='red', linestyle='--', linewidth=2, alpha=1, label=f'min: {qc_cuts["min_total_counts"]}')
            axes[i].axhline(y=qc_cuts['max_total_counts'], color='red', linestyle='--', linewidth=2, alpha=0.7, label=f'max: {qc_cuts["max_total_counts"]}')
            axes[i].legend()
        elif key in qc_cuts:
            threshold = qc_cuts[key]
            axes[i].axhline(y=threshold, color='red', linestyle='--', linewidth=2, alpha=0.8, label=f'threshold: {threshold}')
            axes[i].legend()
        
        axes[i].set_title(f'{key}', fontsize=12)
    
    fig.suptitle(f'{tissue}', fontsize=14, fontweight='bold')
    return fig, axes
    
def umap_alpha(adata, color_col='leiden', alpha_col ='total_counts', size = 5, normalize_a= True, embedding = "X_umap"):
    emb = adata.obsm[embedding]
    clusters = adata.obs[color_col]      # or any cluster column
    values = adata.obs[alpha_col]  # column that drives alpha
        # normalize alpha into [0.05, 1]
    alpha = values.values.astype(float)

    if normalize_a:
        alpha = (alpha - alpha.min()) / (alpha.max() - alpha.min() + 1e-9)
        alpha = 0.05 + 0.95 * alpha

    # assign colors per cluster (Scanpy already stores them after sc.pl.umap)
    palette = adata.uns["leiden_colors"]
    # Convert cluster categories to integers for sorting
    cluster_labels = list(clusters.unique())
    cluster_labels_int = sorted(cluster_labels, key=lambda x: int(x))

        # Reorder colors according to numeric order
    color_map = dict(zip(cluster_labels, palette))
    ordered_colors = [color_map[k] for k in cluster_labels_int]

    fig, ax = plt.subplots(figsize=(7, 6))

    # plot scatter in the same numeric order
    for clust in cluster_labels_int:
        color = color_map[clust]
        idx = clusters == clust
        ax.scatter(
            emb[idx, 0],
            emb[idx, 1],
            s=size,
            c=color,
            alpha=alpha[idx],
            linewidth=0,
        )

    # build proxy legend
    handles = [
        Line2D([0], [0], marker="o", color=color, markersize=6, linestyle="",
            markerfacecolor=color, markeredgewidth=0)
        for color in ordered_colors
    ]

    ax.legend(
        handles,
        cluster_labels_int,
        title="Clusters",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"UMAP with per-cell alpha (driven by {alpha_col})")

    plt.tight_layout()
    plt.show()


def describe_obs_by_cluster(adata, cluster_key, columns):
    """
    Compute mean, min, max, std for selected obs columns grouped by clusters.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    cluster_key : str
        Column in adata.obs that defines clusters (e.g. "leiden").
    columns : list of str
        obs columns to summarize.

    Returns
    -------
    pandas.DataFrame
        Multi-index columns (metric, variable) and index = clusters.
    """
    df = adata.obs[[cluster_key] + columns].copy()
    df[cluster_key] = df[cluster_key].astype(str)

    grouped = df.groupby(cluster_key)   # pandas groupby

    stats = grouped.agg(["mean", "min", "max", "std"])
    return stats

def boxplot_obs_by_cluster(adata, cluster_key, columns):
    """
    Plot boxplots for selected obs columns grouped by cluster.
    One figure per obs column.
    """
    df = adata.obs[[cluster_key] + columns].copy()
    df[cluster_key] = df[cluster_key].astype(str)

    # sort cluster labels numerically
    cluster_order = sorted(df[cluster_key].unique(), key=lambda x: int(x))

    for col in columns:
        plt.figure(figsize=(10, 5))
        plt.title(f"Distribution of {col} by cluster")

        data = [df.loc[df[cluster_key] == cl, col].dropna().values
                for cl in cluster_order]

        plt.boxplot(
            data,
            labels=cluster_order,
            notch=False,
            patch_artist=True
        )

        plt.xlabel("Cluster")
        plt.ylabel(col)
        plt.grid(axis="y", linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.show()


def check_pc_correlation(adata, n_pcs=None, pc_col_name="X_pca"):
    """Check correlation between PCs and QC metrics to identify potential confounders."""
    from scipy.stats import spearmanr, pearsonr
    from statsmodels.stats.multitest import multipletests
    import matplotlib.pyplot as plt

    if n_pcs is None:
        n_pcs_check = adata.obsm[pc_col_name].shape[1]
    else: n_pcs_check = min(n_pcs, adata.obsm[pc_col_name].shape[1])
    pcs = pd.DataFrame(
        adata.obsm[pc_col_name][:, :n_pcs_check],
        index=adata.obs_names,
        columns=["PC"+str(i+1).zfill(2) for i in range(n_pcs_check)]
    )

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", 'pct_counts_ribo']
    qc = adata.obs[qc_cols].astype(float)
    qc['species'] = adata.obs['species'].astype('category').cat.codes
    qc_cols.append('species')

    rows = []
    for pc in pcs.columns:
        for cov in qc_cols:
            # Spearman es más robusto; cambia a pearsonr si prefieres lineal
            r, p = spearmanr(pcs[pc].values, qc[cov].values, nan_policy="omit")
            rows.append({"PC": pc, "covariate": cov, "r": r, "p": p})

    corr = pd.DataFrame(rows)
    corr["q"] = multipletests(corr["p"], method="fdr_bh")[1]
    return corr

def plot_pc_qc_correlation(corr, save_path=None):
    """Plot heatmap of PC vs QC metric correlations."""
    mat = corr.pivot(index="covariate", columns="PC", values="r")
    plt.figure(figsize=(10,3))
    im = plt.imshow(mat.values, aspect="auto", interpolation="nearest")
    plt.yticks(range(mat.shape[0]), mat.index)
    plt.xticks(range(mat.shape[1]), mat.columns, rotation=90)
    plt.colorbar(im, label="Spearman r")
    plt.title("Correlation PC vs QC")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()
    