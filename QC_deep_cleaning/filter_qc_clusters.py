import marimo

__generated_with = "0.15.2"
app = marimo.App(
    width="full",
    app_title="QC Clustering",
    layout_file="layouts/filter_qc_clusters.grid.json",
)


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
    # Cell Filtering by QC Clustering Pipeline

    This notebook executes the quality control filtering pipeline for single-cell RNA-seq data. The pipeline:

    1. Clusters cells based on QC gene expression (mitochondrial, ribosomal, hemoglobin, transcription factor genes)

    2. Computes composite QC scores using PERMANOVA statistics

    3. Classifies clusters as "low_QC" or "acceptable"

    4. Saves filtered AnnData objects with low-quality cells removed


    **Author:** Raquel Cossio  
    **Date:** 2026-01-08
    """
    )
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import os
    from sklearn.metrics import pairwise_distances
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.distance import permanova
    from pathlib import Path
    return DistanceMatrix, Path, np, os, pairwise_distances, pd, permanova, sc


@app.cell
def _(Path, mo):
    def path_status(path: str):
        if not path:
            return mo.md("⚠️ No path provided")
        if Path(path).exists():
            return mo.md("✅ Path exists")
        return mo.md("❌ Path does not exist. Please change it for a valid path")
    return (path_status,)


@app.cell
def _(mo):
    # ============================================================================
    # INPUTS
    # ============================================================================
    denoised_input = mo.ui.text(value="/home/ratopin/big_data/cellbender/full_denosied",
                                label="PATH_TO_DENOISED_DATA",
                                full_width=True)
    mito_input = mo.ui.text(value= 'NMR-snRNA-seq/CrossTransfer/gene_lists/mito_genes.txt', 
                                label = 'PATH_TO_MITO_GENES',
                               full_width=True)
    tf_input = mo.ui.text(value = 'NMR-snRNA-seq/CrossTransfer//gene_lists/tf_genes.txt', 
                                label = 'PATH_TO_TF_GENES',
                               full_width=True)
    hemo_input = mo.ui.text(value = 'NMR-snRNA-seq/CrossTransfer//gene_lists/hemo_genes.txt', 
                                label = 'PATH_TO_HEMO_GENES',
                               full_width=True)
    output_dir = mo.ui.text(value = '/home/ratopin/big_data/filtered', 
                                label = 'OUTPUT_DIR',
                               full_width=True)
    return denoised_input, hemo_input, mito_input, output_dir, tf_input


@app.cell
def _(
    denoised_input,
    hemo_input,
    mito_input,
    mo,
    output_dir,
    path_status,
    tf_input,
):
    mo.vstack([mo.hstack([denoised_input, path_status(denoised_input.value)], widths= 'equal'),
               mo.hstack([mito_input, path_status(mito_input.value)], widths= 'equal'),
               mo.hstack([tf_input, path_status(tf_input.value)], widths= 'equal'), 
               mo.hstack([hemo_input, path_status(hemo_input.value)], widths= 'equal'),
               mo.hstack([output_dir, path_status(hemo_input.value)], widths= 'equal')], )
    return


@app.cell
def _(denoised_input, select_data_ui):
    ui, file_checkboxes = select_data_ui(denoised_input.value)
    return file_checkboxes, ui


@app.cell
def _(ui):
    ui
    return


@app.cell
def _(file_checkboxes):
    selected_input = [str(file) for use_cb, file in file_checkboxes if use_cb.value]
    return (selected_input,)


@app.cell
def _(selected_input):
    selected_input
    return


@app.cell
def _(denoised_input, hemo_input, mito_input, output_dir, tf_input):
    PATH_TO_DENOISED_DATA = denoised_input.value
    PATH_TO_MITO_GENES = mito_input.value
    PATH_TO_TF_GENES = tf_input.value
    PATH_TO_HEMO_GENES = hemo_input.value
    OUTPUT_DIR = output_dir.value
    return (
        PATH_TO_DENOISED_DATA,
        PATH_TO_HEMO_GENES,
        PATH_TO_MITO_GENES,
        PATH_TO_TF_GENES,
    )


@app.cell
def _(mo):
    # ============================================================================
    # PARAMETERS
    # ============================================================================
    leiden_res = mo.ui.slider(0.1,2.6, 0.05, label='Leiden Resolution', value = 0.2)
    leiden_iter = mo.ui.slider(-1,10, 1, label='Leiden Iterations', value = -1)
    pca_solver = mo.ui.multiselect(options = ['arpack','auto', 'full', 'randomized'], value = ['arpack'], label='PCA solver')
    umap_metric = mo.ui.multiselect(['euclidean','cosine'], value = ['cosine'],label='UMAP metric')
    permanova_permutations = mo.ui.slider(0,200, 10, label='Number of Permanova permutations', value = 20)
    return (
        leiden_iter,
        leiden_res,
        pca_solver,
        permanova_permutations,
        umap_metric,
    )


@app.cell
def _(mo):
    pseudo_F_weight = mo.ui.slider(-1,5, 0.2, value=2, label= 'Pseudo-F weight for qc_composite')
    return (pseudo_F_weight,)


@app.cell
def _(
    leiden_iter,
    leiden_res,
    mo,
    pca_solver,
    permanova_permutations,
    umap_metric,
):
    mo.vstack([leiden_res,
              leiden_iter,
              pca_solver, 
               umap_metric,
              permanova_permutations])
    return


@app.cell
def _(
    leiden_iter,
    leiden_res,
    pca_solver,
    permanova_permutations,
    umap_metric,
):
    # Heavy computation parameters
    LEIDEN_RESOLUTION = leiden_res.value
    LEIDEN_ITERATIONS = leiden_iter.value
    PCA_SOLVER = pca_solver.value[0]
    UMAP_METRIC = umap_metric.value[0]
    PERMANOVA_PERMUTATIONS = permanova_permutations.value
    return (
        LEIDEN_ITERATIONS,
        LEIDEN_RESOLUTION,
        PCA_SOLVER,
        PERMANOVA_PERMUTATIONS,
        UMAP_METRIC,
    )


@app.cell
def _(PCA_SOLVER, UMAP_METRIC):
    PCA_SOLVER
    UMAP_METRIC
    return


@app.cell
def _(mo):
    # Sample-specific threshold adjustments (std multiplier)
    SAMPLE_THRESHOLD_TOLERANCES = {
        'NMR1': mo.ui.slider(-0.5,1.5,.05,value=1.0, label='Threshold Tolerance for NMR1'),
        'NMR2': mo.ui.slider(-0.5,1.5,.05,value=1.0, label='Threshold Tolerance for NMR2'),
        'NMR3': mo.ui.slider(-0.5,1.5,.05,value=0.7, label='Threshold Tolerance for NMR3'),
        'NMR4': mo.ui.slider(-0.5,1.5,.05,value=0.7, label='Threshold Tolerance for NMR4'),
        'NMR5': mo.ui.slider(-0.5,1.5,.05,value=-0.15, label='Threshold Tolerance for NMR5'),
        'NMR6': mo.ui.slider(-0.5,1.5,.05,value=-0.15, label='Threshold Tolerance for NMR6'),
    }
    DEFAULT_THRESHOLD_TOLERANCE_COEF = mo.ui.slider(0.1,1.5,.1,value=1.0, label='Default Threshold Tolerance')
    SAMPLE_THRESHOLD_TOLERANCES
    return DEFAULT_THRESHOLD_TOLERANCE_COEF, SAMPLE_THRESHOLD_TOLERANCES


@app.cell
def _(pseudo_F_weight):
    pseudo_F_weight
    return


@app.cell
def _(pseudo_F_weight):
    # Light computation parameters

    PSEUDO_F_WEIGHT = pseudo_F_weight.value
    return (PSEUDO_F_WEIGHT,)


@app.cell
def _(mo, os):

    def select_data_ui(input_dir: str):
        """
        Opens a UI display to select which files from input_dir will be further processed.

        Parameters
        ----------
        input_dir : str
            Directory containing cellbender output H5 files

        Returns
        -------
        tuple of for UI and list of checkboxes

        """
        checkboxes = []

        for file in sorted(os.listdir(input_dir)):
            if not file.endswith(".h5"):
                continue

            # Expected: {SAMPLE}_{TISSUE}_denoised.h5
            parts = file.split('_',1)
            sample = parts[0]
            tissue = parts[1].split('_denoised.')[0]

            use_cb = mo.ui.checkbox(
                value=True,
                label=f"Use {sample} from {tissue}",
            )
            checkboxes.append((use_cb, file))

        return mo.vstack([use_cb for use_cb, _ in checkboxes]), checkboxes
    return (select_data_ui,)


@app.function
def load_gene_lists(path_mito, path_tf, path_hemo):
    """
    Load QC gene lists from text files.

    Parameters
    ----------
    path_mito : str
        Path to mitochondrial genes file (one gene per line)
    path_tf : str
        Path to transcription factor genes file
    path_hemo : str
        Path to hemoglobin genes file

    Returns
    -------
    tuple of lists
        (mito_genes, tf_genes, hemo_genes)
    """
    with open(path_mito, 'r', newline='\n') as f:
        mito_genes = f.read().split('\n')[:-1]
    with open(path_tf, 'r') as f:
        tf_genes = f.read().split('\n')[:-1]
    with open(path_hemo, 'r') as f:
        hemo_genes = f.read().split('\n')[:-1]

    return mito_genes, tf_genes, hemo_genes


@app.cell
def _(mo, run_chunk0, selected_input):
    run_result_0 = mo.ui.run_button(label="Load Data")

    def execute_chunk0():
        if run_result_0.value:
            print(f"Loading {selected_input}")
            return run_chunk0(selected_input)
        return None, None

    mo.vstack([mo.md("""
    ## Step 0: Load data
    - Load Data 
    - Calculate QC metrics
                     """), run_result_0])
    return execute_chunk0, run_result_0


@app.cell
def _(execute_chunk0, run_result_0):
    # Execute Chunk 1
    if run_result_0.value:
        adatas, sample_names = execute_chunk0()
    else:
        adatas, sample_names = None, None
    return adatas, sample_names


@app.cell
def _(
    PATH_TO_HEMO_GENES,
    PATH_TO_MITO_GENES,
    PATH_TO_TF_GENES,
    load_cellbender_data,
):
    def run_chunk0(selected_input):
        """
        Load data and compute QC-metrics


        Parameters
        ---------
        list of str:
            selected file names

        Returns
        ---------
        list of adatas
        """
        mito, tf, hemo = load_gene_lists(
            PATH_TO_MITO_GENES,
            PATH_TO_TF_GENES,
            PATH_TO_HEMO_GENES,
        )

        return load_cellbender_data(selected_input,mito, tf, hemo)
    return (run_chunk0,)


@app.cell
def _(adatas, mo, run_chunk1):
    if adatas is not None:
        run_result_1 = mo.ui.run_button(label="Run heavy computations") 
    else:
        run_result_1 = mo.md("*Load data first...*")
    def execute_chunk1(): 
        if run_result_1.value: 
            return run_chunk1(adatas) 
        return None, None, None 
    ui1 = mo.vstack([mo.md(""" 
    ## Step 1: Heavy Computations
    - Subset QC-genes
    - Normalize
    - Scale
    - PCA
    - Compute neighbors graph 
    - UMAP 
    - Leiden clustering """),
               run_result_1])
    ui1
    return execute_chunk1, run_result_1


@app.cell
def _(execute_chunk1, run_result_1):
    # Execute Chunk 1
    if run_result_1 is bool:
        if run_result_1.value:
            qc_adatas, pairwise_dfs = execute_chunk1()
        else:
            qc_adatas, pairwise_dfs = None, None
    return pairwise_dfs, qc_adatas


@app.cell
def _(cluster_qc_genes, compute_permanovas):
    def run_chunk1(adatas):
        qc_adatas = []
        pairwise_dfs = []

        for adata in adatas:
            qc_adata = cluster_qc_genes(adata)
            qc_adatas.append(qc_adata)
            pairwise_dfs.append(compute_permanovas(qc_adata))

        return qc_adatas, pairwise_dfs
    return (run_chunk1,)


@app.cell
def _(mo, qc_adatas, sc):
    if qc_adatas is None:
        mo.md("*Waiting for heavy computations...*")
        figs = None
    else:
        figs = []

        for qc_adata in qc_adatas:
            sample = qc_adata.obs.get('sample', 'unknown').iloc[0]

            fig = sc.pl.umap(qc_adata, color="leiden", show=False, return_fig=True)

            figs.append((sample, fig))
    return figs, qc_adata


@app.cell
def _(figs, mo):
    # Display figures
    def display_clusters(figs):
        if figs is None:
            return mo.md("*Waiting for Step 1 to finish*")
        else:
            return mo.vstack([
                mo.md(f"### {name} – UMAP by Leiden Algorithm"),
                fig,
            ] for name, fig in figs)
    display_clusters(figs)
    return


@app.cell
def _():
    return


@app.cell
def _(mo, pairwise_dfs, qc_adatas, run_chunck2):
    # UI button config

    if qc_adatas:
        run_result_2 = mo.ui.run_button(label="Run light computations (QC classification)")
    else:
        run_result_2 = mo.md("*Complete Step 1 first...*")

    def execute_chunck2(): 
        if run_result_2.value: 
            return run_chunck2(qc_adatas, pairwise_dfs) 
        return None, None
    ui2 = mo.vstack([mo.md("""
    ## Step 2: Light Computations
    - Calculate QC composite scores
    - Classify clusters as low_QC or acceptable
    - Add QC labels to cells
                     """), run_result_2])
    ui2
    return execute_chunck2, run_result_2


@app.cell
def _(execute_chunck2, run_result_2):
    if run_result_2.value:
        final_qc_adatas, summaries = execute_chunck2()
    else:
        final_qc_adatas, summaries = None, None
    return (final_qc_adatas,)


@app.cell
def _(process_light_computations):
    def run_chunck2(qc_adatas, pairwise_dfs):
        summaries = []
        final_qc_adatas = []

        for qc_adata, pairwise_df in zip(qc_adatas, pairwise_dfs):
            sample = qc_adata.obs.get('sample', 'unknown').iloc[0]
            final_qc_adata, summary = process_light_computations(qc_adata, pairwise_df, sample)
            summaries.append(summary)
            final_qc_adatas.append(final_qc_adata)
        return final_qc_adatas, summaries
    return (run_chunck2,)


@app.cell
def _(final_qc_adatas):
    final_qc_adatas
    return


@app.cell
def _(qc_adata, sc):
    def plot_qc_scores(final_qc_adatas):
        if final_qc_adatas is None:
            return None, None

        umap_figs_composite = []
        umap_figs_label = []

        for final_qc_adata in final_qc_adatas:
            sample = final_qc_adata.obs.get('sample', 'unknown').iloc[0]
            fig_composite = sc.pl.umap(qc_adata, color="qc_composite", show=False, return_fig=True, title=f"QC composite score - {sample}")
            umap_figs_composite.append(fig_composite)
            fig_label = sc.pl.umap(qc_adata, color="qc_label", show=False, return_fig=True, title=f"QC label - {sample}")
            umap_figs_label.append(fig_label)
        return umap_figs_composite, umap_figs_label

    return (plot_qc_scores,)


@app.cell
def _(final_qc_adatas, plot_qc_scores):
    if final_qc_adatas is not None:
        umap_figs_composite, umap_figs_label = plot_qc_scores(final_qc_adatas)
    else:
        umap_figs_composite, umap_figs_label = None, None
    return umap_figs_composite, umap_figs_label


@app.cell
def _(mo, sample_names, umap_figs_composite, umap_figs_label):
    # Display UMAP QC scores:
    def display_qc_scores(umap_figs_composite, umap_figs_label):
        if umap_figs_composite is None:
            return mo.md("*No UMAPs yet*")
        return mo.vstack([[
                mo.md(f"### {sample} – UMAP for QC Scores"),
                fig_composite, fig_label
        ] for fig_composite, fig_label, sample in zip(umap_figs_composite, umap_figs_label, sample_names)])
    display_qc_scores(umap_figs_composite, umap_figs_label)
    return


@app.cell
def _(final_qc_adatas, mo):
    if final_qc_adatas:
        run_result_3 = mo.ui.run_button(label="Save filtered data")
    else:
        run_result_3 = mo.md("*Complete Step 2 first...*")

    mo.vstack([mo.md("""
    ## Step 3: Save Results
    - Filter out low-QC cells
    - Save filtered AnnData objects
                     """), run_result_3])
    return (run_result_3,)


@app.cell
def _(
    adatas,
    final_qc_adatas,
    mo,
    run_chunk3,
    run_result_3,
    umap_figs_composite,
    umap_figs_label,
):
    if run_result_3.value if isinstance(run_result_3, object) and hasattr(run_result_3, 'value') else False:
        result_display = run_chunk3(adatas, final_qc_adatas, umap_figs_composite, umap_figs_label)
        mo.vstack([result_display])
    else:
        mo.md("*Waiting for Chunk 3 to complete...*")
    return


@app.cell
def _(PATH_TO_DENOISED_DATA, mo, os, sc):
    def load_cellbender_data(selected_input, mito_genes, tf_genes, hemo_genes):
        """
        Load and annotate cellbender-denoised data.

        Parameters
        ----------
        input_dir : str
            Directory containing cellbender output H5 files
        mito_genes : list
            Mitochondrial gene names
        tf_genes : list
            Transcription factor gene names
        hemo_genes : list
            Hemoglobin gene names

        Returns
        -------
        list of AnnData
            Annotated data objects with QC metrics calculated
        """
        adatas = []
        sample_names = []
        mo.md("Loading data...")

        for file in selected_input:
            # Parse filename: {SAMPLE}_{TISSUE}_denoised.h5
            parts = file.split('_')
            sample = parts[0]
            sample_names.append(sample)
            tissue = parts[1].split('_denoised.')[0]

            print(f'Loading {sample} {tissue}...')
            mo.md(f'**Loading {sample} {tissue}...**')

            # Read H5 file
            filepath = os.path.join(PATH_TO_DENOISED_DATA, file)
            adata = sc.read_10x_h5(filepath)

            # Add metadata
            adata.obs['species'] = 'nmr'
            adata.obs['sample'] = sample
            adata.obs['tissue'] = tissue
            adata.var['genome'] = 'HetGla1.0'

            # Flag QC genes
            adata.var['mt'] = [gene in mito_genes for gene in adata.var_names]
            adata.var['ribo'] = adata.var_names.str.startswith(("RPL", "RPS"))
            adata.var['hb'] = [gene in hemo_genes for gene in adata.var_names]
            adata.var['tf'] = [gene in tf_genes for gene in adata.var_names]

            # Calculate QC metrics
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", 'hb', 'tf'], inplace=True)
            sc.pp.filter_cells(adata, min_counts=1, inplace=True)

            print(f'  Total cells: {len(adata)}')
            mo.md(f'✅ Loaded: {len(adata)} cells')
            adatas.append(adata)

        return adatas, sample_names
    return (load_cellbender_data,)


@app.cell
def _(LEIDEN_ITERATIONS, LEIDEN_RESOLUTION, PCA_SOLVER, UMAP_METRIC, sc):
    def cluster_qc_genes(adata):
        """
        Cluster cells based on QC gene expression.

        Parameters
        ----------
        adata : AnnData
            Input data object
        leiden_resolution : float
            Leiden clustering resolution

        Returns
        -------
        AnnData
            QC-clustered data (subset to QC genes only)
        """
        qc_adata = adata.copy()

        # Subset to QC genes
        qc_adata.var['qc_gene'] = qc_adata.var[['mt', 'ribo', 'hb', 'tf']].any(axis=1)
        qc_genes = qc_adata.var_names[qc_adata.var['qc_gene']].to_list()

        if len(qc_genes) == 0:
            raise ValueError("No QC genes found; PCA would fail")

        qc_adata = qc_adata[:, qc_genes]

        # Standard preprocessing
        sc.pp.normalize_total(qc_adata)
        sc.pp.scale(qc_adata)
        sc.tl.pca(qc_adata, svd_solver=PCA_SOLVER, use_highly_variable=False)
        sc.pp.neighbors(qc_adata, metric=UMAP_METRIC)
        sc.tl.umap(qc_adata)
        sc.tl.leiden(qc_adata, resolution=LEIDEN_RESOLUTION, n_iterations=LEIDEN_ITERATIONS)

        return qc_adata
    return (cluster_qc_genes,)


@app.cell
def _(
    DistanceMatrix,
    PERMANOVA_PERMUTATIONS,
    mo,
    pairwise_distances,
    pd,
    permanova,
):
    def compute_permanovas(qc_adata, cluster_col = 'leiden', n_pc_used=5, metric = 'euclidean'):
        """
        Compute global PREMANOVA and pairwise PERMANOVA tests to asess cluster similarities between clusters. **PERMANOVA** tests whether centroids in multivariate space differ. The **pairwise PERMANOVA** extensions let you assess which clusters are most distant.

        Parameters
        ----------
        qc_adata : AnnData object
            Data with leiden clusters and UMAP coordinates. Depends on previosly running `sc.tl.neighbors()` and a clustering algorithm such as leiden

        Returns
        -------
        pd.DataFrame
            Parwise Dataframe of PERMANOVA results 
        """
        sample = qc_adata.obs['sample'].iloc[0]
        tissue = qc_adata.obs['tissue'].iloc[0]

        print(f'\nComputing QC scores for {sample} {tissue}...')

        mo.md(f'### Computing QC scores for {sample} {tissue}')
        locations = pd.DataFrame(
            qc_adata.obsm['X_pca'],
            index=qc_adata.obs_names,
            columns=['PC'+str(i+1) for i in range(n_pc_used)]
        )
        pairwise_dist = pairwise_distances(locations[['PC'+str(i+1) for i in range(n_pc_used)]], metric=metric)
        distance_matrix = DistanceMatrix(
            pairwise_dist,
            ids=qc_adata.obs_names,
            validate=True,
            condensed=False
        )
        res = permanova(distance_matrix, grouping=qc_adata.obs[cluster_col], 
                        permutations=PERMANOVA_PERMUTATIONS)
        print(f'Global PERMANOVA p-value: {res["p-value"]:.4f}')
        mo.md(f'✅ Global PERMANOVA p-value: {res["p-value"]:.4f}')

        # ========== Pairwise PERMANOVA ==========
        clusters = sorted(qc_adata.obs[cluster_col].unique())
        pairwise_results = []

        for g1 in clusters:
            for g2 in clusters:
                if g1 >= g2:
                    continue

                ids = qc_adata.obs_names[(qc_adata.obs['leiden'] == g1) | 
                                          (qc_adata.obs['leiden'] == g2)]
                sub_dm = distance_matrix.filter(ids, strict=True)
                sub_group = qc_adata.obs.loc[ids, 'leiden']

                result = permanova(sub_dm, grouping=sub_group, 
                                 permutations=PERMANOVA_PERMUTATIONS)

                pairwise_results.append({
                    'cluster1': g1,
                    'cluster2': g2,
                    'pseudo_F': result['test statistic'],
                    'p_value': result['p-value']
                })
                print(f"From cluster {g1} to {g2}: \tPseudo-F={result['test statistic']} \tp-value={result['p-value']}")
                mo.md(f"From cluster {g1} to {g2}: \tPseudo-F={result['test statistic']} \tp-value={result['p-value']}")

        return pd.DataFrame(pairwise_results)
    return (compute_permanovas,)


@app.cell
def _(PSEUDO_F_WEIGHT, pd):
    def calcule_qc_composite(qc_adata, pairwise_df):
        """
        Compute composite QC score using PERMANOVA and biological markers

        Parameters
        ----------
        qc_adata : AnnData object
            Data with leiden clusters and UMAP coordinates. Depends on previosly running `sc.tl.neighbors()` and a clustering algorithm such as leiden.
        parwise_df : pd.DataFrame
            Parwise Dataframe of PERMANOVA results  obtained with compute_permanovas(qc_adata)

        Returns
        -------
        pd.DataFrame
            Parwise Dataframe of PERMANOVA results 


        """
        # ========== Compute pseudo-F QC score per cluster ==========
        qc_scores = {}
        clusters = qc_adata.obs['leiden'].unique()
        for g in clusters:
            mask = (pairwise_df['cluster1'] == g) | (pairwise_df['cluster2'] == g)
            qc_scores[g] = pairwise_df.loc[mask, 'pseudo_F'].mean()

        qc_scores = pd.Series(qc_scores, name='qc_score')

        # ========== Combine metrics ==========
        qc_markers = qc_adata.obs[['total_counts', 'total_counts_tf',
                                   'pct_counts_mt', 'pct_counts_ribo',
                                   'pct_counts_hb', 'leiden']].copy()
        summary = qc_markers.groupby('leiden', observed=True).mean()
        summary['pseudoF_qc_score'] = qc_scores

        # Normalize each metric (higher = worse QC)
        summary['qc_composite'] = (
            -summary['total_counts'].rank() +                    # higer counts = better
            -summary['total_counts_tf'].rank() +                 # higher TF genes = better
            summary['pct_counts_mt'].rank() +                    # higher mito = worse
            summary['pct_counts_ribo'].rank() +                  # higher ribo = worse
            summary['pct_counts_hb'].rank() +                    # higher Hb = worse
            PSEUDO_F_WEIGHT * summary['pseudoF_qc_score'].rank() # pseudo-F weighted 2x
        )
        return summary
    return (calcule_qc_composite,)


@app.cell
def _(DEFAULT_THRESHOLD_TOLERANCE_COEF, SAMPLE_THRESHOLD_TOLERANCES, mo, np):
    def classify_and_label(qc_adata, summary, sample, treshold_tolerance_coef=None):
        """
        Classify clusters as "low_QC" or "acceptable" and add labels to adata.

        Parameters
        ----------
        qc_adata : AnnData
            Data with leiden clusters
        summary : pd.DataFrame
            QC summary with composite scores
        sample : str
            Sample name (for sample-specific thresholds)
        threshold_multiplier : float, optional
            Std multiplier for threshold. If None, uses sample-specific or default.

        Returns
        -------
        AnnData
            Data with qc_composite and qc_label added to .obs
        """
        if treshold_tolerance_coef is None:
            treshold_tolerance_coef = SAMPLE_THRESHOLD_TOLERANCES.get(sample, DEFAULT_THRESHOLD_TOLERANCE_COEF).value

        # Compute threshold
        threshold = summary['qc_composite'].mean() + treshold_tolerance_coef * summary['qc_composite'].std()

        # Classify
        summary['QC_label'] = np.where(summary['qc_composite'] >= threshold, 'low_QC', 'acceptable')

        # Add to adata
        qc_adata.obs['qc_composite'] = qc_adata.obs['leiden'].map(summary['qc_composite'].to_dict())
        qc_adata.obs['qc_label'] = qc_adata.obs['leiden'].map(summary['QC_label'].to_dict())

        print(f'  Threshold: {threshold:.2f} (multiplier={treshold_tolerance_coef})')
        print(f'  Low-QC clusters: {summary[summary["QC_label"]=="low_QC"].index.tolist()}')
        print(f'  Acceptable clusters: {summary[summary["QC_label"]=="acceptable"].index.tolist()}')
        mo.md(f'✅ Threshold: {threshold:.2f} | Low-QC: {summary[summary["QC_label"]=="low_QC"].index.tolist()} | Acceptable: {summary[summary["QC_label"]=="acceptable"].index.tolist()}')

        return qc_adata, summary
    return (classify_and_label,)


@app.cell
def _(mo, os):
    def save_filtered_data(qc_adata, adata, output_dir, save_filtered=True):
        """
        Filter out low-QC cells and save cleaned AnnData object.

        Parameters
        ----------
        qc_adata : AnnData
            Data with QC labels
        adata : AnnData
            Original full data object
        output_dir : str
            Directory to save output file

        Returns
        -------
        str
            Path to saved file
        file
            H5AD file containing all genes and cells with an added qc_composite and qc_label column.
        """
        sample = adata.obs['sample'].iloc[0]
        tissue = adata.obs['tissue'].iloc[0]

        # Transfer QC labels to original data
        adata.obs['qc_composite'] = qc_adata.obs['qc_composite']
        adata.obs['qc_label'] = qc_adata.obs['qc_label']

        # Filter
        adata_filtered = adata[adata.obs['qc_label'] == 'acceptable', :].copy()
        if save_filtered:
            # Save only acceptable cells with qc_scores
            output_path = os.path.join(output_dir, f'{sample}_{tissue}_filtered.h5ad')
            adata_filtered.write_h5ad(output_path)

        else: 
            # Save all cells with qc_scores
            output_path = os.path.join(output_dir, f'{sample}_{tissue}_unfiltered.h5ad')
            adata.write_h5ad(output_path)

        n_removed = len(adata) - len(adata_filtered)
        print(f'  Cells filtered: {n_removed} → remaining: {len(adata_filtered)}')
        print(f'  Saved to: {output_path}')
        mo.md(f'✅ **Saved:** {n_removed} cells removed | {len(adata_filtered)} cells remaining | {output_path}')    

        return output_path
    return


@app.cell
def _(cluster_qc_genes, compute_permanovas):
    def process_heavy_computations(adata):
        sample = adata.obs['sample'].iloc[0]
        tissue = adata.obs['tissue'].iloc[0]

        print(f'\n{"="*70}')
        print(f'Processing {sample} {tissue}')
        print(f'{"="*70}')

        # Cluster QC genes
        qc_adata = cluster_qc_genes(adata)
        print(f'  QC clusters: {len(qc_adata.obs["leiden"].unique())}')

        # Compute  global and pairwise PERMANOVAs
        pairwise_df = compute_permanovas(qc_adata)

        return qc_adata, pairwise_df
    return


@app.cell
def _(calcule_qc_composite, classify_and_label):
    def process_light_computations(qc_adata, pairwise_df, sample):
        # Compute qc_composite score
        summary = calcule_qc_composite(qc_adata, pairwise_df)

        # Classify and label
        final_qc_adata, summary = classify_and_label(qc_adata, summary, sample)

        return final_qc_adata, summary
    return (process_light_computations,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
