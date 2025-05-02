CMVerify: CMV Serostatus Predictor
===================================

CMVerify is a tool designed to predict Cytomegalovirus (CMV) serostatus based on single-cell RNA sequencing data. It uses a machine learning model to analyze blood transcriptomic data and provide predictions on CMV status. The tool is built to be easy to use, with automatic data preprocessing, model annotation, and prediction generation.
- [Ucar Lab at The Jackson Laboratory](https://ucarlab.github.io/)

Installation
------------

To install CMVerify, you can use `pip`:

    pip install cmverify

Alternatively, if you'd like to install the package from source:

    git clone https://github.com/luketrinityjax/CMVerify_v0.git
    cd CMVerify
    pip install .

Requirements
------------

- Python 3.8+
- `scanpy`
- `matplotlib`
- `numpy`
- `pandas`
- `seaborn`
- `anndata`
- `celltypist`
- `ipython`
- `scikit-learn`
- `joblib`

These dependencies are automatically installed when you install CMVerify.
[For R users please refer to this vignette for converting from Seurat to Anndata](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

Usage
-----

Here is how you can use CMVerify in your Python environment:

### 1. Importing CMVerify

To use CMVerify, start by importing the necessary module in your Python script:

    from cmverify import predict, visualize

### 2. Data Preparation

Before running predictions, ensure your data is prepared correctly. CMVerify requires **raw counts** for accurate predictions.

#### If raw counts are in `.raw`, convert `.raw` to `.X`:

You can convert your AnnData object from `.raw` to `.X` as follows:

    import scanpy as sc

    # Read the .h5ad file containing the raw data
    adata_pre = sc.read_h5ad('path_to_data.h5ad')

    # Use the raw counts for analysis (and retain metadata in .obs and .var)
    adata = sc.AnnData(X=adata_pre.raw.X, obs=adata_pre.obs, var=adata_pre.raw.var)

This ensures that the data used for predictions is based on raw counts, which is essential for the analysis.

### 3. Making Predictions

To make predictions, you need to load your single-cell RNA-seq data (`AnnData` object) and provide the relevant parameters. The `predict` function handles data normalization, model annotation, and CMV status prediction.

#### Example:

    import scanpy as sc
    from cmverify import predict

    # Load your single-cell RNA-seq data (AnnData object)
    adata = sc.read('path_to_data.h5ad')

    # Specify the columns for donor and longitudinal data
    donor_obs_column = 'donor_id'
    longitudinal_obs_column = 'timepoint'

    # Predict CMV status
    results = predict(adata, donor_obs_column, longitudinal_obs_column)

    # Output the predictions
    print(results)

### 4. Return Fractions

You can also return the calculated cell type fractions along with the predictions by setting the `return_frac` parameter to `True`.

    results, fractions_df = predict(adata, donor_obs_column, longitudinal_obs_column, return_frac=True)

    # Display the first 5 rows of the cell type fractions
    print(fractions_df.head())

## 4b. Append Ground Truth CMV Status (Optional)

If you have true CMV status in a separate metadata file (not in adata.obs), you can use the append_status function to match and append it to the prediction output.

    from cmverify import append_status
    
    # Add true labels to predictions
    updated_predictions = append_status(results, cmv_metadata_df, patient_col='patientID', cmv_col='CMV')


### 5. Visualize Longitudinal Predictions

You can visualize longitudinal CMV prediction probabilities across timepoints using the `visualize` function.

```python
from cmverify import visualize

visualize(results, visit_order=['Baseline', 'Month3', 'Month6'], save=True, filename='example_filename.png')
```

This will generate a figure connecting donor predictions across visits and mark the decision threshold.

Functions
---------

### `predict(adata, donor_obs_column, longitudinal_obs_column=None, verbose=1, return_frac=False, true_status=None, norm=True)`

This function predicts CMV status using pre-trained models on single-cell RNA-seq data. It normalizes, transforms, annotates, and calculates cell type fractions per donor before predicting CMV status.

**Parameters**:
- `adata`: An AnnData object containing your single-cell RNA-seq data.
- `donor_obs_column`: The column in `adata.obs` that contains the donor ID.
- `longitudinal_obs_column` (optional): The column in `adata.obs` for timepoints (e.g., for longitudinal data).
- `verbose` (optional): Verbosity level for progress messages. Default is 1 (show output), set to 0 to suppress.
- `return_frac` (optional): If `True`, returns the fractions DataFrame along with predictions.
- `true_status` (optional): The column in `adata.obs` for true donor serostatus (ground truth) for evaluation (default is None).
- `norm` (optional): bool, default = True; if raw counts are unavailable, or an error occurs during execution, user may turn normalization off by setting norm = False

**Returns**:
- A list of dictionaries containing donor IDs, predictions, and probabilities (CMV status).
- Optionally, the DataFrame of cell type fractions if `return_frac=True`.

### `visualize(results, visit_order=None, figWidth=8, figHeight=3, dpi_param=100, save=False, filename='cmverify_viz.png', metrics=False)`

This function visualizes longitudinal CMV prediction probabilities per donor.

**Parameters**:
- `results`: A list of dictionaries with keys `'donor_id_timepoint'` (tuple), `'probability'` (float), and optionally `true_label`.
- `visit_order` (optional): list specifying the order of visit labels.
- `figWidth`, `figHeight` (optional): Figure dimensions in inches.
- `dpi_param` (optional): Dots-per-inch resolution for the plot.
- `save` (optional): Whether to save the plot as an image file.
- `filename` (optional): Output filename if saving.
- `metrics` (optional): If True, outputs additional metrics like confusion matrix, roc-curve.

### `append_status(intermed_cmv_predictions, cmv_df, patient_col, cmv_col)`

This utility function appends true CMV status to the intermediate prediction output by matching donor IDs with a reference CMV status DataFrame.
Use this if you have CMV status but it is not in the adata.

**Parameters**:
- `intermed_cmv_predictions`: A list of dictionaries, each containing a `'donor_id_timepoint'` tuple as returned by `predict`.
- `cmv_df` (DataFrame or dict): A `pandas.DataFrame` or dict containing CMV serostatus for each donor.
- `patient_col`: The name of the column in `cmv_df` that contains donor/patient IDs.
- `cmv_col`: The name of the column in `cmv_df` that contains CMV status values (e.g., 0 for negative, 1 for positive).

**Returns**:
- An updated list of dictionaries with a new key `'true_label'` representing the normalized CMV status for each donor. If a donor is not found in `cmv_df`, `'true_label'` is set to `None`.

Model Training
--------------

CMVerify uses a random forest classifier (`rf_best_estimator`) and a corresponding scaler (`rf_scaler`). These models have been trained on relevant single-cell RNA-seq data and are used to predict CMV serostatus based on cell type composition.

Many thanks to the Allen Institute for sharing their cohort data and enhancing reproducibility in science. Below are some helpful links to their preprint, analysis pipeline, scRNA-seq data, and celltypist models.

- [Gong et al. 2024 preprint via bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.10.612119v1)
- [Immune Health Atlas Analysis](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/analysis/)
- [scRNA-seq Downloads – Dynamics of Immune Health with Age](https://apps.allenimmunology.org/aifi/insights/dynamics-imm-health-age/downloads/scrna/)
- [Model Downloads – Immune Health Atlas](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/)

Contributing
------------

If you'd like to contribute to the development of CMVerify, please fork the repository and submit pull requests with proposed changes. All contributions are welcome!

License
-------

CMVerify is released under the AGPL-3.0 License. See LICENSE for more information.

---

Thank you for using CMVerify! If you have any issues or questions, please feel free to open an issue on the repository or contact me at luke.trinity@jax.org
