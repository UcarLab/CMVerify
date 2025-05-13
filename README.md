CMVerify: CMV Serostatus Predictor
===================================

CMVerify is a tool designed to predict Cytomegalovirus (CMV) serostatus based on single-cell RNA sequencing data. It uses a machine learning model to analyze blood transcriptomic data and provide predictions on CMV status. The tool is built to be easy to use, with automatic data preprocessing, model annotation, and prediction generation.
- [Ucar Lab at The Jackson Laboratory](https://ucarlab.github.io/)

Installation
------------

To install CMVerify, you can use `pip`:

     pip install CMVerify

Alternatively, if you'd like to install the package from source, the repo will be made public upon publication:

    git clone https://github.com/UcarLab/CMVerify.git
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

R Users
-------
[For R users please refer to this vignette for converting from Seurat to Anndata](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)
Further support for R users coming in version 2.

Usage
-----

Here is how you can use CMVerify in your Python environment:

### 1. Importing CMVerify

To use CMVerify, start by importing the necessary module in your Python script:

```python
from cmverify import predict, visualize
```

### 2. Data Preparation

Before running predictions, ensure your data is prepared correctly. CMVerify requires **raw counts** for accurate predictions.

#### If raw counts are in `.raw`, convert `.raw` to `.X`:

You can convert your AnnData object from `.raw` to `.X` as follows:

```python
import scanpy as sc

# Read the .h5ad file containing the raw data
adata_pre = sc.read_h5ad('path_to_data.h5ad')

# Use the raw counts for analysis (and retain metadata in .obs and .var)
adata = sc.AnnData(X=adata_pre.raw.X, obs=adata_pre.obs, var=adata_pre.raw.var)
```

Please ensure that the data used for predictions is based on raw counts, which is essential for the analysis.

### 3. Making Predictions

To make predictions, you need to load your single-cell RNA-seq data (`AnnData` object) and provide the relevant parameters. The `predict` function handles data normalization, model annotation, and CMV status prediction.

#### Example:

```python
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
```

### 4. Return Fractions

You can also return the calculated cell type fractions along with the predictions by setting the `return_frac` parameter to `True`.

```python
results, fractions_df = predict(adata, donor_obs_column, longitudinal_obs_column, return_frac=True)

# Display the first 5 rows of the cell type fractions
print(fractions_df.head())
```

#### 4b. Append Ground Truth CMV Status (Optional)

If you have true CMV status in a separate metadata file (not in `adata.obs`), you can use the append_status function to match and append it to the prediction output. This method supports dataframe with patient/cmv status column or dict input type (key=patient, value=CMV status).

```python
from cmverify import append_status

# Add true labels to predictions
append_status(results, cmv_metadata)
```


### 5. Visualize Longitudinal Predictions (and optionally test the model with ground truth)

You can visualize longitudinal CMV prediction probabilities across timepoints using the `visualize` function.

```python
from cmverify import visualize

visualize(results)
```

This will generate a figure connecting donor predictions across visits and mark the decision threshold.
You can also set `metrics=True` if you have ground truth CMV serostatus and wish to evaluate this model.

Functions
---------

### `predict(adata, donor_obs_column, longitudinal_obs_column=None, verbose=1, return_frac=False, true_status=None, norm=True, force_norm=False)`

Predict CMV serostatus from single-cell RNA-seq data using CMVerify. This function handles normalization, transformation, annotation, cell type fraction calculation, and model inference.

**Parameters**:
- `adata` (`AnnData`): Single-cell RNA-seq data object.
- `donor_obs_column` (`str`): The column in `adata.obs` that contains the donor or sample identifiers.
- `longitudinal_obs_column` (`str`, optional): Column in `adata.obs` for timepoint or longitudinal visit labels, if applicable. Default is `None`.
- `verbose` (`int`, optional): Verbosity level for progress messages. (0 = silent, 1 = standard output). Default is 1.
- `return_frac` (`bool`, optional): If `True`, return the cell type fraction `DataFrame` along with predictions. Default is `False`.
- `true_status` (`str`, optional): Column in `adata.obs` for true donor serostatus (ground truth) for evaluation. Default is `None`.
- `norm` (`bool`, optional): Whether to normalize to 10,000 counts per cell and log1p-transform. Disable only if raw counts are in `adata.X`. Default is `True`.
- `force_norm` (`bool`, optional): If the adata has the log1p layer but has not been normalized, user may encounter error from celltypist annotation step. Set `force_norm=True` to force the normalization and resolve the issue. Default is False.

**Returns**:
- `List[Dict]`: List of dictionaries containing donor ID, timepoint (if applicable), predicted label, and probability of CMV seropositivity from CMVerify.
- If `return_frac=True`, also returns a DataFrame with  cell type fractions and donor_id-timepoint metadata.

### `predict_from_frac(fractions_df, verbose=1)`

Apply CMVerify to a precomputed DataFrame of cell type fractions.

**Parameters**:
- `fractions_df` (`DataFrame`): A DataFrame where rows are donor-timepoints and columns are cell type fractions. The last column must contain the donor-timepoint info as a tuple (`donor_id`,`timepoint`).
- `verbose` (`int`, optional): Verbosity level for progress messages. (0 = silent, 1 = standard output). Default is 1.

**Returns**:
- `List[Dict]`: List of dictionaries containing donor ID, timepoint (if applicable), predicted label, and probability of CMV seropositivity from CMVerify.


### `visualize(results, visit_order=None, figWidth=8, figHeight=3, dpi_param=100, save=False, filename='cmverify_viz.png', metrics=False)`

This function visualizes CMV prediction probabilities per donor, and if applicable for each timepoint.

**Parameters**:
- `results` (`List[Dict]`): A list of dictionaries with keys `'donor_id_timepoint'`, i.e., tuple: (`donor_id`,`timepoint`), `'probability'` (`float`), and optionally `true_label`.
- `visit_order` (`List[str]`, optional): List specifying the order of visit labels (e.g., `["Baseline", "Visit 1", "Visit 2"]`). Default is `None`.
- `figWidth`, `figHeight` (`float`, optional): Figure dimensions in inches. Default is 8 x 3.
- `dpi_param` (`int`, optional): Dots-per-inch resolution for the plot. Default is 100.
- `save` (`bool`, optional): If `True`, saves the plot as an image file. Default is `False`.
- `filename` (`str`, optional): Output filename if saving. Default is ending with `cmverify_viz.png` (multiple figures will be generated with different prefix).
- `metrics` (`bool`, optional): If `True`, outputs additional metrics like confusion matrix, roc-curve (requires `true_label` in results). Default is `False`.

### `append_status(intermed_cmv_predictions, cmv_df, patient_col='patientID', cmv_col='CMV')`

This utility function appends true CMV status to the intermediate prediction output by matching donor IDs with a reference CMV status DataFrame.
Use this if you have CMV status but it is not in the adata.

**Parameters**:
- `intermed_cmv_predictions` (`List[Dict]`): Output from `predict`, contains `'donor_id_timepoint'` (tuple).
- `cmv_df` (`DataFrame` or `dict`): A `pandas.DataFrame` or dict containing known CMV serostatus for each donor.
- `patient_col` (`str`, optional): The name of the column in `cmv_df` that contains donor/patient IDs.
- `cmv_col` (`str`, optional): The name of the column in `cmv_df` that contains CMV status values (e.g., 0 for negative, 1 for positive).

**Returns**:
- Updates intermediate predictions in place with new key `true_label` added to each dictionary entry.

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
