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
- `scikit-learn`
- `matplotlib`
- `numpy`
- `pandas`
- `seaborn`
- `anndata`
- `celltypist`
- `ipython`

These dependencies are automatically installed when you install CMVerify.

Usage
-----

Here is how you can use CMVerify in your Python environment:

### 1. Importing CMVerify

To use CMVerify, start by importing the necessary module in your Python script:

    from cmverify import predict, visualize

### 2. Data Preparation

Before running predictions, ensure your data is prepared correctly. CMVerify requires **raw counts** for accurate predictions.

#### Converting `.raw` to `.X`:

You can convert your AnnData object from `.raw` to `.X` as follows:

    import scanpy as sc

    # Read the .h5ad file containing the raw data
    adata_pre = sc.read_h5ad('path_to_data.h5ad')

    # Use the raw counts for analysis
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

Functions
---------

### `predict(adata, donor_obs_column, longitudinal_obs_column=None, verbose=1, return_frac=False)`

This is the main function used to predict CMV status. It performs the following steps:
1. Normalizes the data to 10k reads per cell.
2. Applies log1p transformation if needed.
3. Annotates the data with the specified model.
4. Calculates the cell type fractions per donor.
5. Loads the models and scaler.
6. Makes predictions based on the cell type fractions.

Parameters:
- **adata**: An AnnData object containing your single-cell RNA-seq data.
- **donor_obs_column**: The column in `adata.obs` that contains the donor ID.
- **longitudinal_obs_column** (optional): The column in `adata.obs` for timepoints (e.g., for longitudinal data).
- **verbose**: Set to 1 for progress messages.
- **return_frac**: If `True`, returns the fractions DataFrame along with predictions.

Returns:
- A list of dictionaries containing donor IDs, predictions, and probabilities (CMV status).
- Optionally, the DataFrame of cell type fractions.

### 5. Error Handling

If a required column is missing in your `adata` object, an error will be raised:

    ValueError: {donor_obs_column} is not a valid column in adata.obs.

Model Training
--------------

CMVerify uses a pre-trained random forest model (`rf_best_estimator`) and a corresponding scaler (`rf_scaler`). These models have been trained on relevant single-cell RNA-seq data and are used to predict CMV serostatus based on cell type composition.

Many thanks to the Allen Institute for sharing their cohort data and enhancing reproducability in science. Below are some helpful links to their analysis pipeline, scRNA-seq data, and celltypist models.
- [Gong et al. 2024 preprint via bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.10.612119v1) – bioRxiv preprint
- [Immune Health Atlas Analysis](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/analysis/)
- [scRNA-seq Downloads – Dynamics of Immune Health with Age](https://apps.allenimmunology.org/aifi/insights/dynamics-imm-health-age/downloads/scrna/)
- [Model Downloads – Immune Health Atlas](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/)

Contributing
------------

If you'd like to contribute to the development of CMVerify, please fork the repository and submit pull requests with proposed changes. All contributions are welcome!

License
-------

CMVerify is released under the MIT License. See LICENSE for more information.

---

Thank you for using CMVerify! If you have any issues or questions, please feel free to open an issue on the repository or contact me at luke.trinity@jax.org
