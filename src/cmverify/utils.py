# src/cmverify/utils.py
import scanpy as sc
import pandas as pd
import numpy as np

def normalize_total_10k(adata, verbose, force_norm):
    """
    Normalize the data in `adata.X` to a total count of 10,000 reads per cell.
    This function assumes that `adata.X` contains the raw counts (no `.raw` check).
    If the data is already log-transformed (as indicated by `log1p` in `adata.uns`), 
    skip normalization to 10k reads per cell.
    """
    # Check if data is log-transformed
    if "log1p" in adata.uns and not force_norm:
        # Skip normalization if already log-transformed
        print("WARNING! Data looks to be already log-transformed (log1p layer detected in adata.uns). Skipping normalization to 10k reads per cell. Double check your pipeline.", flush=True)
    elif adata.X is not None:
        if verbose == 1:
            print("Normalizing .X to 10k reads per cell.", flush=True)
        # Normalize .X to 10k reads per cell
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        # If .X is not available, raise an error
        raise ValueError("No expression data (.X) found in the AnnData object.")

def log1p_if_needed(adata, verbose, force_norm):
    """
    Apply log1p transformation to `adata.X` if not already log-transformed.
    The function assumes that log1p is applied if 'log1p' exists in `.uns`.
    If `.raw` is present, the log1p will be applied to `.raw` and not `.X`.
    """
    # provide some info (best guess) if the data is log already)
    is_log, min_val, max_val = is_log1p_transformed(adata)
    print(f"Current state (best guess) Log1p-transformed: {is_log} (min={min_val}, max={max_val})")
    if "log1p" not in adata.uns or force_norm:
        sc.pp.log1p(adata)
        if verbose == 1:
            print("Applied log1p transformation to the data.", flush=True)
            is_log, min_val, max_val = is_log1p_transformed(adata)
            print(f"Log1p-transformed: {is_log} (min={min_val}, max={max_val})")
    else:
        print("Data looks to be already log-transformed (log1p layer detected in adata.uns), skipping log1p.", flush=True)

def normalize_cmv(value):
    """Convert CMV labels to binary format (0 = negative, 1 = positive)."""
    if pd.isna(value):
        return None
    value_str = str(value).strip().lower()
    if value_str in {'1', '1.0', 'pos', 'positive','p','+'}:
        return 1
    elif value_str in {'0', '0.0', 'neg', 'negative','n','-'}:
        return 0
    else:
        print(f"Unrecognized CMV label: '{value}'. Expected values are variants of [0, 1, 0.0, 1,0, 'positive','negative','pos', 'neg','+','-'].")
        return None  # Or raise an error if strict

def is_log1p_transformed(adata, threshold=20):
    """
    Checks whether the expression matrix in an AnnData object appears to be log1p-transformed.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix (typically from Scanpy).
    threshold : float, optional (default=20)
        The upper value threshold used to infer whether the data has been log-transformed.
        If the maximum value in the matrix is below this threshold, it is assumed to be log1p-transformed.

    Returns
    -------
    is_log : bool
        True if the data is likely log1p-transformed, False otherwise.
    min_val : float
        Minimum value in the expression matrix.
    max_val : float
        Maximum value in the expression matrix.

    Notes
    -----
    This function makes a heuristic guess. It assumes:
    - Raw UMI counts often have max values > 100 or even > 1000.
    - log1p(counts) typically yields values between 0 and ~10.
    """
    X = adata.X

    if hasattr(X, "min"):  # sparse matrix with .min and .max methods
        min_val = X.min()
        max_val = X.max()
    else:  # dense array or needs toarray conversion
        X = X.toarray() if hasattr(X, "toarray") else X
        min_val = np.min(X)
        max_val = np.max(X)

    is_log = max_val < threshold and min_val >= 0
    return is_log, min_val, max_val
