# src/cmverify/utils.py
import scanpy as sc

def normalize_total_10k(adata):
    """
    Normalize the data in `adata.X` to a total count of 10,000 reads per cell.
    This function assumes that `adata.X` contains the raw counts (no `.raw` check).
    If the data is already log-transformed (as indicated by `log1p` in `adata.uns`), 
    skip normalization to 10k reads per cell.
    """
    # Check if data is log-transformed
    if "log1p" in adata.uns:
        # Skip normalization if already log-transformed
        print("WARNING! Data looks to be already log-transformed. Skipping normalization to 10k reads per cell. Double check your pipeline.")
    elif adata.X is not None:
        # Normalize .X to 10k reads per cell
        print("Normalizing .X to 10k reads per cell.")
        sc.pp.normalize_total(adata, target_sum=1e4)
    else:
        # If .X is not available, raise an error
        raise ValueError("No expression data (.X) found in the AnnData object.")


def log1p_if_needed(adata):
    """
    Apply log1p transformation to `adata.X` if not already log-transformed.
    The function assumes that log1p is applied if 'log1p' exists in `.uns`.
    If `.raw` is present, the log1p will be applied to `.raw` and not `.X`.
    """
    if "log1p" not in adata.uns:
        print("Applying log1p transformation to the data.")
        sc.pp.log1p(adata)
    else:
        print("Data already log-transformed, skipping log1p.")
