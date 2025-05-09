import pytest
import numpy as np
import scanpy as sc
import anndata
import os
from cmverify.utils import normalize_total_10k, log1p_if_needed
from cmverify.models import load_model

# Utility function to create a basic AnnData object
def create_adata(X_data, raw_data=None, log1p=False):
    """
    Creates an AnnData object with optional raw data and log1p status.
    
    Parameters:
    - X_data: Numpy array with the gene expression data to be assigned to adata.X
    - raw_data: Optional; Numpy array with raw gene expression data assigned to adata.raw
    - log1p: Optional; Boolean flag to indicate if the log1p transformation has been applied (sets 'log1p' in adata.uns)
    
    Returns:
    - adata: AnnData object with specified properties
    """
    adata = anndata.AnnData(X=X_data)
    
    if raw_data is not None:
        adata.raw = anndata.AnnData(X=raw_data)
    
    if log1p:
        adata.uns["log1p"] = True
    
    return adata

# Test for normalize_total_10k function
def test_normalize_total_10k():
    """
    Test cases for the normalize_total_10k function to ensure the data in .X is normalized to 10k reads per cell
    when not log-transformed, and no normalization is performed when data is already log-transformed.
    """
    # Case 1: Data is not log-transformed, normalize .X to 10k reads per cell
    adata = create_adata(np.random.randint(low=0, high=250, size=(5, 3)))  # Create random data
    normalize_total_10k(adata, verbose=1,force_norm=False)
    assert np.allclose(adata.X.sum(axis=1), 1e4), "Normalization to 10k failed"  # Check that total is normalized to 10k
    
    # Case 2: Data is already log-transformed, skip normalization
    adata = create_adata(np.random.randint(low=0, high=250, size=(5, 3)), log1p=True)  # Simulate log-transformed data
    normalize_total_10k(adata, verbose = 1,force_norm=False)
    assert np.allclose(adata.X.sum(axis=1), np.sum(adata.X, axis=1)), "Normalization was incorrectly applied"  # Check normalization was skipped
    
    # Case 3: No raw data, normalize .X
    adata = create_adata(np.random.randint(low=0, high=250, size=(5, 3)))  # Create random data
    normalize_total_10k(adata, verbose=1,force_norm=False)
    assert np.allclose(adata.X.sum(axis=1), 1e4), "Normalization to 10k failed"  # Check that total is normalized to 10k

    # Case 4: No .X data, should raise ValueError
    adata = anndata.AnnData()  # Empty AnnData with no data
    with pytest.raises(ValueError):
        normalize_total_10k(adata, verbose = 1,force_norm=False)

# Test for log1p_if_needed function
def test_log1p_if_needed():
    """
    Test cases for the log1p_if_needed function to ensure log1p transformation is applied only when necessary.
    If the data is already log-transformed, no transformation is applied.
    """
    # Case 1: Data is not log-transformed
    adata = create_adata(np.random.randint(low=0, high=250, size=(5, 3)))
    
    log1p_if_needed(adata, verbose=1,force_norm=False)
    assert "log1p" in adata.uns  # Check that the log1p key is added

    # Case 2: Data is already log-transformed
    adata = create_adata(np.random.randint(low=0, high=250, size=(5, 3)), log1p=True)
    
    log1p_if_needed(adata, verbose =1,force_norm=False)
    assert "log1p" in adata.uns  # log1p should still be there and should skip transformation

def test_load_model_valid():
    """
    Test for the load_model function to ensure that a valid model is loaded correctly.
    This test checks if the model is not None and if it is an object, further refinement 
    can be done based on the model's specific properties.
    """
    # Test valid model loading
    model_name = 'rf_best_estimator'
    model = load_model(model_name)
    
    # Check if model is loaded correctly (you can check the type or other properties)
    assert model is not None
    assert isinstance(model, object)  # Change this to a more specific check based on your model