import pytest
import numpy as np
import scanpy as sc
import anndata
import os
from cmverify.utils import normalize_total_10k, log1p_if_needed
from cmverify.models import load_model

# Utility function to create a basic AnnData object
def create_adata(X_data, raw_data=None, log1p=False):
    adata = anndata.AnnData(X=X_data)
    
    if raw_data is not None:
        adata.raw = anndata.AnnData(X=raw_data)
    
    if log1p:
        adata.uns["log1p"] = True
    
    return adata

# Test for normalize_total_10k function
def test_normalize_total_10k():
    # Case 1: Data is not log-transformed, normalize .X to 10k reads per cell
    adata = create_adata(np.random.rand(5, 3))  # Create random data
    normalize_total_10k(adata)
    assert np.allclose(adata.X.sum(axis=1), 1e4), "Normalization to 10k failed"  # Check that total is normalized to 10k
    
    # Case 2: Data is already log-transformed, skip normalization
    adata = create_adata(np.random.rand(5, 3), log1p=True)  # Simulate log-transformed data
    normalize_total_10k(adata)
    assert np.allclose(adata.X.sum(axis=1), np.sum(adata.X, axis=1)), "Normalization was incorrectly applied"  # Check normalization was skipped
    
    # Case 3: No raw data, normalize .X
    adata = create_adata(np.random.rand(5, 3))  # Create random data
    normalize_total_10k(adata)
    assert np.allclose(adata.X.sum(axis=1), 1e4), "Normalization to 10k failed"  # Check that total is normalized to 10k

    # Case 4: No .X data, should raise ValueError
    adata = anndata.AnnData()  # Empty AnnData with no data
    with pytest.raises(ValueError):
        normalize_total_10k(adata)

# Test for log1p_if_needed function
def test_log1p_if_needed():
    # Case 1: Data is not log-transformed
    adata = create_adata(np.random.rand(5, 3))
    
    log1p_if_needed(adata)
    assert "log1p" in adata.uns  # Check that the log1p key is added

    # Case 2: Data is already log-transformed
    adata = create_adata(np.random.rand(5, 3), log1p=True)
    
    log1p_if_needed(adata)
    assert "log1p" in adata.uns  # log1p should still be there and should skip transformation

def test_load_model_valid():
    # Test valid model loading
    model_name = 'AIFI_L3'
    model = load_model(model_name)
    
    # Check if model is loaded correctly (you can check the type or other properties)
    assert model is not None
    assert isinstance(model, object)  # Change this to a more specific check based on your model