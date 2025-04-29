# src/CMVerify/core.py
import scanpy as sc
from .utils import normalize_total_10k, log1p_if_needed
from .models import load_model

def analyze_data(adata):
    """Normalize to 10k, apply log1p, and load the model."""
    
    # Normalize the data to 10k reads per cell
    print("Normalizing the data to 10k reads per cell...")
    normalize_total_10k(adata)
    
    # Apply log1p transformation if needed
    print("Applying log1p transformation if necessary...")
    log1p_if_needed(adata)
    
    # Load the pre-trained model (for example 'AIFI_L3')
    print("Loading the pre-trained cell typist model...")
    model = load_model('AIFI_L3')
    
    # Return processed data and model for further use
    return adata, model
