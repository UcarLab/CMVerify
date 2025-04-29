# src/CMVerify/core.py
import scanpy as sc
from .utils import normalize_total_10k, log1p_if_needed
from .annotation import annotate_with_model, check_and_add_labels, calculate_cell_type_fractions
from .models import load_model
import pandas as pd
from IPython.display import display

def analyze_data(adata,donor_obs_column):
    """Normalize to 10k, apply log1p, and load the model."""
    
    # Normalize the data to 10k reads per cell
    print("Checking if normalizing the data to 10k reads per cell is needed...", flush=True)
    normalize_total_10k(adata)
    
    # Apply log1p transformation if needed
    print("Checking if log1p transformation if necessary...", flush=True)
    log1p_if_needed(adata)
    

    model_name = 'AIFI_L3'
    # Annotate the data using the loaded model
    print(f"Annotating the data using the {model_name} model...", flush=True)
    label_results = annotate_with_model(adata, model_name)

    # Check and add labels, and print the summary
    print(f"Checking label summary for {model_name}...", flush=True)
    check_and_add_labels(adata, label_results, model_name)

    # Calculate the fraction of cells for each label per patient (person)
    print(f"Calculating the fraction of cells for each label per donor...", flush=True)
    fractions_df = calculate_cell_type_fractions(adata, model_name, donor_obs_column)
    
    # Display the calculated fractions
    display(fractions_df.head().style.hide(axis='index'))
    
def load_models():
    rf_best_model = load_model('rf_best_estimator')
    rf_model = load_model('rf_model')
    scaler = load_model('rf_scaler')