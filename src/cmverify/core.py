# src/CMVerify/core.py
import scanpy as sc
from .utils import normalize_total_10k, log1p_if_needed
from .annotation import annotate_with_model, check_and_add_labels, calculate_cell_type_fractions
from .models import load_model

def analyze_data(adata,donor_obs_column):
    """Normalize to 10k, apply log1p, and load the model."""
    
    # Normalize the data to 10k reads per cell
    print("Checking if normalizing the data to 10k reads per cell is needed...")
    normalize_total_10k(adata)
    
    # Apply log1p transformation if needed
    print("Checking if log1p transformation if necessary...")
    log1p_if_needed(adata)
    

    model_name = 'AIFI_L3'
    # Annotate the data using the loaded model
    print(f"Annotating the data using the {model_name} model...this can take awhile (after the 'Prediction Done!' step)")
    label_results = annotate_with_model(adata, model_name)

    # Check and add labels, and print the summary
    print(f"Checking label summary for {model_name}...")
    check_and_add_labels(adata, label_results, model_name)

    # Calculate the fraction of cells for each label per patient (person)
    print(f"Calculating the fraction of cells for each label per patient using the {model_name} model...")
    fractions_df = calculate_cell_type_fractions(adata, model_name, donor_obs_column)
    
    # Display the calculated fractions
    print(fractions_df)
    
def load_models():
    rf_best_model = load_model('rf_best_estimator')
    rf_model = load_model('rf_model')
    scaler = load_model('rf_scaler')