# src/CMVerify/core.py
from .utils import normalize_total_10k, log1p_if_needed
from .annotation import annotate_with_model, check_and_add_labels, calculate_cell_type_fractions
from .models import load_model
try:
    from IPython.display import display  # Try importing display for Jupyter environments
except ImportError:
    display = print  # Use print as fallback in non-IPython environments

def load_models():
    """Load the models and scaler for use in predictions."""
    print("Loading model...", flush=True)

    # Load the models and scaler
    rf_best_model = load_model('rf_best_estimator')
    scaler = load_model('rf_scaler')

    # Return the models and scaler so they can be used in analysis
    return rf_best_model, scaler

def predict(adata,donor_obs_column, longitudinal_obs_column=None, verbose = False):
    """Normalize to 10k, apply log1p, load the models, annotate and predict."""
    # Confirm required parameters
    if donor_obs_column not in adata.obs.columns:
        raise ValueError(f"{donor_obs_column} is not a valid column in adata.obs.")
    
    # Normalize the data to 10k reads per cell
    print("Checking if normalizing the data to 10k reads per cell is needed...", flush=True)
    normalize_total_10k(adata)
    
    # Apply log1p transformation if needed
    print("Checking if log1p transformation is necessary...", flush=True)
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
    fractions_df, donor_ids = calculate_cell_type_fractions(adata, model_name, donor_obs_column, longitudinal_obs_column)
    
    # Display the calculated fractions
    if verbose:
        print(f"Displaying first 5 rows of donor level peripheral blood mononuclear cell composition:", flush=True)
        display(fractions_df.head().style.hide(axis='index'))

    # Load models and scaler
    rf_best_model, scaler = load_models()

    # Scale the fractions data using the pre-loaded scaler
    print("Scaling the calculated fractions...", flush=True)
    fractions_df_scaled = scaler.transform(fractions_df)
    
    # Get the predictions (CMV status)
    print("Making predictions using the CMVerify model...", flush=True)
    cmv_pred = rf_best_model.predict(fractions_df_scaled)

    # Get the predicted probabilities for CMV status
    print("Getting predicted probabilities for CMV status...", flush=True)
    cmv_pred_probs = rf_best_model.predict_proba(fractions_df_scaled)[:, 1]  # Probability of the positive class
    
    # Combine the donor ID, prediction, and probability into a list of dictionaries
    results = []
    for donor_id, pred, prob in zip(donor_ids, cmv_pred, cmv_pred_probs):
        results.append({
            'donor_id_timepoint': donor_id,
            'prediction': pred,
            'probability': round(prob,3)
        })

    if verbose:
        print("Outputting predictions", flush=True)
        print(results)
        
    print("All done. Thank you!", flush=True)
    return results