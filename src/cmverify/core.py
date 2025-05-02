# src/CMVerify/core.py
from .utils import normalize_total_10k, log1p_if_needed, normalize_cmv
from .annotation import annotate_with_model, check_and_add_labels, calculate_cell_type_fractions
from .models import load_model
import pandas as pd
try:
    from IPython.display import display  # Try importing display for Jupyter environments
except ImportError:
    display = print  # Use print as fallback in non-IPython environments

def load_models(verbose):
    """Load the models and scaler for use in predictions."""
    if verbose ==1:
        print("Loading the CMVerify model...", flush=True)

    # Load the models and scaler
    rf_best_model = load_model('rf_best_estimator')
    scaler = load_model('rf_scaler')

    # Return the models and scaler so they can be used in analysis
    return rf_best_model, scaler

def predict(adata,donor_obs_column, longitudinal_obs_column=None, verbose = 1,return_frac=False, true_status=None, norm = True):
    """
    Predicts donor classification from an AnnData object using pre-trained models.
    
    This function performs the following steps:
    1. Normalizes expression data to 10,000 counts per cell.
    2. Applies log1p transformation if not already applied.
    3. Loads the appropriate model(s) for prediction.
    4. Annotates the AnnData object with predicted donor-level outcomes.
    
    Parameters:
    - adata : AnnData
        Single-cell AnnData object containing expression data in `.X` and metadata in `.obs`.
    - donor_obs_column : str
        Column name in `adata.obs` that uniquely identifies each donor (used for aggregation).
    - longitudinal_obs_column : str, optional
        Column name in `adata.obs` indicating longitudinal timepoints for donors, if applicable.
    - verbose : int, default=1
        Verbosity level. Set to 0 for silent mode.
    - return_frac : bool, default=False
        Whether to return the fraction of predictive cell types used for classification.
    - true_status : str or None, default=None
        Optional ground truth donor status column name for evaluation or model validation.
    - norm : bool, default = True
        We highly encourage passing raw counts into this method, however, if raw counts are unavailable, or an error occurs during execution, user may turn normalization off by setting norm = False

    Returns:
    - AnnData
        The input AnnData object with added predictions in `adata.obs` or `adata.uns`.
        If `return_frac` is True, also returns a DataFrame with the fractions of predictive cell types.
    """
    # Confirm required parameters
    if donor_obs_column not in adata.obs.columns:
        raise ValueError(f"{donor_obs_column} is not a valid column in adata.obs.")
    if longitudinal_obs_column is not None and longitudinal_obs_column not in adata.obs.columns:
            raise ValueError(f"{longitudinal_obs_column} is not a valid column in adata.obs.")
    if true_status is not None and true_status not in adata.obs.columns:
            raise ValueError(f"{true_status} is not a valid column in adata.obs.")

    if norm:
        if verbose == 1:
            print("Checking if normalizing the data to 10k reads per cell is needed...", flush=True)
        # Normalize the data to 10k reads per cell
        normalize_total_10k(adata,verbose)
       
        if verbose == 1:
            print("Checking if log1p transformation is necessary...", flush=True)
        # Apply log1p transformation if needed
        log1p_if_needed(adata, verbose)
    else:
        print("User turned off normalization, data should already be normalized to 10k reads and log1p...", flush=True)
    
    model_name = 'AIFI_L3'
    if verbose == 1:
        print(f"Annotating the data using the {model_name} model...", flush=True)
    # Annotate the data using the loaded model
    label_results = annotate_with_model(adata, model_name)

    
    if verbose == 1:
        print(f"Checking label summary for {model_name}...", flush=True)
    # Check and add labels, and print the summary
    check_and_add_labels(adata, label_results, model_name, verbose)

    if verbose == 1:
        print(f"Calculating the fraction of cells for each label per donor...", flush=True)
    # Calculate the fraction of cells for each label per patient (person)
    fractions_df, donor_ids = calculate_cell_type_fractions(adata, model_name, donor_obs_column, longitudinal_obs_column, verbose)

    # Display the calculated fractions
    if verbose == 1:
        print(f"Displaying first 5 rows of donor level peripheral blood mononuclear cell composition:", flush=True)
        display(fractions_df.head().style.hide(axis='index'))

    # Load models and scaler
    rf_best_model, scaler = load_models(verbose)

    if verbose == 1:
        print("Scaling the calculated fractions...", flush=True)
    # Scale the fractions data using the pre-loaded scaler
    fractions_df_scaled = scaler.transform(fractions_df)
    
    if verbose == 1:
        print("Making predictions using the CMVerify model...", flush=True)
    # Get the predictions (CMV status)
    cmv_pred = rf_best_model.predict(fractions_df_scaled)

    if verbose == 1:
        print("Getting predicted probabilities for CMV status...", flush=True)
    # Get the predicted probabilities for CMV status
    cmv_pred_probs = rf_best_model.predict_proba(fractions_df_scaled)[:, 1]  # Probability of the positive class
    
    # Combine the donor ID, prediction, and probability into a list of dictionaries
    results = []
    for donor_id, pred, prob in zip(donor_ids, cmv_pred, cmv_pred_probs):
        results.append({
            'donor_id_timepoint': donor_id,
            'prediction': pred,
            'probability': round(prob,3)
        })

    # Return fractions df as well, with donor_id_timepoint, pred and prob
    fractions_df["donor_id_timepoint"] = donor_ids
    
    # Reorder columns
    cols = ["donor_id_timepoint"] + [col for col in fractions_df.columns if col not in ["donor_id_timepoint"]]

    if true_status is not None:
        if verbose:
            print("Standardizing CMV labels", flush=True)
        # Create the CMV dictionary
        df = adata.obs[[donor_obs_column, true_status]].copy()
        df['cmv'] = df[true_status].apply(normalize_cmv)
        df = df.dropna().drop_duplicates(subset=donor_obs_column)
        true_labels = df.set_index(donor_obs_column)['cmv'].astype(int).to_dict()   
        # Assuming 'true_labels' has been already created as a dictionary
        for result in results:
            donor_id = result['donor_id_timepoint'][0]  # Assuming donor_id_timepoint is a tuple (donor_id, timepoint)
            result['true_label'] = true_labels.get(donor_id, None)
    if verbose:
        print("Outputting predictions", flush=True)
        print(results)
        print("All done. Thank you!", flush=True)
    if return_frac:
        return results, fractions_df
    else:
        return results

def append_status(intermed_cmv_predictions, cmv_df, patient_col, cmv_col):
    """
    Appends normalized CMV status to intermed_cmv_predictions list of dictionaries.

    Parameters:
    - intermed_cmv_predictions (list of dict): List of dictionaries, each containing 'donor_id_timepoint'.
    - cmv_df (DataFrame or dict): DataFrame or dictionary containing 'patientID' and 'CMV' columns.
    - patient_col (str): Name of the column in cmv_df that contains donor ID.
    - cmv_col (str): Name of the column in cmv_df that contains CMV status.

    Returns:
    - None: The function updates intermed_cmv_predictions in place.
    """
    # Check if cmv_df is a dictionary
    if isinstance(cmv_df, dict):
        # Convert the dictionary to a DataFrame
        cmv_df = pd.DataFrame(list(cmv_df.items()), columns=['key', 'value'])
        cmv_df.columns = [patient_col, cmv_col]  # Rename columns to match expected structure
    
    # Loop through each dictionary in the predictions list
    for d in intermed_cmv_predictions:
        # Extract the donor_id from the dictionary (first item in donor_id_timepoint tuple)
        donor_id = d['donor_id_timepoint'][0]
        
        # Look up the CMV value from the DataFrame based on the donor_id (patientID)
        cmv_value = cmv_df.loc[cmv_df[patient_col] == donor_id, cmv_col].values
        if len(cmv_value) > 0:
            # Normalize the CMV value and add it as 'true_label'
            d['true_label'] = normalize_cmv(cmv_value[0])
        else:
            print("Warning, there is a donor with no CMV status. Metrics may not run correctly.")
            # If no match is found, you can choose to add None or handle the error
            d['true_label'] = None
    return intermed_cmv_predictions
