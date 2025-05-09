# src/CMVerify/core.py
from .utils import normalize_total_10k, log1p_if_needed, normalize_cmv, is_log1p_transformed
from .annotation import annotate_with_model, check_and_add_labels, calculate_cell_type_fractions, check_for_missing_columns
from .models import load_model
from .config import EXPECTED_COLUMNS
import pandas as pd
import numpy as np
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

def predict(adata,donor_obs_column, longitudinal_obs_column=None, verbose = 1,return_frac=False, true_status=None, norm = True, force_norm=False):
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
    - force_norm : bool, default = False
        If the adata has the log1p layer but has not been normalized, user may encounter error from celltypist annotation step. Set force_norm=True to force the normalization and resolve the issue.

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
        normalize_total_10k(adata,verbose,force_norm)
       
        if verbose == 1:
            print("Checking if log1p transformation is necessary...", flush=True)
        # Apply log1p transformation if needed
        log1p_if_needed(adata, verbose,force_norm)
    else:
        print("User turned off normalization, data should already be normalized to 10k reads and log1p...", flush=True)
        # provide some info (best guess) if the data is log already)
        is_log, min_val, max_val = is_log1p_transformed(adata)
        print(f"Log1p-transformed: {is_log} (min={min_val}, max={max_val})")
    
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
    fractions_df, donor_ids_timepoints = calculate_cell_type_fractions(adata, model_name, donor_obs_column, longitudinal_obs_column, verbose)

    # Display the calculated fractions
    if verbose == 1:
        print(f"Displaying first 5 rows of donor level peripheral blood mononuclear cell composition:", flush=True)
        display(fractions_df.head().style.hide(axis='index'))

    results = predict_from_model(fractions_df, donor_ids_timepoints, verbose)

    # Return fractions df as well, with donor_id_timepoint, pred and prob
    fractions_df["donor_id_timepoint"] = donor_ids_timepoints
    
    # Reorder columns
    cols = ["donor_id_timepoint"] + [col for col in fractions_df.columns if col not in ["donor_id_timepoint"]]

    if true_status is not None:
        if verbose:
            print("Getting true status from adata and standardizing CMV labels", flush=True)
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

def predict_from_frac(fractions_df, verbose = 1):
    """
    Predicts donor classification (e.g., CMV status) from precomputed cell type fractions.

    This function is useful when you already have a DataFrame of per-donor (or per-donor-timepoint) 
    cell type fractions and want to apply a pre-trained model to get classification predictions 
    and probabilities.

    Parameters:
    - fractions_df : pd.DataFrame
        A DataFrame where rows correspond to donors (or donor-timepoint pairs) and columns are 
        cell type fractions. The last column must be the donor ID or a tuple identifier.

    Returns:
    - List[dict]
        A list of dictionaries containing donor IDs, predicted labels, and predicted probabilities.
    """
    # Extract and save the last column
    donor_ids_timepoints = fractions_df.iloc[:, -1]
    
    # Drop the last column
    fractions_df_clean = fractions_df.iloc[:, :-1]

    fractions_df_clean = check_for_missing_columns(fractions_df_clean)
    
    results = predict_from_model(fractions_df_clean, donor_ids_timepoints, verbose)
    
    if verbose:
        print("Outputting predictions", flush=True)
        print(results)
        print("All done. Thank you!", flush=True)
    return results

def predict_from_model(fractions_df, donor_ids_timepoints, verbose):
    """
    Predicts CMV status using a pre-trained random forest model.

    This function takes in a dataframe of cell type fractions and a list of donor ID/timepoint identifiers,
    scales the input using a pre-loaded scaler, and applies a pre-trained random forest model to predict 
    CMV status. It returns predictions and associated probabilities for each donor-timepoint pair.

    Parameters
    ----------
    fractions_df : pandas.DataFrame
        DataFrame containing immune cell type fractions per donor-timepoint.
    
    donor_ids_timepoints : list of str
        List of donor IDs with associated timepoints, matching the rows in `fractions_df`.

    Returns
    -------
    list of dict
        A list where each element is a dictionary containing:
        - 'donor_id_timepoint': the donor-timepoint identifier
        - 'prediction': predicted CMV status (0 or 1)
        - 'probability': predicted probability of CMV-positive (rounded to 3 decimals)
    """
    # Load pre-trained random forest model and corresponding scaler
    rf_best_model, scaler = load_models(verbose)

    if verbose == 1:
        print("Scaling the fractions...", flush=True)
    # Scale the fractions data using the pre-loaded scaler
    fractions_df_scaled = scaler.transform(fractions_df)
    
    if verbose == 1:
        print("Making predictions using the CMVerify model...", flush=True)
    # Get the predictions (CMV status)
    cmv_pred = rf_best_model.predict(fractions_df_scaled)

    if verbose == 1:
        print("Getting predicted probabilities for CMV status...", flush=True)
    # Get the predicted probabilities for CMV status
    cmv_pred_probs = np.round(rf_best_model.predict_proba(fractions_df_scaled)[:, 1],2)  # Probability of the positive class
    
    # Combine the donor ID, prediction, and probability into a list of dictionaries
    results = []
    for donor_id_tp, pred, prob in zip(donor_ids_timepoints, cmv_pred, cmv_pred_probs):
        results.append({
            'donor_id_timepoint': donor_id_tp,
            'prediction': pred,
            'probability': round(prob,3)
        })
    return results
    

def append_status(intermed_cmv_predictions, cmv_df, patient_col='patientID', cmv_col='CMV'):
    """
    Appends normalized CMV status to intermed_cmv_predictions list of dictionaries.

    Parameters:
    - intermed_cmv_predictions (list of dict): List of dictionaries, each containing 'donor_id_timepoint'.
    - cmv_df (DataFrame or dict): DataFrame or dictionary containing 'patientID' and 'CMV' columns.
    - patient_col (str, optional): Column in cmv_df with donor ID. Default: 'patientID'.
    - cmv_col (str, optional): Column in cmv_df with CMV status. Default: 'CMV'.

    Returns:
    - None: The function updates intermed_cmv_predictions in place.
    """
    # Check if cmv_df is a dictionary
    if isinstance(cmv_df, dict):
        # Convert the dictionary to a DataFrame
        cmv_df = pd.DataFrame(list(cmv_df.items()), columns=[patient_col, cmv_col])

    # Cast donor IDs in DataFrame to string for reliable comparison
    cmv_df[patient_col] = cmv_df[patient_col].astype(str)
    
    # Loop through each dictionary in the predictions list
    for d in intermed_cmv_predictions:
        # Extract the donor_id from the dictionary (first item in donor_id_timepoint tuple)
        donor_id = str(d['donor_id_timepoint'][0])
        
        # Look up the CMV value from the DataFrame based on the donor_id (patientID)
        cmv_value = cmv_df.loc[cmv_df[patient_col] == donor_id, cmv_col].values
        if len(cmv_value) > 0:
            # Normalize the CMV value and add it as 'true_label'
            d['true_label'] = normalize_cmv(cmv_value[0])
        else:
            print("Warning, there is a donor with no CMV status. Metrics may not run correctly.")
            # If no match is found, you can choose to add None or handle the error
            d['true_label'] = None
