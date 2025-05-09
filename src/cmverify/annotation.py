# src/CMVerify/annotation.py
import celltypist
import os
from .config import EXPECTED_COLUMNS

def annotate_with_model(adata, model_name):
    """
    Load a pre-trained model and annotate cells in the AnnData object.
    
    Parameters:
    adata: AnnData
        The AnnData object to annotate.
    model_name: str
        The name of the model to use for annotation.
    
    Returns:
    dict: A dictionary containing annotated labels and their scores.
    """

    label_results = {}

    # Define the columns for predictions and scores
    label_column = f'{model_name}_prediction'

    #resolve the path to the model
    model_file = 'models/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl'
    model_path = os.path.join(os.path.dirname(__file__), model_file)

    try:
        # Annotate the AnnData object using the selected model
        predictions = celltypist.annotate(adata, model=model_path)
    except ValueError as e:
        if "Input X contains NaN" in str(e):
            print("Error: The input contains NaN values. Reload your adata object and try calling `predict(..., norm=False)` to avoid redundant normalization that may introduce NaNs.")
        else:
            raise  # re-raise other ValueErrors not related to NaNs

    # Extract the predicted labels and rename the column
    labels = predictions.predicted_labels
    labels = labels.rename({'predicted_labels': label_column}, axis=1)

    # Store the results in the label_results dictionary
    label_results[model_name] = labels

    return label_results

def check_and_add_labels(adata, label_results, model_name, verbose):
    """
    Check the summary of predicted labels, add them to `adata.obs`, and check for sufficient unique cell types.
    
    Parameters:
    adata: AnnData
        The AnnData object to update.
    label_results: dict
        Dictionary containing the annotation results.
    model_name: str
        The name of the model used for annotation.
    
    Returns:
    None
    """
    # Extract predicted labels from label results
    label_df = label_results[model_name]  # Since we have only one model

    label_column = f'{model_name}_prediction'
    
    # Print a summary of the predicted labels
    label_summary = label_df[label_column].value_counts()

    # Check if there are at least 20 unique cell types in the predictions
    unique_cell_types = label_summary.shape[0]
    if unique_cell_types < 20:
        print(f"Warning: The model '{model_name}' predicted only {unique_cell_types} unique cell types. You may not have enough cells.")
    else:
        if verbose == 1:
            print("Top 20 most frequent cell types detected with cell counts:")
            print(label_summary[0:20], end = '\n\n')

    # Add the predicted labels and score to adata.obs
    adata.obs[label_column] = label_df[label_column]

def calculate_cell_type_fractions(adata, model_name, donor_obs_column, longitudinal_obs_column, verbose):
    """
    Calculate the fraction of cells for each label per patient (person).
    
    Parameters:
    adata: AnnData
        The AnnData object to calculate fractions for.
    model_name: str
        The name of the model to use for annotation.
    
    Returns:
    pd.DataFrame: DataFrame containing the fractions of each predicted label per patient.
    """
    # Extract relevant columns from `adata.obs` to a dataframe
    label_column = f'{model_name}_prediction'

    # modularizing to allow longitudinal prediction
    if longitudinal_obs_column is not None:
        obs_df = adata.obs[[donor_obs_column, longitudinal_obs_column,label_column]]
        # Calculate the fraction of cells for each label per patient per timepoint
        fractions_df = (
            obs_df.groupby([donor_obs_column, longitudinal_obs_column, label_column])
            .size()
            .unstack(fill_value=0)  # Converts to a wide format with labels as columns
        )
    else:
        obs_df = adata.obs[[donor_obs_column, label_column]]
        # Calculate the fraction of cells for each label per patient
        fractions_df = (
            obs_df.groupby([donor_obs_column, label_column])
            .size()
            .unstack(fill_value=0)  # Converts to a wide format with labels as columns
        )

    # output number of donors and visits per donor
    if verbose == 1:
        unique_donors = fractions_df.index.get_level_values(donor_obs_column).unique()
        num_unique_donors = len(unique_donors)
        print(f"We detected {num_unique_donors} unique donors.", flush=True)
        if longitudinal_obs_column is not None:
            #Count how many visits each donor has (using the longitudinal column)
            visits_per_donor = adata.obs.groupby(donor_obs_column)[longitudinal_obs_column].nunique()
        
            # Count how many donors have 1 visit, 2 visits, etc.
            donor_visit_counts = visits_per_donor.value_counts().sort_index()
        
            # Calculate the total number of donor visits (sum of all visit counts)
            total_visits = visits_per_donor.sum()
        
            for visit_count, num_donors in donor_visit_counts.items():
                print(f"{num_donors} donors had {visit_count} sample{'s' if visit_count > 1 else ''}.", flush=True)
            print(f"For a total of {total_visits} donor sample{'s' if total_visits > 1 else ''}.", flush=True)
        else:
            print("You did not provide visit information for multiple samples (per donor). If any donor has multiple samples, you must provide the longitudinal_obs_column in the cmverify 'predict' function call or predictions will be inaccurate.")

    
    # Capture the donor IDs before resetting the index
    fractions_df = fractions_df[(fractions_df.sum(axis=1) > 0)]
    donor_ids_partial = fractions_df.index.get_level_values(donor_obs_column).tolist()
    if longitudinal_obs_column is not None:
        visit_ids = fractions_df.index.get_level_values(longitudinal_obs_column).tolist()
        donor_ids = list(zip(donor_ids_partial, visit_ids))
    else:
        donor_ids = [(donor_id, "Baseline") for donor_id in donor_ids_partial]
    
    # Normalize the values to get fractions
    fractions_df = fractions_df.div(fractions_df.sum(axis=1), axis=0).reset_index()
    # Check for missing columns
    fractions_df = check_for_missing_columns(fractions_df)
    # Return the calculated fractions
    return fractions_df, donor_ids


def check_for_missing_columns(fractions_df):
    """
    Ensures that all expected cell type columns are present in the input dataframe.
    
    If any expected columns are missing, they will be added and initialized with zeroes.
    This typically happens with shallow sequencing or small input datasets.

    Parameters:
    -----------
    fractions_df : pandas.DataFrame
        A dataframe containing cell type fractions as columns.

    Returns:
    --------
    pandas.DataFrame
        The updated dataframe with all EXPECTED_COLUMNS present and in the correct order.
    """
    # Ensure that all expected columns exist in the fractions dataframe.
    # If any columns are missing, they will be initialized with zeroes.
    existing_columns = fractions_df.columns.tolist()
    missing_columns = []

    # Iterate over the expected columns and add missing ones with zeroes
    for column in EXPECTED_COLUMNS:
        if column not in existing_columns:
            fractions_df[column] = 0
            missing_columns.append(column)
            
    # If there are missing columns, issue a warning
    if missing_columns:
        print(f"Note: the following cell types were not detected and have been initialized with zeroes: {', '.join(missing_columns)}", flush=True)
        print("Typically this occurs when using the model with fewer cells or sequencing at lower depth. Results may or may not be affected.", flush=True)

    # Ensure the columns are in the expected order
    fractions_df = fractions_df[EXPECTED_COLUMNS]
    fractions_df.index.name = None
    return fractions_df