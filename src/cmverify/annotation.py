# src/CMVerify/annotation.py
import celltypist
import pandas as pd
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
    dict: Dictionary containing annotated labels and scores.
    """

    label_results = {}

    # Define the columns for predictions and scores
    label_column = f'{model_name}_prediction'

    #resolve the path to the model
    model_file = 'models/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl'
    model_path = os.path.join(os.path.dirname(__file__), model_file)
    
    # Annotate the AnnData object using the selected model
    predictions = celltypist.annotate(adata, model=model_path)

    # Extract the predicted labels and rename the column
    labels = predictions.predicted_labels
    labels = labels.rename({'predicted_labels': label_column}, axis=1)

    # Store the results in the label_results dictionary
    label_results[model_name] = labels

    return label_results

def check_and_add_labels(adata, label_results, model_name):
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
    label_df = label_results[model_name]  # Since we have only one model

    label_column = f'{model_name}_prediction'
    
    # Print a summary of the predicted labels
    label_summary = label_df[label_column].value_counts()

    # Check if there are at least 20 unique cell types in the predictions
    unique_cell_types = label_summary.shape[0]
    if unique_cell_types < 20:
        warnings.warn(f"Warning: The model '{model_name}' predicted only {unique_cell_types} unique cell types. You may not have enough cells.", stacklevel=2)
    else:
        print("Top 20 most frequent cell types detected with cell counts:")
        print(label_summary[0:20], end = '\n\n')

    # Add the predicted labels and score to adata.obs
    adata.obs[label_column] = label_df[label_column]

def calculate_cell_type_fractions(adata, model_name, donor_obs_column, longitudinal_obs_column=None):
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
        # Capture the donor IDs before resetting the index
        fractions_df = fractions_df[(fractions_df.sum(axis=1) > 0)]
        donor_ids_partial = fractions_df.index.get_level_values(donor_obs_column).tolist()
        visit_ids = fractions_df.index.get_level_values(longitudinal_obs_column).tolist()
        donor_ids = list(zip(donor_ids_partial, visit_ids))
    else:
        obs_df = adata.obs[[donor_obs_column, label_column]]
        # Calculate the fraction of cells for each label per patient
        fractions_df = (
            obs_df.groupby([donor_obs_column, label_column])
            .size()
            .unstack(fill_value=0)  # Converts to a wide format with labels as columns
        )
        fractions_df = fractions_df[(fractions_df.sum(axis=1) > 0)]
        # Capture the donor IDs before resetting the index
        donor_ids = fractions_df.index.get_level_values(donor_obs_column).tolist()
    
    # Normalize the values to get fractions
    fractions_df = fractions_df.div(fractions_df.sum(axis=1), axis=0).reset_index()

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
        print(f"The following cell types were not detected and have been initialized with zeroes: {', '.join(missing_columns)}")

    # Ensure the columns are in the expected order
    fractions_df = fractions_df[EXPECTED_COLUMNS]
    
    # Return the calculated fractions
    return fractions_df, donor_ids
