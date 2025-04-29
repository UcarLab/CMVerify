# src/cmverify/models.py
import pkg_resources
import pickle
import os

def load_model(model_name):
    """
    Load a pre-trained model from the package's resources.

    Parameters
    ----------
    model_name : str
        The name of the model to load.

    Returns
    -------
    model : object
        The loaded model.
    """
    model_files = {
        'AIFI_L3': 'models/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19_jl.pkl'
    }

    model_file = model_files.get(model_name)
    
    if model_file is None:
        raise ValueError(f"Model {model_name} not found.")
    
    # Resolve the full path to the model file
    model_path = os.path.join(os.path.dirname(__file__), model_file)
    
    # Check if the model file exists
    if not os.path.exists(model_path):
        raise ValueError(f"Model file {model_path} not found.")
    
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    return model

