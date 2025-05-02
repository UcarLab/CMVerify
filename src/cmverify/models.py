# src/cmverify/models.py
import os
import joblib

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
        'rf_best_estimator': 'models/AIFI_rf_best_estimator.pkl',
        'rf_scaler': 'models/rf_scaler.pkl'
    }

    model_file = model_files.get(model_name)
    
    if model_file is None:
        raise ValueError(f"Model {model_name} not found.")
    
    # Resolve the full path to the model file
    model_path = os.path.join(os.path.dirname(__file__), model_file)
    
    # Check if the model file exists
    if not os.path.exists(model_path):
        raise ValueError(f"Model file {model_path} not found.")
    
    # Load the model using joblib
    model = joblib.load(model_path)
    
    return model
