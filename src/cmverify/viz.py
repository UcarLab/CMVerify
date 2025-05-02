# src/cmverify/utils.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc, roc_auc_score

def visualize(results, visit_order=None,figWidth=8,figHeight=3,  dpi_param=100,save=False,filename='cmverify.png', metrics=False):
    """
    Plot longitudinal model-predicted probabilities per donor over timepoints (visits).
    
    This function generates a line plot showing predicted probabilities for each donor
    across multiple visits, with optional labeling of true outcomes and prediction accuracy.

    Parameters:
    - results (list of dict): List where each dictionary contains prediction info with keys:
        - 'donor_id_timepoint': tuple (donor_id, visit label)
        - 'probability': float, model-predicted probability for the donor at that timepoint
        - Optional: 'true_label': ground truth label (e.g., 0/1)
                     'prediction': binary prediction from model
    - visit_order (list, optional): Custom list specifying the order of visit labels on the x-axis.
    - figWidth (float): Width of the figure in inches. Default is 8.
    - figHeight (float): Height of the figure in inches. Default is 3.
    - dpi_param (int): Dots per inch (resolution) of the figure. Default is 100.
    - save (bool): If True, saves the figure to `filename`. Default is False.
    - filename (str): Name of the output image file if `save` is True. Default is 'cmverify.png'.
    - metrics (bool): If True, overlays accuracy or other metrics if available in `results`.

    Returns:
    - None. Displays and/or saves the plot.
    """
    print("Generating visualization", flush=True)
        
    # Convert list of dicts to DataFrame
    df = pd.DataFrame(results)

    # Split donor_id_timepoint tuple into separate columns
    df['Donor_id'] = df['donor_id_timepoint'].apply(lambda x: x[0])
    df['Visit'] = df['donor_id_timepoint'].apply(lambda x: x[1])

    # Apply categorical ordering to visits
    if visit_order is None:
        df['Visit'] = pd.Categorical(df['Visit'], categories=df['Visit'].unique(), ordered=True)
    else:
        df['Visit'] = pd.Categorical(df['Visit'], categories=visit_order, ordered=True)

    # Set up figure and plot base points as a stripplot
    plt.figure(figsize=(figWidth, figHeight), dpi=dpi_param)

    # Plot stripplot, using hue if true_label is available
    if 'true_label' in df.columns:
        sns.stripplot(
            data=df,
            x='Visit',
            y='probability',
            hue='true_label',
            palette=["#1eb8d4", "#faa31b"],
            alpha=1,
            jitter=False,
            edgecolor='black',
            linewidth=0.1
        )
    else:
        sns.stripplot(
            data=df,
            x='Visit',
            y='probability',
            alpha=1,
            jitter=False,
            edgecolor='black',
            linewidth=0.1,
            palette=["black"]
        )

    # Draw dashed lines connecting points for each donor
    for donor_id in df['Donor_id'].unique():
        donor_data = df[df['Donor_id'] == donor_id]
        
        # Drop any rows with missing data
        donor_data = donor_data.dropna(subset=['Visit', 'probability'])
        
        # Only connect dots if donor has multiple timepoints
        if len(donor_data) > 1:
            sorted_cat = (donor_data['Visit'].cat.codes).sort_values()
            plt.plot(sorted_cat, 
                     donor_data['probability'].loc[sorted_cat.index], 
                     linestyle='--', 
                     linewidth=0.2, 
                     color='black', 
                     alpha=0.5,
                     marker=''
                    )

            # Add donor ID text at final timepoint
            last_x = sorted_cat.index[-1]
            last_y = donor_data['probability'].loc[sorted_cat.index[-1]]
            plt.text(
                sorted_cat.iloc[-1]+.1, 
                last_y, str(donor_id), 
                fontsize=6, 
                verticalalignment='center', 
                horizontalalignment='left'
            )
        else:
            plt.text(
                .1, 
                donor_data['probability'], str(donor_id), 
                fontsize=6, 
                verticalalignment='center', 
                horizontalalignment='left'
            )

    # Add axis labels and formatting
    plt.xlabel('Timepoint')
    plt.ylabel('Model Prediction')
    plt.xticks(fontsize=6)
    
    # Draw horizontal threshold line at 0.5
    threshold_line = plt.axhline(y=0.5, color='red', lw=0.5, linestyle='--')
    
    # Custom legend handling
    handles, labels = plt.gca().get_legend_handles_labels()

    # Replace '0' and '1' with your custom labels
    label_map = {'0': 'True CMV-', '1': 'True CMV+'}
    updated_labels = [label_map.get(label, label) for label in labels]

    
    # Add threshold line label
    handles.append(threshold_line)
    updated_labels.append('Decision Threshold')
    
    plt.legend(handles=handles, labels=updated_labels, loc='best', fontsize=8)
    
    # Fit layout and optionally save
    plt.tight_layout()
    if save:
        plt.savefig('scatter_' + filename, dpi=dpi_param, bbox_inches='tight')
    plt.show()

    if ('true_label' in df.columns) and metrics:
        # Print classification report
        print(classification_report(df['true_label'], df['prediction']))
        # show confusion matrix
        cm = confusion_matrix(df['true_label'], df['prediction'])
        # Create subplots
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=dpi_param)
        
        # --- Confusion Matrix ---
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['0', '1'], yticklabels=['0', '1'], ax=axes[0])
        axes[0].set_xlabel('Predicted')
        axes[0].set_ylabel('Actual')
        axes[0].set_title('Confusion Matrix')
        
        # --- ROC Curve ---
        fpr, tpr, thresholds = roc_curve(df['true_label'], df['probability'])
        auc = roc_auc_score(df['true_label'], df['probability'])
        
        axes[1].plot(fpr, tpr, label=f'ROC Curve (AUC = {auc:.2f})', color='blue')
        axes[1].plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random Guess')
        axes[1].set_xlabel('False Positive Rate (FPR)')
        axes[1].set_ylabel('True Positive Rate (TPR)')
        axes[1].set_title('ROC Curve')
        axes[1].legend(loc='lower right')
        
        # Layout and save/show
        plt.tight_layout()
        if save:
            plt.savefig('metrics_' + filename, dpi=dpi_param, bbox_inches='tight')
        plt.show()
