# src/cmverify/utils.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

def visualize(results, visit_order=None,figWidth=8,figHeight=3,  dpi_param=100,save=False,filename='cmverify_viz.png',show=True):
    """
    Plot longitudinal prediction probabilities per donor over visits.

    Parameters:
        results (list of dict): Each dict must include:
            - 'donor_id_timepoint': a tuple (donor_id, visit)
            - 'probability': model-predicted probability
        visit_order (list, optional): Custom ordering of visit labels.
        figWidth, figHeight (float): Figure dimensions.
        dpi_param (int): DPI for the figure.
        verbose (int): Verbosity level.
        save (bool): Whether to save the figure.
        filename (str): Output filename if saving.
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
        plt.legend(title='True CMV', labels=['CMV-', 'CMV+'], fontsize=8, title_fontsize=9, loc='best')
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

    # Add axis labels and formatting
    plt.xlabel('Timepoint')
    plt.ylabel('Model Prediction')
    plt.xticks(fontsize=6)
    
    # Draw horizontal threshold line at 0.5
    threshold_line = plt.axhline(y=0.5, color='red', lw=0.5, linestyle='--')
    
    # Fit layout and optionally save
    # Custom legend handling
    handles, labels = plt.gca().get_legend_handles_labels()
    
    # If 'true_label' in df, add color patches
    if 'true_label' in df.columns:
        custom_patches = [
            mpatches.Patch(color="#1eb8d4", label='True CMV-'),
            mpatches.Patch(color="#faa31b", label='True CMV+')
        ]
        handles.extend(custom_patches)
    
    # Add threshold line label
    handles.append(threshold_line)
    labels.append('Decision Threshold')
    
    plt.legend(handles=handles, labels=labels, loc='best', fontsize=8)

    plt.tight_layout()
    if save:
        plt.savefig(filename, dpi=dpi_param, bbox_inches='tight')
        
    if show:
        plt.show()
