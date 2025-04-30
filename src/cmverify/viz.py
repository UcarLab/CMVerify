# src/cmverify/utils.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_longitudinal_predictions(results, visit_order=None):
    """
    Plot longitudinal prediction probabilities per donor over visits.

    Parameters:
        results (list of dict): Each dict must include:
            - 'donor_id': a tuple (donor_id, visit)
            - 'probability': model-predicted probability
    """
    # Convert list of dicts to DataFrame
    df = pd.DataFrame(results)

    # Split donor_id tuple
    df['Donor_id'] = df['donor_id'].apply(lambda x: x[0])
    df['Visit'] = df['donor_id'].apply(lambda x: x[1])

    # Ensure Visit is a properly ordered categorical variable
    if visit_order:
        df['Visit'] = pd.Categorical(df['Visit'], categories=visit_order, ordered=True)
    else:
        df['Visit'] = pd.Categorical(df['Visit'], categories=sorted(df['Visit'].unique()), ordered=True)

    # Plotting
    plt.figure(figsize=(8, 3), dpi=100)
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

    # Draw lines per donor
    for donor_id in df['Donor_id'].unique():
        donor_data = df[df['Donor_id'] == donor_id]
        # Filter out rows where 'cohort_timepoint' or 'col' is NaN
        donor_data = donor_data.dropna(subset=['Visit', 'probability'])
        # Plot the line if there is data for at least two timepoints
        if len(donor_data) > 1:
            sorted_cat = (donor_data['Visit'].cat.codes).sort_values()
            c='black'
            width = 0.2
            plt.plot(sorted_cat, donor_data['probability'].loc[sorted_cat.index], 
                     linestyle='--', linewidth=width, color=c, alpha=0.5,marker='')  # Adjust alpha for transparency
            last_x = sorted_cat.index[-1]
            last_y = donor_data['probability'].loc[sorted_cat.index[-1]]
            plt.text(sorted_cat.iloc[-1]+.1, last_y, str(donor_id), 
                    fontsize=6, verticalalignment='center', 
                    horizontalalignment='left')


    # Final formatting
    plt.xlabel('Timepoint')
    plt.ylabel('Model Prediction')
    plt.xticks(fontsize=6)
    plt.axhline(y=0.5, color='red', lw=0.5, linestyle='--')
    plt.tight_layout()
    plt.show()
