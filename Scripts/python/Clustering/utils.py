import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches

def nan_to_max(l): 
    l[np.isnan(l)] = np.nanmax(l)
    return l

def labels_to_color(labels):
    '''Turn a list of labels into a list of colors.
    
    Parameters
    ----------
    labels : list
        List of labels
    
    Returns
    -------
    col : list 
        List of colors
    legend : list
        list of color patches for each label to plot the legend using plt.legend(handles=legend)
    
    '''
    unique_labels = np.unique(labels)
    colors = sns.color_palette('hls', len(unique_labels))

    # Associate each sample with a color indicating its type
    idx = np.zeros(len(labels)).astype(int)
    legend = []
    for i in range(len(unique_labels)):
        idx[labels == unique_labels[i]] = i
        legend.append(mpatches.Patch(color=colors[i], label=unique_labels[i]))
                       
    col = [colors[i] for i in idx]
    return col, legend