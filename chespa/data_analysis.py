# -*- coding: utf-8 -*-
"""
chespa: Streamlining chemical space evaluation of molecular sets and
assessing trends of ionizability

@author: Jamie Nunez
(C) 2020 - Pacific Northwest National Laboratory
"""

#%% Imports
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial import ConvexHull
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import sys

sns.set_context({'figure.dpi':600})
sns.set(font_scale=0.8)

#%% Globals
METHOD_NAMES = ['(+)-ESI', '(-)-ESI', 'ESI']
PATH = ''
FAX = 12
FLAB = 14
NSTD = 1.95996 # nstd = X standard deviations, 1.95996 = 95% CI

#%% Analysis functions

def normalize(data, avgs_path='avgs.npy',
              stds_path='stds.npy'):
    '''
    Normalize data to have an average of 0 and standard deviation of 1,
    given averages and standard deviations from an existing space.

    Parameters
    ----------
    data: Numpy array (MxN) 
        Matrix with data to normalize
    avg_path: string (optional)
        Path to .npy file with N average values to use
        Default = 'avgs.npy'
    stds: Numpy array or list (N elements) (optional)
        Path to .npy file with N standard deviation values to use.
        Default = 'stds.npy'

    Returns
    ----------
    data: Numpy array (MxN)
        Normalized data
    '''
    avgs = np.load(avgs_path)
    stds = np.load(stds_path)
    data = np.copy(data)

    for j in range(data.shape[1]):
        avg = avgs[j]
        std = stds[j]
        data[:, j] = (data[:, j] - avg) / std

    data[np.where(np.isnan(data))] = 0
        
    return data


def get_pca(data, coeff_path='coeff.npy'):
    '''
    Get PCA transform using coefficients from existing space.

    Parameters
    ----------
    data: Numpy array (MxN) 
        Matrix with data to transform
    coeff_path: string (optional)
        Path to .npy file with NxN coefficients generated from previous PCA.
        Default = 'coeff.npy'

    Returns
    ----------
    data: Numpy array (MxN)
        Transformed data
    '''
    coeff = np.load(coeff_path)
    return np.dot(coeff.T, data.T).T


def get_counts(df, labels, col, uniq_labels=None):
    '''
    Returns the counts of elements in a Pandas DataFrame column. If any
    uniq_labels are not present, adds with a count of 0. If NaN elements
    exists, removes them from the final summary. Sorts counts elements at
    the very end.

    Parameters
    ----------
    df: Pandas DataFrame
        DataFrame to pull from.
    col: string
        Name of column to get counts for
    uniq_labels: iterable (optional)
        Labels to ensure are included in the final count. If any are not,
        they are added to the count summary with a count of 0.

    Returns
    ----------
    res: Pandas Series
        Final count summary
    '''

    # Remove duplicates then join with target
    cols = ['CpdInd', col]
    labels = labels.copy()[cols]
    labels.drop_duplicates(subset=cols, inplace=True)
    df_temp = df.merge(labels.set_index('CpdInd'), on='CpdInd')
    
    # Count
    res = df_temp[col].value_counts()
    
    # Force all labels to be present (if any missing)
    if uniq_labels is not None:
        for label in uniq_labels:
            if label not in res.index:
                res.set_value(label, 0)
    
    # Remove nan labels, if they exist
    res = res.loc[res.index.dropna()]
    
    # Sort
    res = res.sort_index()
    
    return res


def get_method_counts(df, labels, col, uniq_labels=None):
    '''
    Returns count information for each method in the global variable
    METHOD_NAMES. Anywhere the given method is 0 is removed from the DataFrame
    before counting. For more information on how the count itself is performed,
    see get_counts().

    Parameters
    ----------
    df: Pandas DataFrame
        DataFrame to pull from.
    col: string
        Name of column to get counts for
    uniq_labels: iterable (optional)
        Labels to ensure are included in the final count. If any are not,
        they are added to the count summary with a count of 0.

    Returns
    ----------
    res: Pandas Series
        Final count summary
    '''
    results = []
    for method in METHOD_NAMES:
        df_temp = df[df[method] > 0]
        results.append(get_counts(df_temp, labels, col, uniq_labels=uniq_labels))
    return results


#%% Library interaction functions

def get_cpd_labels(fname):
    '''
    Returns group/cluster labels for each compound in the 'ClassyFire',
    'ChemSpace', 'DarkChem', and 'Substructures' tabs of the Excel
    spreadsheet stored at fname. Joins all into one DataFrame.
    '''

    # ClassyFire labels
    cols = ['CpdInd', 'ClassyFire']
    data = pd.read_excel(fname,
                         sheet_name='ClassyFire', usecols=cols, header=0,
                         dtype={cols[1]: str})

    # Add all other Chemical Spaces
    for chemspace in ['ChemSpace', 'DarkChem', 'MACCS', 'SPECTRe']:
        cols = ['CpdInd', chemspace]
        data_temp = pd.read_excel(fname,
                                  sheet_name=chemspace, usecols=cols, header=0,
                                  dtype={cols[1]: str})
        data = data.join(data_temp.set_index('CpdInd'), on='CpdInd')

    return data
    

def get_spiked_results(fname):
    '''
    Returns observed, not observed, and spiked in compound information from
    the Excel spreadsheet stored at fname. Group/cluster labels are added to
    the final DataFrame.
    '''

    # Observed data
    cols = ['CpdInd', 'Mix']
    cols.extend(METHOD_NAMES)
    data_obs = pd.read_excel(fname, sheet_name='Observed', usecols=cols,
                             header=0)

    # Not observed data
    cols = ['CpdInd', 'Mix']
    data_not_obs = pd.read_excel(fname, sheet_name='NotObs', usecols=cols,
                                 header=0)
    
    # Spiked in compounds
    cols = ['CpdInd', 'Mix']
    data_spiked = pd.read_excel(fname, sheet_name='SpikedIn', usecols=cols,
                                header=0)
    
    return data_obs, data_not_obs, data_spiked

#%% Plotting Functions

# Adapted from https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
def silhouette_analysis(data, n_clusters, filename=None):
    '''
    Generates silhouette plot which can be used to assess the distribution
    of elements across the given number of clusters after KMeans is performed.

    Parameters
    ----------
    data: Arrary
        Data to perform clustering on
    n_clusters: int
        Number of clusters to use
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''
    # Create a subplot with 1 row and 2 columns
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(3.5, 3.5)
    ax.set_xlim([-0.1, 1])
    ax.set_ylim([0, len(data) + (n_clusters + 1) * 10])
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(data)
    silhouette_avg = silhouette_score(data, cluster_labels)
    print('For n_clusters =', n_clusters,
          'The average silhouette_score is :', silhouette_avg)
    sys.stdout.flush()
    
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(data, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i + 1))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_xlabel('Coefficient values')
    ax.set_ylabel('Cluster Label')

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color='red', linestyle='--')

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    
    if filename is not None:
        plt.savefig(filename, dpi=600, bbox_inches='tight')
    
    plt.show()


def gen_colors(n):
    '''
    Helper function that returns n evenly spaced colors from the Viridis
    colormap.
    '''
    c = cm.viridis
    return [c(int(255 / n * x)) for x in range(n)]


def _plot(x, y, xlab, ylab, a=1, c='k'):
    '''
    Plot convex hull for the given x and y values.
    '''    
    # Convex hull
    x = np.copy(x)
    y = np.copy(y)
    points = np.vstack((x, y)).T
    hull = ConvexHull(points)
    for simplex in hull.simplices:
        h, = plt.plot(points[simplex, 0], points[simplex, 1], c=c, lw=3)
    plt.fill(points[hull.vertices,0], points[hull.vertices,1], color=c, alpha=a)
    return h


def plot_vars(data, labels, klabels, filename=None):
    '''
    Plots convex hulls for given data, 2 dimensions at a time, separately for
    rows with different associated klabels.

    Parameters
    ----------
    data: Arrary
        Data to plot with convex hulls.
    labels: iterable
        Names for each column in data (used to label axes)
    klabels: iterable
        Group name of each row in data.
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''

    cols = 2
    if data.shape[1] < 4:
        cols = 1
    rows = np.ceil(data.shape[1] / 4)
    clusters = pd.unique(klabels)
    clusters.sort()
    n_clusters = len(clusters)
    c = gen_colors(n_clusters)
    plt.figure(figsize=(4 * cols, 4 * rows))
    plt.subplots_adjust(wspace=0.125, hspace=0.125)

    for j in range(int(np.floor(data.shape[1] / 2))): # If odd num, it'll skip the last one
        num = 2 * j
        plt.subplot(rows * 100 + 20 + j + 1) 
            
        # Plot points
        h = []
        for i in range(n_clusters):
            ind = np.where(np.array(klabels) == clusters[i])
            if len(ind[0]) > 0:
                h.append(_plot(data[ind, num], data[ind, num + 1], labels[num],
                         labels[num + 1], c=c[i], a=0.3))
        if j == 1:
            plt.legend(h, clusters, fontsize=FAX, loc='upper right')
        if data.shape[1] == 2:
            plt.legend(h, clusters, loc='best', fontsize=FAX)

        plt.axis('off')

    if filename is not None:
        plt.savefig(filename, dpi=600, bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_bar_percent(data, ylabel, title=None, palette='viridis',
                     ymin=0, ymax=100, fig_width=1.75, fig_height=1.75,
                     filename=None):
    '''
    Plot histogram for each element in data.

    Parameters
    ----------
    data: Pandas Series
        Data to plot.
    ylabels: string
        Labels for y axis.
    title: string (optional)
        Title to place above plot. If None, not title added. Default = None.
    palette: string (optional)
        Colormap to use for histogram. Options are those included with
        matplotlib. Default = 'viridis'.
    ymin: float (optional)
        Minimum y value
    ymax: float (optional)
        Maximum y value
    fig_width: float (optional)
        Width of figure to generate. Default = 1.75.
    fig_height: float (optional)
        Height of figure to generate. Default = 1.75.
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''

    # Plot
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.barplot(data.index, data.values, palette=palette)
    
    # Format
    ax.set_ylim(ymin, ymax)
    plt.xlabel(data.name)
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
                
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
        plt.show()
    return


def plot_bars_percent_method(datasets, totals, ylabel, names, filename=None):
    '''
    Wrapper for plot_bar_percent() that creates a histogram for each method in
    the global variable METHOD_NAMES for each dataset.

    Parameters
    ----------
    datasets: List of Pandas Series
        List of data to plot
    totals: List of floats
        List of values to divide each dataset by
    ylabel: string
        Label for y axis.
    names: List of strings
        List of substrings to use when saving the figure (if saved).
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''
    for t, l, name in zip(totals, datasets, names):
        for i in range(len(METHOD_NAMES)):
            mname = METHOD_NAMES[i]
            plot_bar_percent(l[i] / t * 100, ylabel, title=mname)
            
            if filename is not None:
                plt.savefig(filename % (name, mname),
                            bbox_inches='tight', dpi=600)
                plt.show()


def plot_bar_counts_multi(datasets, ylabel, titles, ymin=0, ymax=100, 
                          palette='viridis', fig_width=11, fig_height=11, 
                          wspace=0.6, hspace=0.75, rows=5, cols=4, filename=None):
    '''
    Plot histograms for each given dataset.

    Parameters
    ----------
    data: Pandas Series
        Data to plot.
    ylabels: string
        Labels for y axis.
    titles: list of strings
        Titles to place above plots.
    palette: string (optional)
        Colormap to use for histogram. Options are those included with
        matplotlib. Default = 'viridis'.
    ymin: float (optional)
        Minimum y value
    ymax: float (optional)
        Maximum y value
    fig_width: float (optional)
        Width of figure to generate. Default = 1.75.
    fig_height: float (optional)
        Height of figure to generate. Default = 1.75.
    wspace: float (optional)
        Space (width) between subplots. Default = 0.6.
    hspace: float (optional)
        Space (height) between subplots. Default = 0.6.
    rows: int (optional)
        Number of rows to use for subplots. Default = 4.
    cols: int (optional)
        Number of columns to use for subplots. Default = 4.
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''

    # Plot set up
    plt.figure(figsize=(fig_width, fig_height))
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    subplot = 1
    
    # Plot
    for data, title in zip(datasets, titles):
        plt.subplot(rows, cols, subplot)
    
        # Plot
        ax = sns.barplot(data.index, data.values, palette=palette)
    
        # Final formatting
        ax.set_ylim(ymin, ymax)
        plt.xlabel(data.name)
        plt.ylabel(ylabel)
        if title is not None:
            plt.title(title)
        subplot += 1
        
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
        plt.show()
    return


def plot_bar_averages_multi(x, y_list, data, palette='viridis',
                        fig_width=5.5, fig_height=8, wspace=0.4, hspace=0.4,
                        rows=5, cols=2, filename=None):
    '''
    Plot histogram for each element in data.

    Parameters
    ----------
    x: string
        Name of column to set x-axis labels of histogram
    y_list: list of strings
        List of column names to set height of bars
    data: Pandas Dataframe
        DataFrame to pull plottinf information from.
    palette: string (optional)
        Colormap to use for histogram. Options are those included with
        matplotlib. Default = 'viridis'.
    fig_width: float (optional)
        Width of figure to generate. Default = 1.75.
    fig_height: float (optional)
        Height of figure to generate. Default = 1.75.
    wspace: float (optional)
        Space (width) between subplots. Default = 0.6.
    hspace: float (optional)
        Space (height) between subplots. Default = 0.6.
    rows: int (optional)
        Number of rows to use for subplots. Default = 4.
    cols: int (optional)
        Number of columns to use for subplots. Default = 4.
    filename: string (optional)
        Path to save figure at. If None, figure is not saved. Default = None.
    '''

    # Plot set up
    plt.figure(figsize=(fig_width, fig_height))
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    subplot = 1
    
    # Begin loop
    for col in y_list:
        plt.subplot(rows, cols, subplot)
        
        # Plot
        sns.barplot(x=x, y=col, data=data, palette=palette)
        
        # Final formatting
        plt.xlabel(None)
        subplot += 1

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
        plt.show()

    return


def draw_cpds(data, label, n=5, filename=None):
    '''
    Randomly selects n representatives from each group (distinguished by
    label) in data and draws their chemical structure. Assumes data includes
    a 'Mols' column, which is an object generated by RDKit.
    '''
    # Start new figure
    plt.figure(figsize=(8, 11))
    plt.subplots_adjust(hspace=0, wspace=0)
    subplot = 1
    data.sort_values(label, inplace=True)
    groups = data[label].unique()
    for group in groups:

        # Choose <=n randomly selected smiles
        ds_data = data[data[label] == group]
        this_n = min(n, len(ds_data))
        smiles = ds_data.sample(n=this_n, random_state=0)['Mols'].tolist()

        # Plot
        for smi in smiles:
            plt.subplot(len(groups), this_n, subplot)
            img = Chem.Draw.MolToImage(smi, size=(600, 600))
            plt.imshow(img)
            plt.axis('off')
            subplot += 1
        
        subplot += (n - this_n)  # Move to next row

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
        plt.show()

    return
