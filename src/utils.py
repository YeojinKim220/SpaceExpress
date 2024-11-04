import numpy as np
import igraph as ig
from tqdm import tqdm
from sklearn.neighbors import kneighbors_graph
import pickle
import matplotlib.pyplot as plt
import os 

def make_complete_graph(pos, k):
    """
    Make a complete graph using the k-nearest neighborhood graph.
    
    Parameters:
    pos (numpy.ndarray): 2D array of node positions
    k (int): Number of nearest neighbors
    
    Returns:
    igraph.Graph: Complete graph
    """
    
    # Check number of connected components
    A = kneighbors_graph(pos, k, mode='connectivity', include_self=False)
    G = ig.Graph.Adjacency((A > 0).toarray().tolist())
    G = G.as_undirected()
    print(f'Number of connected components: {len(G.components())}...')
    
    plt.figure(figsize=(6, 6))
    plt.scatter(pos[:, 0], pos[:, 1], color='black', s=1) 

    components = G.clusters()

    # Connect the closest nodes from different components
    for i in range(len(components)):
        min_dist_idx_list = []
        for j in [x for x in range(len(components)) if x != i]:
            nodes_group1 = components[i]
            nodes_group2 = components[j]
            
            coords_group1 = pos[nodes_group1]
            coords_group2 = pos[nodes_group2]
            
            diff_matrix = np.expand_dims(coords_group1, axis=1) - np.expand_dims(coords_group2, axis=0)
            distances_matrix = np.sqrt(np.sum(np.square(diff_matrix), axis=-1))
            min_dist_idx = np.unravel_index(np.argmin(distances_matrix), distances_matrix.shape)
            min_dist_idx = (nodes_group1[min_dist_idx[0]], nodes_group2[min_dist_idx[1]], np.min(distances_matrix))
            min_dist_idx_list.append(min_dist_idx)
        
        idx = min(min_dist_idx_list, key=lambda x: x[2])
        G.add_edges([(idx[0], idx[1])])  # igraph uses add_edges for adding multiple edges
        plt.plot([pos[idx[0], 0], pos[idx[1], 0]], [pos[idx[0], 1], pos[idx[1], 1]], 'red')
        print(f"Connected {idx[0]} to {idx[1]}")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    return G

def choose_k (adata_list):
    """
    Choose the smallest k that makes the graph connected.

    Parameters:
        adata_list (list): List of anndata.AnnData objects 
    """
    pos_list = [adata.obsm['spatial'] for adata in adata_list]
    
    k = 0
    for pos in pos_list:
        n_cc = float('inf')
        # Find the smallest k that makes the graph connected
        for i in range(4, 21):
            A = kneighbors_graph(pos, i, mode='connectivity', include_self=False)
            G = ig.Graph.Adjacency(A.toarray().tolist())
            G = G.as_undirected()
            if len(G.components()) < n_cc:
                n_cc = len(G.components())
                recommended_k = i
            if G.is_connected():
                print(f'Construct the k-nearest neighborhood graph using k = {i}...')
                break
        k = max(k, recommended_k)
        if k == 20:
            print(f'Graph is disconnected. Recommended k = {recommended_k}...')
    print(f'Choose k = {k} for all datasets...')
    return recommended_k


def shortest_path (pos, out_dir, k = None):
    """
    Calculate the shortest paths between all pairs of nodes in the graph.
    
    Parameters:
    pos (numpy.ndarray): 2D array of node positions
    out_dir (str): Output directory to save the shortest paths
    
    Returns:
    None
    """
    
    if os.path.exists(out_dir):
        print(f'Shortest paths already exist at {out_dir}...')
        return
    
    if k is not None:
        A = kneighbors_graph(pos, k, mode='connectivity', include_self=False)
        G = ig.Graph.Adjacency(A.toarray().tolist())
        G = G.as_undirected()
        
        if G.is_connected() == False:
            print('Graph is not able to connected...')
            print('Making a complete graph...')
            G = make_complete_graph(pos, k)
    
    else:
        n_cc = float('inf')
        # Find the smallest k that makes the graph connected
        for i in range(4, 21):
            A = kneighbors_graph(pos, i, mode='connectivity', include_self=False)
            G = ig.Graph.Adjacency(A.toarray().tolist())
            G = G.as_undirected()
            if len(G.components()) < n_cc:
                n_cc = len(G.components())
                recommended_k = i
            if G.is_connected():
                print(f'Construct the k-nearest neighborhood graph using k = {i}...')
                break
            
        # If the graph is not connected, make a complete graph
        if i == 20:
            print('Graph is not able to connected...')
            print('Making a complete graph...')
            G = make_complete_graph(pos, recommended_k)
            
    # If the graph is connected, save the shortest paths
    with open(out_dir, 'wb') as handle:
        for node in tqdm(range(G.vcount())):
            shortest_paths = G.shortest_paths(source=node, mode=ig.OUT)
            pickle.dump({node: shortest_paths}, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Number of nodes:', G.vcount()), print('Number of edges:', G.ecount())
    print(f'Shortest paths saved to {out_dir}...')
    

def plot_embedding(adata):
    """
    Plots each column of the given matrix as a separate scatter plot in multiple subplots arranged in 2 rows.
    The figure size adjusts based on the number of columns.
    
    Parameters:
    adata (anndata.AnnData): Anndata object containing the SpaceExpress matrix and spatial coordinates
    
    Returns:
    None
    """
    
    matrix = adata.obsm['SpaceExpress']
    pos = adata.obsm['spatial']
    num_dim = matrix.shape[1]
    cols_per_row = (num_dim + 1) // 2

    # Adjust figure size dynamically based on the number of columns
    fig_width = max(5 * cols_per_row, 12)  # Minimum width of 12 inches
    fig_height = 8  # Constant height of 8 inches, adjust as needed

    plt.figure(figsize=(fig_width, fig_height))  # Set the dynamically calculated figure size

    # Iterate over each column in the matrix
    for i in range(num_dim):
        plt.subplot(2, cols_per_row, i + 1)  # Create a subplot for each dimension
        plt.scatter(pos[:,0], pos[:,1], c = matrix[:,i], s = 3)
        plt.title(f'Dimension {i+1}')
        plt.xlabel('Y')
        plt.colorbar()
        plt.ylabel('X')

    plt.tight_layout()
    plt.show()

def jaccard_similarity(adata, k=10):
    """
    Compute the Jaccard similarity between the k-nearest neighborhood graphs constructed using the SpaceExpress embeddings and spatial coordinates.
    
    Parameters:
    adata (anndata.AnnData): Anndata object containing the SpaceExpress matrix and spatial coordinates
    k (int): Number of nearest neighbors
    
    Returns:
    float: Jaccard similarity
    """
    
    emb = adata.obsm['SpaceExpress']
    pos = adata.obsm['spatial']
    # Create adjacency matrices from embeddings and positions using kneighbors_graph
    A = kneighbors_graph(emb, k, mode='connectivity', include_self=False, metric='minkowski').toarray()
    B = kneighbors_graph(pos, k, mode='connectivity', include_self=False, metric='minkowski').toarray()

    # Create igraph graphs from adjacency matrices
    G = ig.Graph.Adjacency((A > 0).tolist())
    G = G.as_undirected()
    H = ig.Graph.Adjacency((B > 0).tolist())
    H = H.as_undirected()

    # Compute edge sets
    G_edges = set(G.get_edgelist())
    H_edges = set(H.get_edgelist())

    # Compute unions and intersections of edges
    union_edges = G_edges.union(H_edges)
    intersection_edges = G_edges.intersection(H_edges)

    # Calculate Jaccard similarity
    jaccard_index = len(intersection_edges) / len(union_edges) if len(union_edges) > 0 else 0

    return round(jaccard_index, 4)

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from pygam import LinearGAM, s
import matplotlib.pyplot as plt    
from sklearn.neighbors import kneighbors_graph

def gam (X, y):
    gam = LinearGAM(s(0, n_splines = 15))
    gam.fit(X, y)

    # Predict on new data
    X_pred = np.linspace(np.min(X), np.max(X), 500).reshape(-1, 1)
    y_pred = gam.predict(X_pred)
    return X_pred, y_pred

def preprocess(exp_1, exp_2):
    exp = np.concatenate([exp_1, exp_2])
    mean_val = np.mean(exp)
    std_val = np.std(exp)

    rm_idx_1 = np.where(exp_1 >= mean_val + 4 * std_val)[0]
    rm_idx_2 = np.where(exp_2 >= mean_val + 4 * std_val)[0]

    if rm_idx_1.size != 0:
        exp_1 = np.delete(exp_1, rm_idx_1, axis=0)
    if rm_idx_2.size != 0:
        exp_2 = np.delete(exp_2, rm_idx_2, axis=0)   
    std_val_1 = np.std(exp_1)
    std_val_2 = np.std(exp_2)

    if std_val_1 != 0:  # Avoid division by zero
        exp_1 = exp_1 / std_val_1
    else:
        print('std == 0')
        
    if std_val_2 != 0:  # Avoid division by zero
        exp_2 = exp_2 / std_val_2
    else:
        print('std == 0')        
        
    return exp_1, exp_2, rm_idx_1, rm_idx_2 

def plot_scatter(fig, ax, loc, color, cmap, title, stream=False, vmin = None, vmax = None, invert = False, color_bar =True):
    plt.rc('font', size=20)

    if cmap == 'cm':
        cmap = LinearSegmentedColormap.from_list('half_coolwarm', colors)
    scatter = ax.scatter(loc[:, 0], loc[:, 1], c=color, cmap=cmap, s=5, alpha=1, vmin = vmin, vmax = vmax)
    ax.set_title(title)
    ax.set_xticks([])    
    ax.set_yticks([])
    
    ax.set_xticks([])    
    ax.set_yticks([])    
    
    if color_bar == True:
        fig.colorbar(scatter, ax=ax)
    if invert == True:
        ax.invert_yaxis()
        
def get_range(data, mode):
    data = np.concatenate(data)
    if mode == 'quantile':
        q1 = np.quantile(data, 0.01)
        q3 = np.quantile(data, 0.99)
    elif mode == 'max':
        q1 = np.min(data)
        q3 = np.max(data)    

    return (q1, q3)  
        
def rm_idx (data, idx):
    for i in range(len(data)):
        data[i] = np.delete(data[i], idx, axis=0)
    return data

def plot_DSE (adata_list, df_fdr, file_list, gene_name, dimension):
    fdr = df_fdr.iloc[dimension,][gene_name]
    dim = dimension

    plt.rc('font', size=20)
    width = (len(adata_list) + 1) * 8
    fig, axs = plt.subplots(2, len(adata_list) + 1, figsize=(width, 11))
    
    coolwarm = plt.colormaps.get_cmap('coolwarm')
    n_colors = coolwarm.N // 2
    colors = coolwarm(np.linspace(0.5, 1, n_colors))
    cm = LinearSegmentedColormap.from_list('half_coolwarm', colors)

    loc_1, emb_1, exp_1 = adata_list[0].obsm['spatial'], adata_list[0].obsm['SpaceExpress'][:, dim], adata_list[0][:, gene_name].X
    loc_2, emb_2, exp_2 = adata_list[1].obsm['spatial'], adata_list[1].obsm['SpaceExpress'][:, dim], adata_list[1][:, gene_name].X

    gene_idx = np.where(adata_list[0].var_names == gene_name)[0][0]
    pred_1, inter_1 = adata_list[0].obsm['DSE-pred'][:, gene_idx, dim], adata_list[0].obsm['DSE-inter'][:, gene_idx, dim]
    pred_2, inter_2 = adata_list[1].obsm['DSE-pred'][:, gene_idx, dim], adata_list[1].obsm['DSE-inter'][:, gene_idx, dim]
    
    exp_1, exp_2, rm_idx_1, rm_idx_2 = preprocess (exp_1, exp_2)
    
    [loc_1, emb_1, pred_1, inter_1] = rm_idx([loc_1, emb_1, pred_1, inter_1], rm_idx_1)
    [loc_2, emb_2, pred_2, inter_2] = rm_idx([loc_2, emb_2, pred_2, inter_2], rm_idx_2)
    
    vmin, vmax = get_range([exp_1, exp_2], 'quantile')

    plot_scatter(fig, axs[1, 0], loc_1, exp_1, 'cm', f'({file_list[0]}) Expression of {gene_name}', vmin = vmin, vmax = vmax)
    plot_scatter(fig, axs[1, 1], loc_2, exp_2, 'cm', f'({file_list[1]}) Expression of {gene_name}', vmin = vmin, vmax = vmax)

    vmin, vmax = get_range([emb_1, emb_2], 'max')

    plot_scatter(fig, axs[0, 0], loc_1, emb_1, 'gist_rainbow', f'({file_list[0]}) Dimension {dim+1}', vmin = vmin, vmax = vmax)
    plot_scatter(fig, axs[0, 1], loc_2, emb_2, 'gist_rainbow', f'({file_list[1]}) Dimension {dim+1}', vmin = vmin, vmax = vmax)
    axs[0, 2].text(0.5, 0.5, f'Gene: {gene_name}\nDimension: {dimension}\nFDR = {fdr:.2e}', ha='center', va='center', fontsize=24, transform=axs[0, 2].transAxes)
    axs[0, 2].set_xticks([])    
    axs[0, 2].set_yticks([])  
    
    x1, y1 = gam(emb_1, pred_1)
    x2, y2 = gam(emb_2, pred_2)
    
    axs[1,2].plot(x1, y1, label='Naive', color='green', linewidth=4)
    axs[1,2].plot(x2, y2, label='Agg. to Pup', color='orange', linewidth=4)
    axs[1,2].set_ylabel(f'Exp. of ${gene_name}$')
    axs[1,2].set_xlabel(f'Dimension {dim+1}')
    axs[1,2].xaxis.set_label_coords(0.5, -0.2)  
    axs[1,2].yaxis.set_label_coords(-0.12, 0.5) 
    axs[1,2].legend(loc='best')

    norm = plt.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap='gist_rainbow', norm=norm)
    sm.set_array([])
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axs[1,2])
    cax = divider.append_axes("bottom", size="4%", pad=0.5)
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_ticks([])  
    cbar.set_label('') 

    sm2 = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm2.set_array([])

    cax2 = divider.append_axes("left", size="2%", pad=0.70)
    cbar2 = fig.colorbar(sm2, cax=cax2, orientation='vertical')
    cbar2.set_ticks([])  
    cbar2.set_label('') 

    plt.tight_layout()
    return fig