import pickle, random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.neighbors import kneighbors_graph
from tqdm import tqdm
import copy
import scipy
from .spaceexpress_dse import *

def get_index (data):
    n_cells = data.shape[0]
    mask = torch.triu(torch.ones((n_cells, n_cells)), diagonal=1)        
    index = torch.nonzero(mask, as_tuple=False)
    return index

def loaded_paths(file_path):
    loaded_paths = {}

    with open(file_path, 'rb') as handle:
        while True:
            try:
                data = pickle.load(handle)
                loaded_paths.update(data)
            except EOFError:
                break  
    shortest_path = [loaded_paths[i][0] for i in range(len(loaded_paths))]
    shortest_path = np.array(shortest_path)
    return shortest_path

class SpaceExpress(nn.Module):
    def __init__(self, input_dim, hidden_dim, latent_dim):
        super(SpaceExpress, self).__init__()

        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim), 
            nn.ReLU(),
            nn.Linear(hidden_dim, latent_dim),
        )

    def forward(self, x):
        encoded = self.encoder(x)
        return encoded

class KKLoss:
    def __init__(self, device, shortest_path):
        self.device = device
        self.shortest_path = torch.tensor(shortest_path, dtype=torch.float32).to(device)

    def __call__(self, p, idx):
        d = self.shortest_path
        mask = d != 0
        d_inv = torch.where(mask, 1 / d**2, torch.zeros_like(d))
        p = p.unsqueeze(1) - p.unsqueeze(0)
        p = torch.sum(torch.abs(p),dim=-1)
        diff = p - d
        out = d_inv * mask * diff.pow(2)
        out = out[idx[:, 0], idx[:, 1]]
        out = torch.sum(out)   
        return out

def get_avg_neighbor(pos, count, k):
    A = kneighbors_graph(pos, k, mode='connectivity', include_self=False)
    neighbor = A.dot(count) / k
    neighbor_dense = neighbor.toarray() if hasattr(neighbor, "toarray") else neighbor
    count_dense = count.toarray() if hasattr(count, "toarray") else count
    out = np.concatenate((count_dense, neighbor_dense), axis=1)
    return out

def train_SpaceExpress(adata, shortest_file_path, device = None, epochs = 10000, lr = 0.01, hid_dim = 32, emb_dim = 8, 
                       patience = 100, random_seed = 42, batch_size = 256, num_hvg = 1000, save_model = False):
    """
    Train SpaceExpress model.
    
    Parameters:
    adata (AnnData): Anndata object
    shortest_path_PATH (str): Path to the shortest path file
    device (str): Device to use
    lr (float): Learning rate
    hid_dim (int): Hidden dimension
    emb_dim (int): Embedding dimension
    patience (int): Patience for early stopping
    random_seed (int): Random seed
    batch_size (int): Batch size
    
    Returns:
    AnnData: Anndata object with SpaceExpress embedding
    """    

    # Set random seed
    torch.manual_seed(random_seed)
    torch.cuda.manual_seed(random_seed)
    torch.cuda.manual_seed_all(random_seed) 
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(random_seed)
    random.seed(random_seed)

    # Set device
    if device == 'mps':
        device = torch.device("mps")
    else:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    # Preprocessing
    percentile_95 = [np.percentile(adata.X[:,i], 95) for i in range(adata.shape[1])] 
    adata.X = np.array([np.clip(adata.X[:,i], a_min=None, a_max=percentile_95[i]) for i in range(adata.shape[1])]).T 
    
    is_count_data = np.all(np.equal(np.mod(adata.X, 1), 0))
    if adata.shape[1] > num_hvg:
        if is_count_data == True:
            sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat_v3', subset = True) 
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, subset = True) 
    sc.pp.scale(adata) 

    # Get data
    data = adata.X 
    pos = adata.obsm['spatial']
    shortest_path = loaded_paths(shortest_file_path)      
    
    data = torch.tensor(data, dtype=torch.float).to(device) 
    model = SpaceExpress(data.shape[1], hid_dim, emb_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    index = get_index(data)
    kkloss = KKLoss(device, shortest_path)

    best_train_loss = float('inf')
    patience_counter = 0
    
    # Training loop
    print('Start training...')
    with tqdm(total=epochs) as pbar:
        for epoch in tqdm(range(epochs)):
            total_loss = 0
            
            idx = index[torch.randperm(index.size(0))[:batch_size**2]]
            optimizer.zero_grad()
            embeddings = model(data.to(device))
            loss = kkloss(embeddings, idx)        
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
            pbar.set_description(f"Epoch {epoch+1}/{epochs}")
            pbar.set_postfix(loss=f"{total_loss:.4f}")
            pbar.update(1)  
                
            if total_loss < best_train_loss:
                best_train_loss = total_loss
                patience_counter = 0
                best_model = copy.deepcopy(model)
            else:
                patience_counter += 1

            if patience_counter >= patience:
                break

    # Switch to evaluation mode
    best_model.eval()
    emb = best_model(data)
    emb = emb.cpu().detach().numpy()
    
    if save_model == True:
        return emb, best_model
    
    return emb

def train_SpaceExpress_multi(adata_list_input, shortest_file_path_list, device = None, epochs = 10000, lr = 0.01, hid_dim = 32, 
                             emb_dim = 4, patience = 100, random_seed = 42, batch_size = 256, num_hvg = 1000, save_model = False):
    """
    Train SpaceExpress model.
    
    Parameters:
    adata_list (list): List of AnnData objects
    shortest_path_PATH (str): Path to the shortest path file
    device (str): Device to use
    lr (float): Learning rate
    hid_dim (int): Hidden dimension
    emb_dim (int): Embedding dimension
    patience (int): Patience for early stopping
    random_seed (int): Random seed
    batch_size (int): Batch size
    
    Returns:
    adata_list (list): List of AnnData objects with SpaceExpress embeddings
    """    
    
    # Set random seed
    torch.manual_seed(random_seed)
    torch.cuda.manual_seed(random_seed)
    torch.cuda.manual_seed_all(random_seed) 
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(random_seed)
    random.seed(random_seed)
    
    adata_list = [i.copy() for i in adata_list_input]
    
    # Set device
    if device == 'mps':
        device = torch.device("mps")
    else:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    # Preprocessing
    for i, adata in enumerate(adata_list):
        # Check if adata.X is sparse and convert if necessary
        if isinstance(adata.X, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
            adata.X = adata.X.toarray()    

        percentile_95 = [np.percentile(adata.X[:,i], 95) for i in range(adata.shape[1])] 
        adata.X = np.array([np.clip(adata.X[:,i], a_min=None, a_max=percentile_95[i]) for i in range(adata.shape[1])]).T 
        is_count_data = np.all(np.equal(np.mod(adata.X, 1), 0))
        if adata.shape[1] > num_hvg:
            if is_count_data == True:
                sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat_v3', subset = True) 
            else:
                sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, subset = True) 
        sc.pp.scale(adata)
        adata_list[i] = adata
    
    intersecting_genes = set.intersection(*(set(adata.var_names) for adata in adata_list))
    adata_list = [adata[:, list(intersecting_genes)] for adata in adata_list]
    print(f'Size of the input data: {[adata_list[i].shape for i in range(len(adata_list))]}')
    
    # Get data
    num_data = len(adata_list)
    data_list = [adata.X for adata in adata_list]
    pos_list = [adata.obsm['spatial'] for adata in adata_list]
    shortest_path_list = [loaded_paths(i) for i in shortest_file_path_list]
    
    num_gene = data_list[0].shape[1]
    print('number of genes:', num_gene)
    for i in data_list: 
        assert num_gene == i.shape[1], "The number of genes should be the same across all datasets"
    
    data_list = [torch.tensor(data, dtype=torch.float).to(device) for data in data_list]    
    model = SpaceExpress(num_gene, hid_dim, emb_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
      
    index_list = [get_index(data) for data in data_list]
    kkloss_list = [KKLoss(device, shortest_path) for shortest_path in shortest_path_list]

    best_train_loss = float('inf')
    patience_counter = 0
    
    # Training loop
    print('Start training...')
    with tqdm(total=epochs) as pbar:
        for epoch in tqdm(range(epochs)):
            total_loss = 0
            reg_loss = 0
            
            loss = 0
            for i in range(num_data):
                data = data_list[i]
                index = index_list[i]
                kkloss = kkloss_list[i]
                
                optimizer.zero_grad()
                embeddings = model(data.to(device))
                idx = index[torch.randperm(index.size(0))[:batch_size**2]]
                loss += kkloss(embeddings, idx)      
                
            loss.backward()
            optimizer.step()
            total_loss += loss.item()        

            pbar.set_description(f"Epoch {epoch+1}/{epochs}")
            pbar.set_postfix(loss=f"{total_loss:.4f}")
            pbar.update(1)  # Update the progress bar once per epoch
            
            if total_loss < best_train_loss:
                best_train_loss = total_loss
                patience_counter = 0
                best_model = copy.deepcopy(model)
            else:
                patience_counter += 1

            if patience_counter >= patience:
                break

    # Switch to evaluation mode
    best_model.eval()
    out = []
    for i in range(num_data):
        data = data_list[i]
        emb = best_model(data)
        emb = emb.cpu().detach().numpy()
        out.append(emb)
        
    if save_model == True:
        return out, best_model
    
    return out
