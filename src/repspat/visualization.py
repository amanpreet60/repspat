import pandas as pd
import networkx as nx

import numpy as np
import matplotlib.pyplot as plt

def plot_spatial_clusters(x, y, labels, figsize=(4, 4), point_size=10, alpha=1.0,
                          title="Spatial Plot of Cells with Cluster Colors", show_legend=True):
    x, y, labels = np.asarray(x), np.asarray(y), np.asarray(labels)
    fig, ax = plt.subplots(figsize=figsize)

    for cid in np.unique(labels):
        m = labels == cid
        ax.scatter(x[m], y[m], s=point_size, alpha=alpha, label=f"Cluster {cid + 1}")

    ax.set(xlabel="X coordinate", ylabel="Y coordinate", title=title)
    ax.grid(False)

    if show_legend:
        ax.legend(title="Clusters", markerscale=1.5,
                  bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)

    return fig, ax


def pairwise_results_to_matrix(df, plot=True):
    # Link clusters that are NOT significantly different (adj_p >= 0.05 = similar spatial distributions)
    df = df.copy()
    df['link'] = (df['adj_p'] >= 0.05).astype(int)
    
    # Nodes and edges
    all_nodes = pd.unique(df[['region_1','region_2']].values.ravel())
    edges = df[df['link'] == 1][['region_1','region_2','obs_mmd_sq']]
    
    # Build graph
    G = nx.Graph()
    G.add_nodes_from(all_nodes)
    for _, row in edges.iterrows():
        G.add_edge(row['region_1'], row['region_2'], weight=row['obs_mmd_sq'])
    
    # Adjacency matrix
    result_matrix = nx.to_pandas_adjacency(G, weight='weight')
    
    if plot:
        fig, ax = plt.subplots()
        pos = nx.circular_layout(G)  # nodes in circle
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw(G, pos, with_labels=True, node_color='white', edge_color='lightgrey',
                width=2, node_size=800, font_size=12, ax=ax)
        nx.draw_networkx_edge_labels(G, pos, edge_labels={k: round(v,2) for k,v in edge_labels.items()})
        plt.show()

    return result_matrix