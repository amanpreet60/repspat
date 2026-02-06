import pandas as pd
import networkx as nx


def pairwise_results_to_matrix(df, plot=True):
    # Link if adjusted p-value >= 0.05
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
        pos = nx.circular_layout(G)  # nodes in circle
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw(G, pos, with_labels=True, node_color='white', edge_color='lightgrey',
                width=2, node_size=800, font_size=12)
        nx.draw_networkx_edge_labels(G, pos, edge_labels={k: round(v,2) for k,v in edge_labels.items()})
    
    return result_matrix