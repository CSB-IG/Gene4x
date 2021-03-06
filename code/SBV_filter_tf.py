
import networkx as nx
#import numpy as np
from scipy import integrate

def extract_backbone(g, alpha):
  backbone_graph = nx.Graph()
  for node in g:
      k_n = len(g[node])
      if k_n > 1:
          sum_w = sum( g[node][neighbor]['weight'] for neighbor in g[node] )
          for neighbor in g[node]:
              edgeWeight = g[node][neighbor]['weight']
              pij = float(edgeWeight)/sum_w
              f = lambda x: (1-x)**(k_n-2) 
              alpha_ij =  1 - (k_n-1)*integrate.quad(f, 0, pij)[0] 
              if alpha_ij < alpha: 
                  backbone_graph.add_edge( node,neighbor, weight = edgeWeight)
  #print len(backbone_graph)
  return backbone_graph



import networkx as nx
fh=open("./temp/net_tf_unique_nodes_int.txt", 'r')
# ------------------------|----- note 'r' not 'rb'
G=nx.read_weighted_edgelist(fh, nodetype=str)
fh.close()

vet=[0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5]

for i,alpha in enumerate(vet):
  nx.write_edgelist(extract_backbone(G, alpha), "./temp/test.edgelist_tf" + str(i) + ".txt")
  #print(extract_backbone(G, alpha))
