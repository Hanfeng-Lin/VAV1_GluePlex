import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import umap
import os
import pandas as pd
import umap.plot

"""
Use with pdbs2csv.py
"""

# os.chdir("C:\\Data\OneDrive - Baylor College of Medicine\\Hanfeng Modeling PC\\VAV1\\pesto\\HADDOCK_defined-active_CRBNautopassive\\599676-CRBN_VAV1_yds_summary\\combined_pdbs\\CRBNaligned_SH3_CA")
coordinates = pd.read_csv("atom.csv")
coordinates_data = coordinates.drop(columns=["ipTM", "file", "cluster", "pair_chains_ipTM_CRBN_VAV1_avg"]).values

reducer = umap.UMAP(n_neighbors=5, metric="euclidean", min_dist=0.7, random_state=42)
embedding = reducer.fit(coordinates_data)

hover_data = pd.DataFrame({'conformation': coordinates["file"].values,
                           'cluster': coordinates["cluster"].values,
                           'ipTM': coordinates["ipTM"].values,
                           'chainpair_ipTM': coordinates["pair_chains_ipTM_CRBN_VAV1_avg"].values})

# Color by conformation

umap.plot.output_file("UMAP_cluster.html")
cluster_plot = umap.plot.interactive(embedding, labels=coordinates["cluster"], hover_data=hover_data,
                          color_key_cmap="Spectral", background="white", point_size=5)

# color_kry_cmap/labels flag for categorical data
umap.plot.show(cluster_plot)

plt.figure(figsize=(6, 6))
ax = plt.gca()
ax = umap.plot.points(embedding, labels="cluster " + coordinates["cluster"].astype(str), color_key_cmap="Spectral", ax=ax)
if len(ax.collections) > 0:
    ax.collections[0].set_sizes([20])
ax.set_title("VAV1 - NGT-201-12 Co-folding from Cluster ", pad=16, fontsize=16)
ax.set_xlabel("UMAP_1", fontsize=14)
ax.set_ylabel("UMAP_2", fontsize=14)
ax.xaxis.set_major_locator(ticker.AutoLocator())
ax.yaxis.set_major_locator(ticker.AutoLocator())
plt.tight_layout()
plt.savefig("UMAP_cluster.png", dpi=300)
plt.show()


# Color by energy
cmap = plt.cm.RdYlBu.reversed()
umap.plot.output_file("UMAP_ipTM.html")
ipTM_plot = umap.plot.interactive(embedding, values=coordinates["ipTM"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)  # cmap/value flag for continuous data

umap.plot.show(ipTM_plot)

embedding = reducer.fit_transform(coordinates_data)

plt.figure(figsize=(7, 6))  
sc = plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=coordinates["ipTM"].values,
    cmap=cmap,
    s=20,
    vmin=0,
    vmax=1.05
)

# Create the colorbar and customize it
cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
cbar.set_label('ipTM', rotation=270, labelpad=16)
cbar.set_ticks(np.arange(0, 1.01, 0.2))

plt.xlabel("UMAP_1", fontsize=14)
plt.ylabel("UMAP_2", fontsize=14)
plt.title("VAV1 - NGT-201-12 Co-folding Score", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("UMAP_ipTM.png",dpi=300)
plt.show()

# Color by chainpair_ipTM

umap.plot.output_file("UMAP_chainpair_ipTM.html")
chainpair_ipTM_plot = umap.plot.interactive(reducer, values=coordinates["pair_chains_ipTM_CRBN_VAV1_avg"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)  # cmap/value flag for continuous data

umap.plot.show(chainpair_ipTM_plot)

embedding = reducer.fit_transform(coordinates_data)

plt.figure(figsize=(7, 6))  
cmap = plt.cm.RdYlBu.reversed()
sc = plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=coordinates["pair_chains_ipTM_CRBN_VAV1_avg"].values,
    cmap=cmap,
    s=20,
    vmin=0,
    vmax=1.05
)

# Create the colorbar and customize it
cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
cbar.set_label('chainpair_ipTM', rotation=270, labelpad=16)
cbar.set_ticks(np.arange(0, 1.01, 0.2))

plt.xlabel("UMAP_1", fontsize=14)
plt.ylabel("UMAP_2", fontsize=14)
plt.title("VAV1 - NGT-201-12 Co-folding Score", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("UMAP_chainpair_ipTM.png",dpi=300)
plt.show()