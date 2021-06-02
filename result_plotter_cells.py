import scanpy.api as sc
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import json
import glob
import re
from matplotlib import cm
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import os
import seaborn as sns
import pandas as pd
import scipy
plt.rcParams["figure.figsize"] = (15,12)
import json



import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input', required=True,
                    help='Path to the generated sc files folder')
parser.add_argument('--ref', required=True,
                    help='Path to full original dataset')
parser.add_argument('--dataset', required=True,
                    help='Name of the dataset ')
parser.add_argument('--plot_ref', required=False,
                    help='Whether the full dataset should be also plotted')
a = parser.parse_args()

refdata = sc.read(a.ref)
input_folder=a.input
dataset = a.dataset
plot_all = False
if a.plot_ref:
    plot_all = True

def array_euc_distance(A, B):
    M = A.shape[0]
    N = B.shape[0]
    
    A_sq = (A*A).sum(axis=1).reshape((M,1))*np.ones(shape=(1,N))
    B_sq = (B*B).sum(axis=1)*np.ones(shape=(M,1))
    AB2 = -2*A.dot(B.T)
    
    dist_array = A_sq + B_sq + AB2
    
    return dist_array

sc.pp.neighbors(refdata) # UMAP is based on the neighbor graph; we'll compute this first
sc.tl.umap(refdata, n_components=2)
os.chdir(input_folder)
full = sc.pl.umap(refdata, 
                  color='cluster', 
                  alpha=0.5, 
                  size= 25, 
                  legend_loc='on data', 
                  save="All_real_cells.pdf", 
                  title="UMAP plot of full reference dataset", 
                  show=False)
plt.close('all')
#full.savefig("".join((unput folder, "All_real_cells.pdf")))
os.chdir("../..")

regex=re.compile(r'\d+')

plt.rcParams["figure.figsize"] = (30,20)

generated_files = glob.glob("/".join((input_folder, "*.h5ad")))
generated_files = np.sort(generated_files)

save_path="".join((input_folder, 'figures/UMAP__plots.pdf'))
RMSE = {}

with PdfPages(save_path) as pages:
    
    for file in generated_files:
        adata = sc.read(file)
        #print(adata.var_names)
        number=np.array(regex.findall(file))[-2].astype(int)
        print("Plotting sample %d" % (number))
        with open("".join(["/ngc/people/viklav/bulk2sc/Bulk2SC/Datasets/",
                           dataset,
                           "/Deconvolution/reference/cells.json"]), 'r') as fp:
            refcells = np.array(json.load(fp))[number-1, :].tolist()
        refdata_s = refdata[refcells]
        adata.var.index = refdata.var.index
        joint_data = adata.concatenate(refdata, join='inner')
        
        
        real_props = pd.read_table("".join(["/ngc/people/viklav/bulk2sc/Bulk2SC/Datasets/",
                           dataset,
                           "/Deconvolution/reference/props.txt"]), sep=" ").iloc[:, number-1]
                                   
        est_props = pd.read_table("".join([input_folder, "cell_proportions.txt"]),
                                  index_col=0).transpose().iloc[:, number-1]
        MSE = np.square(np.subtract(real_props,est_props)).mean() 
 
        RMSE_t = np.sqrt(MSE)
        RMSE[str(number)] = RMSE_t
        
        sc.pp.neighbors(joint_data) # UMAP is based on the neighbor graph; we'll compute this first
        sc.tl.umap(joint_data, n_components=2)
        
        new_real = []
        for item in refcells:
            new_real.append("".join([item, "-1"]))
        real_cells = joint_data[new_real]
        
        clusters = np.array(joint_data.obs["cluster"]).astype(int)
        real_clusters = np.array(real_cells.obs["cluster"]).astype(int)
        batches = np.array(joint_data.obs["batch"]).astype(int)
        
        embedded_cells_train = joint_data.obsm['X_umap'][adata.shape[0]:, :]
        embedded_cells_real = real_cells.obsm['X_umap']
        embedded_cells_fake = joint_data.obsm['X_umap'][:adata.shape[0], :]
        
        clusters_no=np.unique(clusters)
        colormap = cm.nipy_spectral
        colors = [colormap(i) for i in np.linspace(0, 1, np.max(clusters_no)+3)]

        fig, ax = plt.subplots(2, 2, sharey=False, sharex=False)
        #ax1 = fig.gca()
        if plot_all == True:
            for i in clusters_no:
                mask = clusters[adata.shape[0]:] == i
                if sum(mask)==0:
                    continue
                ax[0, 0].scatter(embedded_cells_train[mask, 0],
                            embedded_cells_train[mask, 1],
                            color=colors[i], marker='x', edgecolors="white",
                            alpha=1000/embedded_cells_train.shape[0],
                            label='training_' + str(i))
                ax[0,1].scatter(embedded_cells_train[mask, 0],
                            embedded_cells_train[mask, 1],
                            color=colors[i], marker='x', edgecolors="white",
                            alpha=1000/embedded_cells_train.shape[0],
                            label='training_' + str(i))
        for i in clusters_no:
            mask = clusters[:adata.shape[0]] == i
            if sum(mask)==0:
                continue
            ax[0, 0].scatter(embedded_cells_fake[mask, 0],
                        embedded_cells_fake[mask, 1],
                        color=colors[i], marker='+', edgecolors="white",
                        alpha=1,
                        label='fake_' + str(i))
            ax[1,0].scatter(embedded_cells_fake[mask, 0],
                        embedded_cells_fake[mask, 1],
                        color=colors[i], marker='+', edgecolors="white",
                        alpha=1,
                        label='fake_' + str(i))


        for i in clusters_no:
            mask = real_clusters == i
            if sum(mask)==0:
                continue

            ax[0,1].scatter(embedded_cells_real[mask, 0],
                        embedded_cells_real[mask, 1],
                        color=colors[i], marker='o',
                        alpha=0.4,
                        #edgecolor='black',
                        #linewidths=0.05,
                        label='real_' + str(i))
            ax[1,0].scatter(embedded_cells_real[mask, 0],
                        embedded_cells_real[mask, 1],
                        color=colors[i], marker='o',
                        alpha=0.4,
                        #edgecolor='black',
                        #linewidths=0.05,
                        label='real_' + str(i))

        ax[1,1].set(xlim=(-1, np.max(real_props.index.tolist())+1),
                    ylim=(0, np.max([np.max(real_props), np.max(est_props)])+0.1))
        ax[1,1].bar(real_props.index-0.125, real_props, label= "True", color="green", width= 0.25)
        ax[1,1].bar(real_props.index+0.125, est_props, width= 0.25, label="Estimated", color="red")
        ax[1,1].legend()
        ax[1,1].set_xticks(real_props.index)
        ax[1,1].set_title("True and estimated cell-type proportions")
                                   
                                   
        ax[0,0].grid(True)
        ax[0,1].grid(True)
        ax[1,0].grid(True)
        ax[1,1].axis(True)
        ax[0,0].set_title(" ".join(('UMAP plot of fake cells of sample', str(number))))
        ax[0,1].set_title(" ".join(('UMAP plot of real cells of sample', str(number))))
        ax[1,0].set_title(" ".join(('Merged (sample', str(number), ')')))
        
        ax[1,0].legend(loc='lower left',
                    numpoints=1, ncol=3,
                    fontsize=8, bbox_to_anchor=(0, 0))
        ax[0,0].legend(loc='lower left',
                    numpoints=1, ncol=3,
                    fontsize=8, bbox_to_anchor=(0, 0))
        ax[0,1].legend(loc='lower left',
                    numpoints=1, ncol=3,
                    fontsize=8, bbox_to_anchor=(0, 0))
        
        

        
        canvas = FigureCanvasPdf(fig)
        canvas.print_figure(pages)
        
        
with open("/".join((input_folder,'RMSE.json')), 'w') as fp:
    json.dump(RMSE, fp)
        
save_path_heat="".join((input_folder, 'figures/Heatmaps/'))

try:
    os.mkdir(save_path_heat)
except:
    pass


for file in generated_files:
        adata = sc.read(file)
        number=np.array(regex.findall(file))[-2].astype(int)
        print("Plotting sample %d" % (number))
        with open("".join(["/ngc/people/viklav/bulk2sc/Bulk2SC/Datasets/",
                           dataset, 
                           "/Deconvolution/reference/cells.json"]), 'r') as fp:
            refcells = np.array(json.load(fp))[number-1, :].tolist()
        refdata_s = refdata[refcells]
        adata.var.index = refdata_s.var.index
        
        ref_ind = np.array(np.argsort(refdata_s.obs["cluster"]))
        a_ind= np.array(np.argsort(adata.obs["cluster"]))
        refdata_s = refdata_s[ref_ind]
        adata_s = adata[a_ind]
     
        adf = np.array(scipy.sparse.csr_matrix.todense(adata_s.X))
        refdf = np.array(refdata_s.X)
        
        distance = pd.DataFrame(np.log(array_euc_distance(refdf, adf)))
        
        
        clusters_present= np.unique(np.concatenate([refdata_s.obs["cluster"], adata_s.obs["cluster"]])).astype(int)
        #print(clusters_present)
        cluster_max=len(clusters_present)
        
        
        cmap = sns.color_palette("viridis", n_colors=2000)
        labels1= np.array(refdata_s.obs["cluster"]).astype(int).tolist()
        labels2= np.array(adata.obs["cluster"]).astype(int).tolist()

        lut = dict(zip(set(clusters_present), sns.color_palette("Spectral", cluster_max)))
        col_colors = pd.DataFrame(labels2, columns=["fake_cluster"])["fake_cluster"].map(lut)
        row_colors = pd.DataFrame(labels1, columns=["real_cluster"])["real_cluster"].map(lut)



        g =sns.clustermap(distance.reset_index(drop=True),
                          cmap=cmap,
                          row_colors = row_colors,
                          col_colors=col_colors,
                         row_cluster=False, col_cluster=False)
        for label in clusters_present:
            g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                            label=label, linewidth=0)
            g.ax_col_dendrogram.legend(loc="center", ncol=cluster_max)
            
        g.ax_row_dendrogram.set_visible(False)
        g.ax_row_dendrogram.set_xlim([0,0])
        g.cax.set_position([.97, .2, .03, .45])
        g.savefig("".join((save_path_heat, "sample_%d" % (number))))

    
    

