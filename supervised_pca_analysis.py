import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['lines.linewidth'] = 3
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy.spatial.distance import pdist

def compute_feature_importance(expression, resistance_data, run_RF=False):
    """Predict resistance profiles from expression profiles and 
    raise important genes for the Random Forest regression.
    When run_RF=False, the already filtered expression data is used."""
    if run_RF:
        X = expression.T.iloc[:-4, :] # exclude the parent strains from the feature matrix X
        y = resistance_data.reindex(X.index)

        # Random Forest regression using hyperparameters defined by grid search
        RF_reg = RandomForestRegressor(n_estimators=300, max_depth=18, random_state=42) 
        RF_reg.fit(X, y)

        df = expression.iloc[np.argsort(RF_reg.feature_importances_)[::-1],:].T # sort genes based on its importance 
        df.to_pickle('./filtered_expression.pkl')
    else:
        df = pd.read_pickle('./data/filtered_expression.pkl')
    return df 

def create_dendrogram(df_pca, plot_figure=False):
    """Hierarchical clustering in the supervised PCA space."""

    linked = linkage(pdist(df_pca[:-4],metric='euclidean'), 'ward', optimal_ordering=False)
    leave_array = hierarchy.leaves_list(linked)
    leave_array = leave_array[leave_array<192]
    strain_names = df_pca.index.tolist()[:-4]
    strain_h = [strain_names[leave_array[i]] for i in range(len(leave_array))]
    hierarchy2 = fcluster(linked, t=15, criterion='maxclust') # cluster number is set to 15
    
    if plot_figure:
        color_list = ['#035ea2']*11
        hierarchy.set_link_color_palette(color_list)
        plt.figure(figsize=(40,3))
        ax = plt.axes()
        dn = dendrogram(linked,color_threshold=4.4835,labels=df_pca.index[:-4],#distance_sort='descending',
                        above_threshold_color='grey',leaf_rotation=90,get_leaves=True,distance_sort=True,
                        leaf_font_size=8
                    )
        plt.yticks([])
        for j,strain_ in enumerate(strain_h):
                    plt.text((j+0.15)/len(strain_h),-0.008,s=strain_,fontsize=13,
                            color ='black',
                            transform=ax.transAxes,verticalalignment='top',rotation=90,
                            weight='normal')
        plt.axis('off')
        #plt.savefig('FigS4_A.pdf', dpi=400, bbox_inches='tight')
        plt.show()
    return strain_h, hierarchy2

def plot_resistance(resistance_data, strain_h):
    # hierarchical clustering of stresses
    custom_cmap = sns.diverging_palette(252,18,s=99,l=52,sep=10,center='light',as_cmap=True)

    # sorting of strains based on the hierarchical clustering in the supervised PCA space
    stress_order = list(np.load('./data/stress_list_MAGE_order.npy')) # order of stresses based on mutant strain analysis
    cl_resistance = resistance_data.reindex(strain_h).T
    cl_resistance = cl_resistance.reindex(stress_order)

    fig = plt.figure(figsize=(40,14))
    ax = plt.axes()
    mic_map = ax.imshow(cl_resistance, aspect=1.36, cmap=custom_cmap, clim=(-3,3))
    cb = fig.colorbar(mic_map, extend='both', orientation='horizontal',ticks=[-3,3],aspect=3,
                    shrink=0.04,use_gridspec=False,anchor=(0.855,1.45))
    cb.ax.set_xticklabels(cb.ax.get_xticklabels(), fontsize=25)

    for j,stress_ in enumerate(stress_order):
        plt.text(-0.005, (len(stress_order)-j-0.8)/len(stress_order), s=stress_,fontsize=13,
                color ='black', transform=ax.transAxes, horizontalalignment='right', weight='normal')
    plt.xlim(-0.5,191.6)
    plt.axis('off')
    for j,strain_ in enumerate(strain_h):
            plt.text((j+0.15)/len(strain_h), -0.005, s=strain_, fontsize=10,
                    color ='black', transform=ax.transAxes,verticalalignment='top',rotation=90,
                    weight='normal')
    #plt.savefig('FigS4_C.pdf', dpi=400, bbox_inches='tight')
    plt.show()


expression = pd.read_excel('/PATH_TO_TABLE_S3/Table S3. Transcriptome data of evolved strains.xlsx', index_col=0, skiprows=1)
resistance_data = pd.read_csv('./data/resistance_norm.csv', index_col=0) # normalized IC50 values

filtered_expression = compute_feature_importance(expression, resistance_data)
filtered_expression = filtered_expression.iloc[:,:213] # threshold for the genes

# perform supervised PCA
pca = PCA(n_components=36, svd_solver='full')
df_pca = pca.fit_transform(filtered_expression)
df_pca = pd.DataFrame(df_pca,index=filtered_expression.index)

strain_h, hiearchy2 = create_dendrogram(df_pca)

plot_resistance(resistance_data, strain_h)