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

def compute_feature_importance(expression, resistance_data):
    """Predict resistance profiles from expression profiles and 
    raise important genes for the Random Forest regression."""
    X = expression.T.iloc[:-4, :] # exclude the parent strains from the feature matrix X
    y = resistance_data.reindex(X.index)

    # Random Forest regression using hyperparameters defined by grid search
    RF_reg = RandomForestRegressor(n_estimators=300, max_depth=18, random_state=42) 
    RF_reg.fit(X, y)

    df = expression.iloc[np.argsort(RF_reg.feature_importances_)[::-1],:].T # sort genes based on its importance 
    df.to_pickle('./supervised_expression.pkl')
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
    return hierarchy2

expression = pd.read_excel('/PATH_TO_TABLE_S3/Table S3. Transcriptome data of evolved strains.xls', \
    index_col=0, skiprows=1)
resistance_data = pd.read_csv('./resistance_norm.csv', index_col=0) # normalized IC50 values


#supervised_expression = compute_feature_importance(expression, resistance_data)
supervised_expression = pd.read_pickle('./supervised_expression.pkl')
supervised_expression = supervised_expression.iloc[:,:213] # threshold for the genes

pca = PCA(n_components=36, svd_solver='full')
df_pca = pca.fit_transform(supervised_expression)
df_pca = pd.DataFrame(df_pca,index=supervised_expression.index)

hiearchy2 = create_dendrogram(df_pca)