# evolved-strain-analysis
Computational analysis for Maeda, Iwasawa et al., 2020.

#### supervised_pca_analysis.py
- `compute_feature_importance` performs random forest regression to predict resistance profiles from gene expression levels, and returns a list of filtered genes that were informative for the prediction task.

- `create_dendrogram` performs hierarchical clustering in the supervised PCA space. If `plot_figure=True`, this function generates a figure of the dendrogram (Fig. S4 (A)).

- `plot_resistance` plots the resistance profiles for all evolved strains in the order defined by hierarchical clustering. This figure corresponds to Fig. S4 (C).