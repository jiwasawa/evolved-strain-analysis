# evolved-strain-analysis
Computational analysis for Maeda, Iwasawa et al., 2020.

#### supervised_pca_analysis.py
The program can be executed by
```
python supervised_pca_analysis.py
```
`supervised_pca_analysis.py` generates a heatmap for the stress resistance profiles of the 192 evolved strains sorted based on the hierarchical clustering in the supervised PCA space. Please specify the PATH for the transcriptome data before running the program.

- `compute_feature_importance` performs random forest regression to predict resistance profiles from gene expression levels, and returns a list of filtered genes that were informative for the prediction task.

- `create_dendrogram` performs hierarchical clustering in the supervised PCA space. If `plot_figure=True`, this function generates a figure of the dendrogram (Fig. S4 (A)).

- `plot_resistance` plots the resistance profiles for all evolved strains in the order defined by hierarchical clustering. This figure corresponds to Fig. S4 (C).
