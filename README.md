# High-throughput laboratory evolution and evolutionary constraints in *Escherichia coli*

This repository contains the codes used for the following paper:
> **High-throughput laboratory evolution and evolutionary constraints in *Escherichia coli***  
>Tomoya Maeda*, Junichiro Iwasawa*, Hazuki Kotani, Natsue Sakata, Masako Kawada, Takaaki Horinouchi, Aki Sakai, Kumi Tanabe, and Chikara Furusawa (*equal contribution)  
>[preprint on bioRxiv](https://www.biorxiv.org/)
>
> **Abstract:** *Understanding the constraints that shape the evolution of antibiotic resistance is critical for predicting and controlling drug resistance. Despite its importance, however, a systematic investigation for evolutionary constraints is lacking. Here, we performed a high-throughput laboratory evolution of Escherichia coli under the addition of 95 antibacterial chemicals and quantified the transcriptome, resistance, and genomic profiles for the evolved strains. Using interpretable machine learning techniques, we analyzed the phenotype-genotype data and identified low dimensional phenotypic states among the evolved strains. Further analysis revealed the underlying biological processes responsible for these distinct states, leading to the identification of novel trade-off relationships associated with drug resistance. We also report a novel constraint that leads to decelerated evolution. These findings bridge the genotypic, gene expression, and drug resistance space and lead to a better understanding of evolutionary constraints for antibiotic resistance.*

## supervised_pca_analysis.py

`supervised_pca_analysis.py` generates a heatmap for the stress resistance profiles of the 192 evolved strains sorted based on the hierarchical clustering in the supervised PCA space. The program can be executed by

```console
python supervised_pca_analysis.py
```

- `compute_feature_importance` performs random forest regression to predict resistance profiles from gene expression levels, and returns a list of filtered genes that were informative for the prediction task. Please specify the PATH for the transcriptome data if you want to run the random forest model using `run_RF=True`.

- `create_dendrogram` performs hierarchical clustering in the supervised PCA space. If `plot_figure=True`, this function generates a figure of the dendrogram (Fig. S3 (A)).

- `plot_resistance` plots the resistance profiles for all evolved strains in the order defined by hierarchical clustering. This figure corresponds to Fig. S3 (C).
