# High-throughput laboratory evolution and evolutionary constraints in *Escherichia coli*

[![DOI](https://zenodo.org/badge/240467126.svg)](https://zenodo.org/badge/latestdoi/240467126)

This repository contains the codes used for the following paper:
> **High-throughput laboratory evolution and evolutionary constraints in *Escherichia coli***  
>Tomoya Maeda*, Junichiro Iwasawa*, Hazuki Kotani, Natsue Sakata, Masako Kawada, Takaaki Horinouchi, Aki Sakai, Kumi Tanabe, and Chikara Furusawa (*equal contribution)  
>[preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.19.956177v1)
>
> **Abstract:** *Understanding the constraints that shape the evolution of antibiotic resistance is critical for predicting and controlling drug resistance. Despite its importance, however, a systematic investigation for evolutionary constraints is lacking. Here, we performed a high-throughput laboratory evolution of Escherichia coli under the addition of 95 antibacterial chemicals and quantified the transcriptome, resistance, and genomic profiles for the evolved strains. Using interpretable machine learning techniques, we analyzed the phenotype-genotype data and identified low dimensional phenotypic states among the evolved strains. Further analysis revealed the underlying biological processes responsible for these distinct states, leading to the identification of novel trade-off relationships associated with drug resistance. We also report a novel constraint that leads to decelerated evolution. These findings bridge the genotypic, gene expression, and drug resistance space and lead to a better understanding of evolutionary constraints for antibiotic resistance.*

## Dependencies

We have tested this code using:

- Windows 10 / macOS 10.14.6
- Python 3.6.9

The full list of Python packages for the code is given in `requirements.txt`. These can be installed using:

```bash
pip install -r requirements.txt
```
This might take a few minutes.

## Usage

### supervised_pca_analysis.py

`supervised_pca_analysis.py` generates a heatmap for the stress resistance profiles of the 192 evolved strains sorted based on the hierarchical clustering in the supervised PCA space. The program can be executed by

```bash
python supervised_pca_analysis.py
```
Using the default options, the program would be executed within a few seconds.

- `compute_feature_importance` performs random forest regression to predict resistance profiles from gene expression levels, and returns a list of filtered genes that were informative for the prediction task. Please specify the PATH for the transcriptome data if you want to run the random forest model using `run_RF=True`. Re-running the random forest model would take about 6-7 minutes.

- `create_dendrogram` performs hierarchical clustering in the supervised PCA space. If `plot_figure=True`, this function generates a figure of the dendrogram (Fig. S3 (A)).

- `plot_resistance` plots the resistance profiles for all evolved strains in the order defined by hierarchical clustering. This figure corresponds to Fig. S3 (C).

## License
[MIT License](LICENSE)
