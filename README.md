# Understanding-the-LLM-based-gene-embeddings

This is the code repository of the paper "Understanding the LLM-based gene embeddings" (Yufei Cai, Dailin Gan, Jun Li).

## Environment

The environment is set up using conda/24.7.1.

```{bash}
conda env create -f environment.yml -y
conda activate interpretation_genept
```

## OpenAI Embeddings

To get OpenAI embeddings (`text-embedding-3-small`, `text-embedding-3-large`, and `text-embedding-ada-002`) for the genes listed in `gene_info_table.csv`. We can use the following commands:

```{bash}
python get_embeddings.py
```

To get the half-length and symbol-only embeddings, we can use:

```{bash}
python get_embeddings_half_length.py
python get_embeddings_gene_name.py
```

After fetching these embeddings, we need to convert them into CSV files using:

```{bash}
python get_unrotated_vectors.py
```

To filter the gene embeddings according to `List of human protein-coding genes.csv`, we need to use:

```{bash}
python filter_dataset.py
```

## Completely randomized embedding matrix

To get the randomly-permuted embedding matrix, we need to use:

```{bash}
python get_permuted_unrotated_vectors.py
```

## Factor rotation

To do factor rotation on these embeddings, we can use the following command:

```{bash}
python get_rotated_vectors.py
```

The parameters of the `main` function can be modified:
```{python}
def main(
    prefix="original_text-embedding-3-small",
    permuted=False,
    permuted_index=1,
    method="Parsimax_orthogonal",
    tol=1e-5,
    filtered=True,
):
```
where the `method` parameter can be among `Quartimax_orthogonal`, `Varimax_orthogonal`, `Parsimax_orthogonal`, `FacParsim_orthogonal`, `Quartimax_oblique`, `Varimax_oblique`, `Parsimax_oblique`, and `FacParsim_oblique`.

The factor rotation depends on the repository [factor_rotation](https://github.com/mvds314/factor_rotation).

To summarize the overall distance ratio and related statistics, we can use:

```{bash}
python calculate_d_a.py
python dr_summary.py
python dr_overall.py
```

## FGSEA

The code for FGSEA is inside the folder `fgsea`. To prepare pathways for FGSEA analysis, we need to run the command first:

```{bash}
cd fgsea
Rscript prepare_pathways.R
```

The main code `fgsea.R` can be used by:

```{bash}
Rscript fgsea.R "h" "../original_text-embedding-3-small_filtered.csv" "original_h"
```
where the "h" stands for hallmark pathways, and can be replaced by "c2" representing C2 pathways. The "../original_text-embedding-3-small_filtered.csv" is the path to the CSV embedding file that needs analysis. The "original_h" is the output CSV filename, which means the results will be stored as "original_h.csv".

To summarize the pathway coverages for different embeddings, we need to use:

```{bash}
python summarize.py
```
and the results will be stored in `summary.csv`.