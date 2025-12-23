import os
import pandas as pd

models = (
    "text-embedding-3-small",
    "text-embedding-3-large",
    "text-embedding-ada-002",
    "bge-small-en-v1.5",
    "biobert-base-cased-v1.1",
    "e5-small-v2",
    "e5-small",
    "GIST-all-MiniLM-L6-v2",
    "GIST-small-Embedding-v0",
    "gte-small",
    "gte-tiny",
    "MedEmbed-small-v0.1",
    "NoInstruct-small-Embedding-v0",
    "stella-base-en-v2",
)
files_to_names = {
    "gene_name_to_embedding.pkl": "original",
    "gene_name_to_embedding_gene_name.pkl": "gene_names",
    "gene_name_to_embedding_half_length.pkl": "half_length",
}


def main():
    df_csv = pd.read_csv(
        "List of human protein-coding genes.csv", usecols=[1], skiprows=[0]
    )
    index_list_from_csv = set(df_csv.iloc[:, 0].tolist())
    while "" in index_list_from_csv:
        index_list_from_csv.remove("")
    for name in files_to_names.values():
        for model in models:
            if os.path.isfile(f"{name}_{model}.csv"):
                df = pd.read_csv(f"{name}_{model}.csv", index_col=0)
                print(df.head())
                df_filtered = df[df.index.isin(index_list_from_csv)]
                df_filtered.to_csv(f"{name}_{model}_filtered.csv")


if "__main__" == __name__:
    import fire

    fire.Fire(main)
