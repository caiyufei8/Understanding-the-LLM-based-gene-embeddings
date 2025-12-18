import pickle
import pandas as pd

models = (
    "text-embedding-3-small",
    "text-embedding-3-large",
    "text-embedding-ada-002",
)
files_to_names = {
    "gene_name_to_embedding.pkl": "original",
    "gene_name_to_embedding_gene_name.pkl": "gene_names",
    "gene_name_to_embedding_half_length.pkl": "half_length",
}


def main():
    for file in files_to_names:
        gene_name_to_embedding = pickle.load(open(file, "rb"))
        for model in models:
            df = pd.DataFrame.from_dict(gene_name_to_embedding[model], orient="index")
            df.to_csv(f"{files_to_names[file]}_{model}.csv", index=True)


if "__main__" == __name__:
    main()
