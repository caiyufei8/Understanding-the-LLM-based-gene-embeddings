import pickle
import numpy as np
import pandas as pd

models = (
    "text-embedding-3-small",
    # "text-embedding-3-large",
    # "text-embedding-ada-002",
)
files_to_names = {
    "gene_name_to_embedding.pkl": "original",
    # "gene_name_to_embedding_gene_name.pkl": "gene_names",
    # "gene_name_to_embedding_half_length.pkl": "half_length",
}


def permute_df(df: pd.DataFrame) -> pd.DataFrame:
    flat_array = df.to_numpy().flatten()
    np.random.shuffle(flat_array)
    permuted_df = pd.DataFrame(
        flat_array.reshape(df.shape), index=df.index, columns=df.columns
    )
    return permuted_df


def main(seed: int = 1027, times: int = 10):
    for name in files_to_names.values():
        for model in models:
            np.random.seed(seed)
            df = pd.read_csv(f"{name}_{model}_filtered.csv", index_col=0)
            permuted_df = permute_df(df)
            permuted_df.to_csv(f"{name}_{model}_filtered_permuted.csv", index=True)
            for i in range(2, times + 1):
                np.random.seed(seed + i)
                permuted_df = permute_df(permuted_df)
                permuted_df.to_csv(
                    f"{name}_{model}_filtered_permuted{i}.csv", index=True
                )


if "__main__" == __name__:
    main()
