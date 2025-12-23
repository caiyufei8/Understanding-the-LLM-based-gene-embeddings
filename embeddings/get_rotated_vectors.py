import numpy as np
import pandas as pd
import factor_rotation as fr
import torch, os

models = (
    "text-embedding-3-small",
    "text-embedding-3-large",
    "text-embedding-ada-002",
)
methods = {
    "Quartimax_orthogonal": "quartimax",
    "Varimax_orthogonal": "varimax",
    "Parsimax_orthogonal": "parsimax",
    "FacParsim_orthogonal": "parsimony",
    "Quartimax_oblique": "quartimin",
    "Varimax_oblique": "covarimin",
    "Parsimax_oblique": "parsimax_oblique",
    "FacParsim_oblique": "parsimony_oblique",
}
files_to_names = {
    "gene_name_to_embedding.pkl": "original",
    "gene_name_to_embedding_gene_name.pkl": "gene_names",
    "gene_name_to_embedding_half_length.pkl": "half_length",
}


def main(
    prefix="original_text-embedding-3-small",
    permuted=False,
    permuted_index=1,
    method="Parsimax_orthogonal",
    tol=1e-5,
    filtered=True,
):
    df = pd.read_csv(
        prefix
        + ("_filtered" if filtered else "")
        + ("_permuted"+(str(permuted_index) if permuted_index!=1 else "") if permuted else "")
        + ".csv",
        index_col=0,
    )
    print(df.head())
    if method in methods:
        if not os.path.isfile(
            prefix
            + ("_filtered" if filtered else "")
            + ("_permuted"+(str(permuted_index) if permuted_index!=1 else "") if permuted else "")
            + f"_{method}_rotated_vectors.csv"
        ):
            df_train = df.to_numpy()
            df_train /= np.max(np.abs(df_train))
            L, T = fr.rotate_factors(
                df_train,
                methods[method],
                dtype=torch.float64,
                device=torch.device("cuda"),
                tol=tol,
            )
            rotated_df = pd.DataFrame(
                L.cpu().numpy(),
                index=df.index,
                columns=df.columns,
            )
            rotated_df.to_csv(
                prefix
                + ("_filtered" if filtered else "")
                + ("_permuted"+(str(permuted_index) if permuted_index!=1 else "") if permuted else "")
                + f"_{method}_rotated_vectors.csv"
            )
            axis_df = pd.DataFrame(
                T.cpu().numpy(), index=df.columns, columns=df.columns
            )
            axis_df.to_csv(
                prefix
                + ("_filtered" if filtered else "")
                + ("_permuted"+(str(permuted_index) if permuted_index!=1 else "") if permuted else "")
                + f"_{method}_axis.csv"
            )


if "__main__" == __name__:
    import fire

    fire.Fire(main)
