import os
import pandas as pd

models = (
    "text-embedding-3-small",
    # "text-embedding-3-large",
    # "text-embedding-ada-002",
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
tols = {1e-5}
# filtereds = {True, False}
filtereds = {True}


def calculate_dr_overall(df: pd.DataFrame, times: int) -> float:
    sum_dr_a = 0
    for i in range(times):
        sum_dr_a += (
            df[f"group_{i}_top_d_a_inter"] / df[f"group_{i}_top_d_a_intra"]
        ).sum() + (
            df[f"group_{i}_bottom_d_a_inter"] / df[f"group_{i}_bottom_d_a_intra"]
        ).sum()
        # sum_dr_a += (
        #     df[f"group_{i}_top_d_a_inter"] / df[f"group_{i}_top_d_a_intra"]
        # ).sum()
    return sum_dr_a / (2 * df.shape[0] * times)
    # return sum_dr_a / (df.shape[0] * times)


def main(dir: str = ".", prefix="original", times: int = 10):
    dr = {}
    for model in models:
        for filtered in filtereds:
            file_prefix = prefix + "_" + model + ("_filtered" if filtered else "")
            csv_path = os.path.join(dir, file_prefix + "_d_a.csv")
            if os.path.isfile(csv_path):
                df = pd.read_csv(csv_path, index_col=0)
                dr[file_prefix] = calculate_dr_overall(df, times)
            for method in methods:
                rotated_file_prefix = file_prefix + f"_{method}_rotated_vectors"
                csv_path = os.path.join(dir, rotated_file_prefix + "_d_a.csv")
                if os.path.isfile(csv_path):
                    df = pd.read_csv(csv_path, index_col=0)
                    dr[rotated_file_prefix] = calculate_dr_overall(df, times)
    dr_df = pd.DataFrame(dr, index=[0]).T
    dr_df.to_csv("dr_overall.csv")


if "__main__" == __name__:
    import fire

    fire.Fire(main)
