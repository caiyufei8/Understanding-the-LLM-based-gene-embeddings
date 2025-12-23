import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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


def calculate_dr_summary(df: pd.DataFrame, times: int) -> float:
    top_dr_a = 0
    bottom_dr_a = 0
    top_dr_a_inter = 0
    top_dr_a_intra = 0
    bottom_dr_a_inter = 0
    bottom_dr_a_intra = 0
    for i in range(times):
        top_dr_a += df[f"group_{i}_top_d_a_inter"] / df[f"group_{i}_top_d_a_intra"]
        bottom_dr_a += (
            df[f"group_{i}_bottom_d_a_inter"] / df[f"group_{i}_bottom_d_a_intra"]
        )
        top_dr_a_inter += df[f"group_{i}_top_d_a_inter"]
        top_dr_a_intra += df[f"group_{i}_top_d_a_intra"]
        bottom_dr_a_inter += df[f"group_{i}_bottom_d_a_inter"]
        bottom_dr_a_intra += df[f"group_{i}_bottom_d_a_intra"]
    top_dr_a = top_dr_a.to_numpy().ravel()
    bottom_dr_a = bottom_dr_a.to_numpy().ravel()
    top_dr_a_inter = top_dr_a_inter.to_numpy().ravel()
    top_dr_a_intra = top_dr_a_intra.to_numpy().ravel()
    bottom_dr_a_inter = bottom_dr_a_inter.to_numpy().ravel()
    bottom_dr_a_intra = bottom_dr_a_intra.to_numpy().ravel()
    top_dr_a /= times
    bottom_dr_a /= times
    top_dr_a_inter /= times
    top_dr_a_intra /= times
    bottom_dr_a_inter /= times
    bottom_dr_a_intra /= times
    top_abs = df["top_abs"].to_numpy().ravel()
    bottom_abs = df["bottom_abs"].to_numpy().ravel()
    # print(
    #     top_dr_a,
    #     bottom_dr_a,
    #     top_dr_a_inter,
    #     top_dr_a_intra,
    #     bottom_dr_a_inter,
    #     bottom_dr_a_intra,
    #     top_abs,
    #     bottom_abs,
    # )
    # exit()
    summary = {}
    summary["A"] = np.mean(top_dr_a)
    summary["B"] = np.mean(bottom_dr_a)
    summary["C"], summary["C_p_value"] = spearmanr(top_dr_a, bottom_dr_a)
    summary["D"], summary["D_p_value"] = spearmanr(
        np.concatenate([top_dr_a, bottom_dr_a]), np.concatenate([top_abs, bottom_abs])
    )
    summary["E"], summary["E_p_value"] = spearmanr(
        np.concatenate([top_dr_a_intra, bottom_dr_a_intra]),
        np.concatenate([top_abs, bottom_abs]),
    )
    summary["F"], summary["F_p_value"] = spearmanr(
        np.concatenate([top_dr_a_inter, bottom_dr_a_inter]),
        np.concatenate([top_abs, bottom_abs]),
    )
    selected_dr_a = np.array(
        [
            top_dr_a[i] if top_abs[i] >= bottom_abs[i] else bottom_dr_a[i]
            for i in range(len(top_dr_a))
        ]
    )
    summary["G"] = np.mean(selected_dr_a)
    return summary


def main(dir: str = ".", prefix="original", times: int = 10):
    summary = {}
    for model in models:
        for filtered in filtereds:
            file_prefix = prefix + "_" + model + ("_filtered" if filtered else "")
            csv_path = os.path.join(dir, file_prefix + "_d_a.csv")
            if os.path.isfile(csv_path):
                df = pd.read_csv(csv_path, index_col=0)
                summary[file_prefix] = calculate_dr_summary(df, times)
            for method in methods:
                rotated_file_prefix = file_prefix + f"_{method}_rotated_vectors"
                csv_path = os.path.join(dir, rotated_file_prefix + "_d_a.csv")
                if os.path.isfile(csv_path):
                    df = pd.read_csv(csv_path, index_col=0)
                    summary[rotated_file_prefix] = calculate_dr_summary(df, times)
    summary_df = pd.DataFrame(summary).T
    summary_df.to_csv("dr_summary.csv")


if "__main__" == __name__:
    import fire

    fire.Fire(main)
