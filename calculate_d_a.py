import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine

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


def find_top_bottom_k(df: pd.DataFrame, k: int) -> pd.DataFrame:
    results = {}
    for col_name in df.columns:
        col = df[col_name]
        top_k_indices = col.nlargest(k).index
        for i in range(k):
            if f"top_{i+1}_index" not in results:
                results[f"top_{i+1}_index"] = {}
            if f"top_{i+1}_int" not in results:
                results[f"top_{i+1}_int"] = {}
            results[f"top_{i+1}_index"][col_name] = top_k_indices[i]
            results[f"top_{i+1}_int"][col_name] = df.index.get_loc(top_k_indices[i])
        bottom_k_indices = col.nsmallest(k).index
        for i in range(k):
            if f"bottom_{i+1}_index" not in results:
                results[f"bottom_{i+1}_index"] = {}
            if f"bottom_{i+1}_int" not in results:
                results[f"bottom_{i+1}_int"] = {}
            results[f"bottom_{i+1}_index"][col_name] = bottom_k_indices[i]
            results[f"bottom_{i+1}_int"][col_name] = df.index.get_loc(
                bottom_k_indices[i]
            )
        # top_k_genes = df.loc[top_k_indices]
        # bottom_k_genes = df.loc[bottom_k_indices]
        # top_k_mean_norm = np.mean(np.linalg.norm(top_k_genes, axis=1))
        # bottom_k_mean_norm = np.mean(np.linalg.norm(bottom_k_genes, axis=1))
        top_k_values = col.loc[top_k_indices]
        bottom_k_values = col.loc[bottom_k_indices]
        top_k_mean_norm = np.mean(np.abs(top_k_values))
        bottom_k_mean_norm = np.mean(np.abs(bottom_k_values))
        if "top_abs" not in results:
            results["top_abs"] = {}
        results["top_abs"][col_name] = top_k_mean_norm
        if "bottom_abs" not in results:
            results["bottom_abs"] = {}
        results["bottom_abs"][col_name] = bottom_k_mean_norm
    results_df_transposed = pd.DataFrame(results)
    return results_df_transposed


def find_intruder(df: pd.DataFrame, top: float, upper: bool = True) -> pd.DataFrame:
    intruder_results = {}
    for current_col_name in df.columns:
        current_col = df[current_col_name]
        lower_half_threshold = current_col.median()
        intruder_found = False
        row_indices = df.index.to_list()
        np.random.shuffle(row_indices)
        for row_index in row_indices:
            row_data = df.loc[row_index]
            if (upper and row_data[current_col_name] <= lower_half_threshold) or (
                not upper and row_data[current_col_name] >= lower_half_threshold
            ):
                column_indices = df.columns.to_list()
                np.random.shuffle(column_indices)
                for other_col_name in column_indices:
                    if other_col_name != current_col_name:
                        other_col = df[other_col_name]
                        top_percent_threshold = other_col.quantile(
                            1 - top if upper else top
                        )
                        if (
                            upper and row_data[other_col_name] >= top_percent_threshold
                        ) or (
                            not upper
                            and row_data[other_col_name] <= top_percent_threshold
                        ):
                            intruder_results[current_col_name] = {
                                f"{"upper" if upper else "lower"}_intruder_index": row_index,
                                f"{"upper" if upper else "lower"}_intruder_int": df.index.get_loc(
                                    row_index
                                ),
                            }
                            intruder_found = True
                            break
                if intruder_found:
                    break
        if not intruder_found:
            raise ValueError(
                f"No {"upper" if upper else "lower"} intruder found for column {current_col_name} in the DataFrame."
            )
    result_df = pd.DataFrame(intruder_results).T
    return result_df


def d_a(
    df: pd.DataFrame,
    top_bottom_k: pd.DataFrame,
    upper_intruder: pd.DataFrame,
    lower_intruder: pd.DataFrame,
    k: int,
) -> pd.DataFrame:
    d_a_dict = {
        "top_d_a_intra": {},
        "top_d_a_inter": {},
        "bottom_d_a_intra": {},
        "bottom_d_a_inter": {},
    }
    for index in top_bottom_k.index:
        top_intra_sum = 0
        for i in range(k):
            for j in range(k):
                if i != j:
                    top_intra_sum += cosine(
                        df.loc[top_bottom_k.loc[index][f"top_{i+1}_index"]],
                        df.loc[top_bottom_k.loc[index][f"top_{j+1}_index"]],
                    )
        top_inter_sum = 0
        for i in range(k):
            top_inter_sum += cosine(
                df.loc[top_bottom_k.loc[index][f"top_{i+1}_index"]],
                df.loc[upper_intruder.loc[index]["upper_intruder_index"]],
            )
        bottom_intra_sum = 0
        for i in range(k):
            for j in range(k):
                if i != j:
                    bottom_intra_sum += cosine(
                        df.loc[top_bottom_k.loc[index][f"bottom_{i+1}_index"]],
                        df.loc[top_bottom_k.loc[index][f"bottom_{j+1}_index"]],
                    )
        bottom_inter_sum = 0
        for i in range(k):
            bottom_inter_sum += cosine(
                df.loc[top_bottom_k.loc[index][f"bottom_{i+1}_index"]],
                df.loc[lower_intruder.loc[index]["lower_intruder_index"]],
            )
        d_a_dict["top_d_a_intra"][index] = top_intra_sum / (k * (k - 1))
        d_a_dict["top_d_a_inter"][index] = top_inter_sum / k
        d_a_dict["bottom_d_a_intra"][index] = bottom_intra_sum / (k * (k - 1))
        d_a_dict["bottom_d_a_inter"][index] = bottom_inter_sum / k
    d_a_pd = pd.DataFrame(d_a_dict)
    return d_a_pd


def calculate_d_a(df: pd.DataFrame, k: int, top: float, times: int) -> pd.DataFrame:
    top_bottom_k = find_top_bottom_k(df, k)
    # upper_intruder = find_intruder(df, top, upper=True)
    # lower_intruder = find_intruder(df, top, upper=False)
    # d_a_pd = d_a(df, top_bottom_k, upper_intruder, lower_intruder, k)
    upper_intruders = []
    lower_intruders = []
    d_a_pds = []
    for i in range(times):
        upper_intruder = find_intruder(df, top, upper=True)
        lower_intruder = find_intruder(df, top, upper=False)
        d_a_pd = d_a(df, top_bottom_k, upper_intruder, lower_intruder, k)
        upper_intruder.columns = [f"group_{i}_{col}" for col in upper_intruder.columns]
        lower_intruder.columns = [f"group_{i}_{col}" for col in lower_intruder.columns]
        d_a_pd.columns = [f"group_{i}_{col}" for col in d_a_pd.columns]
        upper_intruders.append(upper_intruder)
        lower_intruders.append(lower_intruder)
        d_a_pds.append(d_a_pd)
    d_a_df = pd.concat(
        d_a_pds + upper_intruders + lower_intruders + [top_bottom_k], axis=1
    )
    return d_a_df


def main(
    dir: str = ".",
    prefix: str = "original",
    k: int = 5,
    top: float = 0.1,
    seed: int = 1027,
    times: int = 10,
):
    np.random.seed(seed)
    for model in models:
        for filtered in filtereds:
            file_prefix = prefix + "_" + model + ("_filtered" if filtered else "")
            csv_path = os.path.join(dir, file_prefix + ".csv")
            if os.path.isfile(csv_path):
                df = pd.read_csv(csv_path, index_col=0)
                d_a = calculate_d_a(df, k, top, times)
                d_a.to_csv(file_prefix + "_d_a.csv")
            for method in methods:
                rotated_file_prefix = file_prefix + f"_{method}_rotated_vectors"
                csv_path = os.path.join(dir, rotated_file_prefix + ".csv")
                if os.path.isfile(csv_path):
                    df = pd.read_csv(csv_path, index_col=0)
                    d_a = calculate_d_a(df, k, top, times)
                    d_a.to_csv(rotated_file_prefix + "_d_a.csv")


if "__main__" == __name__:
    import fire

    fire.Fire(main)
