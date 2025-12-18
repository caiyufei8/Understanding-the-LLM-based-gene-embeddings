import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from statsmodels.stats.multitest import multipletests

model_dim = {
    "original": 1536,
    "original_permuted": 1536,
    "original_permuted2": 1536,
    "original_permuted3": 1536,
    "original_permuted4": 1536,
    "original_permuted5": 1536,
    "original_permuted6": 1536,
    "original_permuted7": 1536,
    "original_permuted8": 1536,
    "original_permuted9": 1536,
    "original_permuted10": 1536,
    "gene_names": 1536,
    "half_length": 1536,
    "bge-small-en-v1.5": 384,
    "gene_names_bge-small-en-v1.5": 384,
    "half_length_bge-small-en-v1.5": 384,
    "biobert-base-cased-v1.1": 768,
    "gene_names_biobert-base-cased-v1.1": 768,
    "half_length_biobert-base-cased-v1.1": 768,
    "e5-small-v2": 384,
    "gene_names_e5-small-v2": 384,
    "half_length_e5-small-v2": 384,
    "e5-small": 384,
    "gene_names_e5-small": 384,
    "half_length_e5-small": 384,
    "GIST-all-MiniLM-L6-v2": 384,
    "gene_names_GIST-all-MiniLM-L6-v2": 384,
    "half_length_GIST-all-MiniLM-L6-v2": 384,
    "GIST-small-Embedding-v0": 384,
    "gene_names_GIST-small-Embedding-v0": 384,
    "half_length_GIST-small-Embedding-v0": 384,
    "gte-small": 384,
    "gene_names_gte-small": 384,
    "half_length_gte-small": 384,
    "gte-tiny": 384,
    "gene_names_gte-tiny": 384,
    "half_length_gte-tiny": 384,
    "MedEmbed-small-v0.1": 384,
    "gene_names_MedEmbed-small-v0.1": 384,
    "half_length_MedEmbed-small-v0.1": 384,
    "NoInstruct-small-Embedding-v0": 384,
    "gene_names_NoInstruct-small-Embedding-v0": 384,
    "half_length_NoInstruct-small-Embedding-v0": 384,
    "stella-base-en-v2": 768,
    "gene_names_stella-base-en-v2": 768,
    "half_length_stella-base-en-v2": 768,
}
methods = {
    # "Quartimax_orthogonal": "quartimax",
    # "Varimax_orthogonal": "varimax",
    "Parsimax_orthogonal": "parsimax",
    # "FacParsim_orthogonal": "parsimony",
    # "Quartimax_oblique": "quartimin",
    # "Varimax_oblique": "covarimin",
    # "Parsimax_oblique": "parsimax_oblique",
    # "FacParsim_oblique": "parsimony_oblique",
}
pathway_nums = {"h": 50, "c2": 3077}
permuteds = [False, True]
# tols = {1e-5}
# filtereds = {True, False}
# num_genes_list = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]


def find_label(key):
    for label in methods:
        if label in key:
            return label
    return "Original"


def count_elements(total: int, input_list: list) -> dict:
    element_counts = Counter(input_list)
    count_dict = {}
    count_dict[0] = total - len(set(input_list))
    for element, count in element_counts.items():
        if count not in count_dict:
            count_dict[count] = 0
        count_dict[count] += 1
    return count_dict


def calculate_average_counts(count_dict: dict) -> float:
    if not count_dict:
        return 0.0
    total_counts = sum(count * freq for count, freq in count_dict.items())
    total_elements = sum(freq for freq in count_dict.values())
    if total_elements == 0:
        return 0.0
    return total_counts / total_elements


def align_counts(count_dicts: dict) -> dict:
    total_max_i = 0
    for key in count_dicts:
        max_i = 0
        for i in count_dicts[key]:
            max_i = max(max_i, i)
        count_dicts[key] = {
            i: (count_dicts[key][i] if i in count_dicts[key] else 0)
            for i in range(max_i + 1)
        }
        total_max_i = max(total_max_i, max_i)
    for key in count_dicts:
        max_i = 0
        for i in count_dicts[key]:
            max_i = max(max_i, i)
        for i in range(max_i + 1, total_max_i + 1):
            count_dicts[key][i] = 0
    return count_dicts


def fdr_filter(df, p, prefix):
    # 1
    df = df[df["padj"] < p / model_dim[prefix]]

    # 2
    # df=df[df['pval'].notna()]
    # fdr_pvals = multipletests(df["pval"].to_numpy().ravel(), method='fdr_bh')[1]
    # df=df[fdr_pvals<p]

    # 3
    # df = df[df["padj"] < p]

    # 4
    # df = df[df["pval"] < p / model_dim[prefix]]
    return df


def main(dir: str = ".", p=0.05):
    output = {}
    for pathway_type in pathway_nums:
        pathway_counts = {}
        average_pathway_counts = {}
        dimension_counts = {}
        average_dimension_counts = {}
        pathway_percentages = {}
        dimension_percentages = {}
        for prefix in model_dim:
            for permuted in permuteds:
                file_prefix = prefix + ("_permuted" if permuted else "")
                csv_path = os.path.join(dir, file_prefix + f"_{pathway_type}.csv")
                if os.path.isfile(csv_path):
                    df = pd.read_csv(csv_path, index_col=0)
                    df = fdr_filter(df, p, prefix)
                    pathway_counts[file_prefix] = count_elements(
                        pathway_nums[pathway_type],
                        df["pathway"].to_numpy().ravel().tolist(),
                    )
                    average_pathway_counts[file_prefix] = calculate_average_counts(
                        pathway_counts[file_prefix]
                    )
                    pathway_percentages[file_prefix] = (
                        len(set(df["pathway"].to_numpy().ravel().tolist()))
                        / pathway_nums[pathway_type]
                    )
                    dimension_counts[file_prefix] = count_elements(
                        model_dim[prefix], df["embed_idx"].to_numpy().ravel().tolist()
                    )
                    average_dimension_counts[file_prefix] = calculate_average_counts(
                        dimension_counts[file_prefix]
                    )
                    dimension_percentages[file_prefix] = (
                        len(set(df["embed_idx"].to_numpy().ravel().tolist()))
                        / model_dim[prefix]
                    )
                    rotated_file_prefix = file_prefix + "_rotated"
                    csv_path = os.path.join(
                        dir, rotated_file_prefix + f"_{pathway_type}.csv"
                    )
                    if os.path.isfile(csv_path):
                        df = pd.read_csv(csv_path, index_col=0)
                        df = fdr_filter(df, p, prefix)
                        pathway_counts[rotated_file_prefix] = count_elements(
                            pathway_nums[pathway_type],
                            df["pathway"].to_numpy().ravel().tolist(),
                        )
                        average_pathway_counts[rotated_file_prefix] = (
                            calculate_average_counts(
                                pathway_counts[rotated_file_prefix]
                            )
                        )
                        pathway_percentages[rotated_file_prefix] = (
                            len(set(df["pathway"].to_numpy().ravel().tolist()))
                            / pathway_nums[pathway_type]
                        )
                        dimension_counts[rotated_file_prefix] = count_elements(
                            model_dim[prefix],
                            df["embed_idx"].to_numpy().ravel().tolist(),
                        )
                        average_dimension_counts[rotated_file_prefix] = (
                            calculate_average_counts(
                                dimension_counts[rotated_file_prefix]
                            )
                        )
                        dimension_percentages[rotated_file_prefix] = (
                            len(set(df["embed_idx"].to_numpy().ravel().tolist()))
                            / model_dim[prefix]
                        )
        pathway_counts = align_counts(pathway_counts)
        dimension_counts = align_counts(dimension_counts)
        pathway_counts_df = pd.DataFrame(pathway_counts).T
        dimension_counts_df = pd.DataFrame(dimension_counts).T
        print(
            df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(pathway_counts) + 1)]),
                caption=f"Pathway Counts ({prefix} - {pathway_type})",
            )
        )
        print(
            dimension_counts_df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(dimension_counts) + 1)]),
                caption=f"Dimension Counts ({prefix} - {pathway_type})",
            )
        )
        pathway_counts_df.to_csv(
            os.path.join(dir, f"pathway_counts_{pathway_type}.csv"), index=True
        )
        dimension_counts_df.to_csv(
            os.path.join(dir, f"dimension_counts_{pathway_type}.csv"), index=True
        )
        average_pathway_counts_df = pd.DataFrame(average_pathway_counts, index=[0]).T
        average_dimension_counts_df = pd.DataFrame(
            average_dimension_counts, index=[0]
        ).T
        print(
            average_pathway_counts_df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(average_pathway_counts) + 1)]),
                caption=f"Average Pathway Counts ({pathway_type})",
            )
        )
        print(
            average_dimension_counts_df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(average_dimension_counts) + 1)]),
                caption=f"Average Dimension Counts ({pathway_type})",
            )
        )
        average_pathway_counts_df.to_csv(
            os.path.join(dir, f"average_pathway_counts_{pathway_type}.csv"), index=True
        )
        average_dimension_counts_df.to_csv(
            os.path.join(dir, f"average_dimension_counts_{pathway_type}.csv"),
            index=True,
        )
        pathway_percentages_df = pd.DataFrame(pathway_percentages, index=[0]).T
        print(
            pathway_percentages_df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(pathway_percentages) + 1)]),
                caption=f"Pathway Percentages ({pathway_type})",
            )
        )
        pathway_percentages_df.to_csv(
            os.path.join(dir, f"pathway_percentages_{pathway_type}.csv"), index=True
        )
        dimension_percentages_df = pd.DataFrame(dimension_percentages, index=[0]).T
        print(
            dimension_percentages_df.to_latex(
                column_format="|"
                + "".join(["c|" for _ in range(len(dimension_percentages) + 1)]),
                caption=f"Dimension Percentages ({pathway_type})",
            )
        )
        dimension_percentages_df.to_csv(
            os.path.join(dir, f"dimension_percentages_{pathway_type}.csv"), index=True
        )
        if pathway_type == "h":
            output["Mean Dimension Counts (H)"] = average_dimension_counts
            output["Dimension Coverage (H)"] = dimension_percentages
            output["Mean Pathway Counts (H)"] = average_pathway_counts
            output["Pathway Coverage (H)"] = pathway_percentages
        elif pathway_type == "c2":
            output["Mean Dimension Counts (C2)"] = average_dimension_counts
            output["Dimension Coverage (C2)"] = dimension_percentages
            output["Mean Pathway Counts (C2)"] = average_pathway_counts
            output["Pathway Coverage (C2)"] = pathway_percentages
    output_df = pd.DataFrame(output)
    output_df.to_csv(os.path.join(dir, "summary.csv"), index=True)


if "__main__" == __name__:
    import fire

    fire.Fire(main)
