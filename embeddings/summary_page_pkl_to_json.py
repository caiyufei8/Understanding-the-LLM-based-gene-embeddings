import pickle, fire, json


def main():
    gene_name_to_summary_page = pickle.load(open("gene_name_to_summary_page.pkl", "rb"))
    about_to_delete = []
    for gene_name in gene_name_to_summary_page:
        if "" == gene_name_to_summary_page[gene_name]:
            about_to_delete.append(gene_name)
    for gene_name in about_to_delete:
        del gene_name_to_summary_page[gene_name]
    json.dump(
        gene_name_to_summary_page, open("gene_name_to_summary_page.json", "w"), indent=4
    )
    gene_name_to_summary_page_halved = pickle.load(
        open("gene_name_to_summary_page_halved.pkl", "rb")
    )
    about_to_delete = []
    for gene_name in gene_name_to_summary_page_halved:
        if "" == gene_name_to_summary_page_halved[gene_name]:
            about_to_delete.append(gene_name)
    for gene_name in about_to_delete:
        del gene_name_to_summary_page_halved[gene_name]
    json.dump(
        gene_name_to_summary_page_halved,
        open("gene_name_to_summary_page_halved.json", "w"),
        indent=4,
    )


if "__main__" == __name__:
    fire.Fire(main)
