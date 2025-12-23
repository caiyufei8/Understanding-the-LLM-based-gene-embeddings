import requests, html2text, mygene, os, pickle, openai, time
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd

mg = mygene.MyGeneInfo()
parts_to_remove = [
    "##  Summary\n",
    "NEW",
    "Try the newGene table",
    "Try the newTranscript table",
    "**",
    "\nGo to the top of the page Help\n",
]


def rough_text_from_gene_name(gene_number):
    # get url
    url = f"https://www.ncbi.nlm.nih.gov/gene/{gene_number}"
    # Send a GET request to the URL
    summary_text = ""
    soup = None
    try:
        response = requests.get(url, timeout=30)
    except requests.exceptions.Timeout:
        print("time out")
        return (summary_text, soup)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the HTML content of the page
        soup = BeautifulSoup(response.content, "html.parser")
        # Find the "summary" tab content by inspecting the page's structure
        summary_tab = soup.find("div", {"class": "rprt-section gene-summary"})
        # Check if the "summary" tab content is found
        if summary_tab:
            # Convert the HTML to plain text using html2text
            html_to_text = html2text.HTML2Text()
            html_to_text.ignore_links = True  # Ignore hyperlinks
            # Extract the plain text from the "summary" tab
            summary_text = html_to_text.handle(str(summary_tab))
            # Remove the specified parts from the original text
            for part in parts_to_remove:
                summary_text = summary_text.replace(part, " ")
                # Replace '\n' with a space
            summary_text = summary_text.replace("\n", " ")
            # Reduce multiple spaces into one space
            summary_text = " ".join(summary_text.split())
            # Print or save the extracted text
        else:
            print("Summary tab not found on the page.")
    else:
        print(f"Failed to retrieve the webpage. Status code: {response.status_code}")
    return (summary_text, soup)


GPT_DIM = 1536
openai.api_key = ""
models = (
    "text-embedding-3-small",
    "text-embedding-3-large",
    "text-embedding-ada-002",
)


def get_gpt_embedding(text, model="text-embedding-ada-002"):
    text = text.replace("\n", " ")
    while True:
        try:
            embedding = openai.Embedding.create(input=[text], model=model)["data"][0][
                "embedding"
            ]
            break
        except Exception as e:
            print("Error:", e)
            time.sleep(1)
    return np.array(embedding)


def find_median_length(gene_name_to_summary_page):
    lengths = []
    for gene_name in gene_name_to_summary_page:
        lengths.append(gene_name_to_summary_page[gene_name].count(" ") + 1)
    lengths.sort()
    if len(lengths) % 2 == 0:
        median_length = (
            lengths[len(lengths) // 2] + lengths[len(lengths) // 2 + 1]
        ) // 2
    else:
        median_length = lengths[len(lengths) // 2 + 1]
    return median_length


def cut_text(text, length):
    words = text.split(" ")
    if len(words) <= length:
        return text
    else:
        return " ".join(words[:length]) + "..."


def main():
    with open("gene_info_table.csv", "r") as f:
        gene_info_table = pd.read_csv(f, header=0, index_col=0)
    print(gene_info_table.head())
    if os.path.isfile("gene_name_to_tax_id.pkl"):
        gene_name_to_tax_id = pickle.load(open("gene_name_to_tax_id.pkl", "rb"))
    else:
        gene_name_to_tax_id = {}
        gene_names = gene_info_table.loc[:, "gene_name"].to_list()
        names = mg.querymany(gene_names, scopes="symbol", species="human")
        for result in names:
            if "_id" in result and "query" in result:
                gene_name_to_tax_id[result["symbol"]] = result["_id"]
        pickle.dump(gene_name_to_tax_id, open("gene_name_to_tax_id.pkl", "wb"))
    print(gene_name_to_tax_id)
    if os.path.isfile("gene_name_to_summary_page.pkl"):
        gene_name_to_summary_page = pickle.load(
            open("gene_name_to_summary_page.pkl", "rb")
        )
    else:
        gene_name_to_summary_page = {}
    for gene_name, page_id in sorted(gene_name_to_tax_id.items()):
        if gene_name not in gene_name_to_summary_page:
            print("gene_name", gene_name)
            parsed_text, unparsed_html = rough_text_from_gene_name(page_id)
            gene_name_to_summary_page[gene_name] = parsed_text
            pickle.dump(
                gene_name_to_summary_page,
                open("gene_name_to_summary_page.pkl", "wb"),
            )
    # median_length = find_median_length(gene_name_to_summary_page)
    # print("Median Length:", median_length)
    if os.path.isfile("gene_name_to_summary_page_halved.pkl"):
        gene_name_to_summary_page_halved = pickle.load(
            open("gene_name_to_summary_page_halved.pkl", "rb")
        )
    else:
        gene_name_to_summary_page_halved = {}
    for gene_name in sorted(gene_name_to_summary_page):
        if gene_name not in gene_name_to_summary_page_halved:
            gene_name_to_summary_page_halved[gene_name] = cut_text(
                gene_name_to_summary_page[gene_name],
                (gene_name_to_summary_page[gene_name].count(" ") + 1) // 2,
            )
    pickle.dump(
        gene_name_to_summary_page_halved,
        open("gene_name_to_summary_page_halved.pkl", "wb"),
    )
    print(gene_name_to_summary_page_halved)
    if os.path.isfile("gene_name_to_embedding_half_length.pkl"):
        gene_name_to_embedding = pickle.load(
            open("gene_name_to_embedding_half_length.pkl", "rb")
        )
    else:
        gene_name_to_embedding = {}
    try:
        for model in models:
            if model not in gene_name_to_embedding:
                gene_name_to_embedding[model] = {}
            for key, text in sorted(gene_name_to_summary_page.items()):
                if key not in gene_name_to_embedding[model]:
                    print("key", key)
                    if text != "":
                        gene_name_to_embedding[model][key] = get_gpt_embedding(
                            gene_name_to_summary_page_halved[key], model
                        )
    except Exception as e:
        pickle.dump(
            gene_name_to_embedding, open("gene_name_to_embedding_half_length.pkl", "wb")
        )
        print(e)
        return
    pickle.dump(
        gene_name_to_embedding, open("gene_name_to_embedding_half_length.pkl", "wb")
    )
    print(gene_name_to_embedding)


if "__main__" == __name__:
    main()
