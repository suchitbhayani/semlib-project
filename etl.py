import csv
import os
import requests
from Bio import Entrez
'''
etl.py contains functions used to get datas.
'''

def download_dialogue():
    '''
    Downloads the abstracts. Saves data in data/MTS-Dialog-TrainingSet.csv
    '''
    url = "https://raw.githubusercontent.com/abachaa/MTS-Dialog/main/Main-Dataset/MTS-Dialog-TrainingSet.csv"
    
    save_path = os.path.join("data", "MTS-Dialog-TrainingSet.csv")
    
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "wb") as f:
            f.write(response.content)
        print(f"File downloaded and saved to {save_path}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

def load_dialogue():
    '''
    Returns list of convos. Each convo is a dictionary containing its ID, section header, section text, and dialogue.
    Must have data/MTS-Dialog-TrainingSet.csv first. Run download_dialogue() if you don't.
    '''
    with open("data/MTS-Dialog-TrainingSet.csv", encoding="latin-1") as f_in:
        csv_file = csv.reader(f_in)
        header = next(csv_file)
        convos = [dict(zip(header, row, strict=False)) for row in csv_file]

    print(f"Loaded {len(convos)} convos\n")

    return convos

def download_abstracts():
    '''
    Downloads the research papers. Saves data in data/PubMed_Conversations.csv.
    '''
    Entrez.email = os.getenv("EMAIL")

    # Define your PubMed query
    query = '"drug repurposing"[Title/Abstract] AND "Alzheimer"[Title/Abstract] AND ("2022"[Date - Publication] : "2025"[Date - Publication])'

    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed",
            term=query,
            usehistory="y"
        )
    )
    count = int(search_results["Count"])

    save_path = os.path.join("data", "PubMed_Conversations.csv")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    header = ["pmid", "title", "abstract"]
    batch_size = 20

    with open(save_path, "w", newline="", encoding="latin-1", errors="ignore") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)

        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)

            # Fetch data in XML for easy parsing
            handle = Entrez.efetch(
                db="pubmed",
                rettype="abstract",
                retmode="xml",
                retstart=start,
                retmax=batch_size,
                webenv=search_results["WebEnv"],
                query_key=search_results["QueryKey"]
            )

            records = Entrez.read(handle)
            handle.close()

            for article in records["PubmedArticle"]:
                pmid = article["MedlineCitation"]["PMID"]
                article_data = article["MedlineCitation"]["Article"]
                title = article_data.get("ArticleTitle", "")
                abstract_parts = article_data.get("Abstract", {}).get("AbstractText", [])
                abstract = " ".join(abstract_parts)
                
                writer.writerow([pmid, title, abstract])

    print(f"File downloaded and saved to {save_path}")

def load_abstracts():
    '''
    Returns list of abstracts. 
    Must have data/PubMed_Conversations.csv. Run download_abstracts() if you don't.
    '''
    with open("data/PubMed_Conversations.csv", encoding="latin-1") as f_in:
        csv_file = csv.reader(f_in)
        header = next(csv_file)
        abstracts = [dict(zip(header, row, strict=False)) for row in csv_file]

    print(f"Loaded {len(abstracts)} abstracts\n")
    return abstracts
