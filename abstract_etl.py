import csv
import requests
from Bio import Entrez
from semlib import Session
import os
from collections import Counter
import asyncio
'''
abstract_etl.py contains for extracting, transforming, and loading the drug repurposing candidates from the abstracts.
'''

MAX_CONCURRENCY = 5
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
session = Session(model="openai/gpt-4.1-mini", max_concurrency=MAX_CONCURRENCY)

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
    Returns structured list of abstracts. 
    Must have data/PubMed_Conversations.csv. Run download_abstracts() if you don't.
    '''
    try:
        with open("data/PubMed_Conversations.csv", encoding="latin-1") as f_in:
            csv_file = csv.reader(f_in)
            header = next(csv_file)
            abstracts = [dict(zip(header, row, strict=False)) for row in csv_file]

        print(f"Loaded {len(abstracts)} abstracts\n")
        return abstracts

    except FileNotFoundError:
        print("Must have data/PubMed_Conversations.csv. Run download_abstracts() if you don't.")
        return []

async def extract_candidates(abstracts):
    extracted_candidates = await session.map(
        abstracts,
        template=lambda r: f"""
        Extract the drug repurposing candidates. Repond only with the candidates, separated by a comma and space. If none were mentioned, respond with 'none'
        {r['abstract']}
        """.strip(),
    )
    print(f'Example candidates: {extracted_candidates[0]}')

    return extracted_candidates