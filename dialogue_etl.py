import csv
import requests
from Bio import Entrez
from semlib import Session
import os
from collections import Counter
import asyncio
'''
dialogue_etl.py contains for extracting, transforming, and loading the data from the patient-doctor conversations.
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

async def extract_reasons(convos):
    '''
    Returns list of extracted reasons for each patient visit.
    '''
    extracted_reasons = await session.map(
        convos,
        template=lambda r: f"""
        Extract the patient's chief complaint or reason for coming to the doctor. 1-2 word response only.
        {r['dialogue']}
        """.strip(),
    )
    print(f'Example reason: {extracted_reasons[0]}')

    return extracted_reasons

async def extract_family_illnesses(convos):
    '''
    Returns list of family illnesses for each patient visit.
    '''
    extracted_family_illnesses = await session.map(
        convos,
        template=lambda r: f"""
        Extract the any illness(es) that the patient and doctor are concerned about. Respond with only the illness(es) separated by commas. If none respond with none.
        {r['section_text']}
        """.strip(),
    )
    print(f'Example family illness: {extracted_family_illnesses[0]}')

    return extracted_family_illnesses

async def extract_symptoms(convos):
    '''
    Returns list of extracted symptoms for each patient visit.
    '''
    extracted_symptoms = await session.map(
        convos,
        template=lambda r: f"""
        Extract the patient's symptom(s). Respond with only the symptom(s) separated by commas. 1-2 words for each symptom. If patient has no symptons, respond with none.
        {r['section_text']}
        """.strip(),
    )
    print(f'Example symptom: {extracted_symptoms[0]}')

    return extract_symptoms