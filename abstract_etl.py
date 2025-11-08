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
    Returns structured list of conversations.
    Must have data/MTS-Dialog-TrainingSet.csv first. Run download_dialogue() if you don't.
    '''
    try:
        with open("data/MTS-Dialog-TrainingSet.csv", encoding="latin-1") as f_in:
            csv_file = csv.reader(f_in)
            header = next(csv_file)
            convos = [dict(zip(header, row, strict=False)) for row in csv_file]

        print(f"Loaded {len(convos)} conversations\n")
        return convos

    except FileNotFoundError:
        print("Must have data/MTS-Dialog-TrainingSet.csv first. Run download_dialogue() if you don't.")
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