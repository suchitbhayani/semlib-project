from semlib import Session
import os
from collections import Counter
import asyncio

MAX_CONCURRENCY = 5
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
session = Session(model="openai/gpt-4.1-mini", max_concurrency=MAX_CONCURRENCY)

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