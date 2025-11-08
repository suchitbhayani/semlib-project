from semlib import Session
import os
from collections import Counter
import asyncio

'''
analysis.py contains functions used to further ETL on the patient-dialgoue data.
'''

MAX_CONCURRENCY = 5
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
session = Session(model="openai/gpt-4.1-mini", max_concurrency=MAX_CONCURRENCY)

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