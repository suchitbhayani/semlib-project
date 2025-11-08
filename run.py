#!/usr/bin/env python
from dialogue_etl import *
from abstract_etl import *
import asyncio
import sys

async def dialogue_approach(targets):
    '''
    Main logic to handle patient-doctor conversation dialogue approach.
    '''
    if 'download_dialogue' in targets:
        download_dialogue()

    all_subargs = ['reasons', 'illnesses', 'symptoms']
    if any(arg in targets for arg in all_subargs):
        convos = load_dialogue()

    # later on do something with these results  
    if 'reasons' in targets:
        visit_reasons = await extract_reasons(convos)
    if 'illnesses' in targets:
        family_illnesses = await extract_reasons(convos)
    if 'symptoms' in targets:
        symptoms = await extract_symptoms(convos)

async def abstracts_approach(targets):
    '''
    Main logic to handle research abstract searching for drug repurposing approach.
    '''
    if 'download_abstracts' in targets:
        download_abstracts()

    all_subargs = ['repurposing']
    if any(arg in targets for arg in all_subargs):
        abstracts = load_abstracts()

    # later on do something with this
    if 'repurposing' in targets:
        repurposing_candidates = await extract_candidates(abstracts)

async def main(targets):
    await dialogue_approach(targets)
    await abstracts_approach(targets)

if __name__ == '__main__':
    targets = sys.argv[1:]
    asyncio.run(main(targets))
