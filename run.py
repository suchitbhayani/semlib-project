#!/usr/bin/env python
from utils import *
from dialogue.dialogue_etl import *
from abstract.abstract_etl import *
import asyncio
import sys
import json

async def dialogue_approach(targets, run_all):
    '''
    Main logic to handle patient-doctor conversation dialogue approach.
    '''
    if run_all or 'download_dialogue' in targets:
        download_dialogue()

    all_subargs = ['reasons', 'illnesses', 'symptoms']
    if run_all or any(arg in targets for arg in all_subargs):
        convos = preprocess_dialogue()

    # later on do something with these results  
    if run_all or 'reasons' in targets:
        visit_reasons = await extract_reasons(convos)
    if run_all or 'illnesses' in targets:
        family_illnesses = await extract_reasons(convos)
    if run_all or 'symptoms' in targets:
        symptoms = await extract_symptoms(convos)

async def abstracts_approach(targets, run_all):
    '''
    Main logic to handle research abstract searching for drug repurposing approach.
    '''
    with open('abstract/abstract_params.json') as fh:
        data_params = json.load(fh)
        
    if run_all or 'download_abstracts' in targets:
        download_abstracts(**data_params)

    all_subargs = ['repurposing']
    if run_all or any(arg in targets for arg in all_subargs):
        abstracts = preprocess_abstracts(**data_params)

    # later on do something with this
    if run_all or 'repurposing' in targets:
        repurposing_candidates = await extract_candidates(abstracts)

async def main(targets):
    if 'clean' in targets:
        clean_files()

    run_all = 'all' in targets

    await dialogue_approach(targets, run_all)
    await abstracts_approach(targets, run_all)

if __name__ == '__main__':
    targets = sys.argv[1:]
    asyncio.run(main(targets))
