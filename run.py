#!/usr/bin/env python
from dialogue_etl import *
from abstract_etl import *
import asyncio
import sys

async def main(targets):
    if 'dialogue' in targets:
        download_dialogue()
        convos = load_dialogue()

    # later on do something with these results  
    if 'reasons' in targets:
        visit_reasons = await extract_reasons(convos)
    if 'illnesses' in targets:
        family_illnesses = await extract_reasons(convos)
    if 'symptoms' in targets:
        symptoms = await extract_symptoms(convos)
    
    if 'abstracts' in targets:
        download_abstracts()
        abstracts = load_abstracts()

        # later on do something with this
        repurposing_candidates = await extract_candidates(abstracts)

if __name__ == '__main__':
    targets = sys.argv[1:]
    asyncio.run(main(targets))
