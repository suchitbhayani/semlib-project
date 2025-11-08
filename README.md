# AI-drug-discovery

This project provides scripts to **download and process patient-doctor dialogue and medical paper abstract datasets** related to drug discovery. The code automatically stores CSV files in the `data/` directory and performs ETL (Extract, Transform, Load) operations on them.  

---

## Setup

1. **Clone the repository**:  
```bash
git clone https://github.com/suchitbhayani/AI-drug-discovery.git
cd AI_drug_discovery
```

2. **Create and activate a virtual environment**

```bash
# Create a virtual environment named 'env'
python -m venv env

# Activate on Windows Git Bash
source env/Scripts/activate
```

3. **Install dependencies**

```bash
pip install -r requirements.txt
```

4. **Create a .env file at the root of the project containing your API key and email**

```ini
OPENAI_API_KEY=your_openai_api_key_here
EMAIL=your_email_here
```

---

## File organization
`abstract` directory has code and parameter configurations for the research abstract data approach. You can configure parameters in `abstract/abstract_params.json`.

`dialogue` directory has code for the patient-doctor conversation data approach. There are currently no parameters to configure for this approach.

---

## Usage 

Run the main script with **targets** to specify what to download and process.

### Targets

Approach 1:
- **`download_dialogue`** – Downloads the patient-doctor dialogue CSV.  
- **`reasons`** – Extracts visit reasons from dialogue data (**requires `python run.py download_dialogue` to have been run before to download data**).  
- **`illnesses`** – Extracts family illnesses from dialogue data (**requires `python run.py download_dialogue` to have been run before to download data**).  
- **`symptoms`** – Extracts symptoms from dialogue data (**requires `python run.py download_dialogue` to have been run before to load downdata**).

Approach 2:
- **`download_abstracts`** – Downloads drug repurposing candidate CSV.  
- **`repurposing`** – Extracts drug repurposing candidates data (**requires `python run.py download_abstracts` to have been run before to download data**).  

Run these by running `py run.py [any combination of targets]`

### Examples

1. **Downloads dialogues data:**
```bash
python run.py download_dialogue
```

2. **Extracts patient reasons for visits, family illnesses, and symptoms from patient-doctor dialogue data:**
```bash
python run.py reasons illnesses symptoms
```
When passing do not have to pass all three of `reasons`, `illnesses` and `symptoms`. Can choose to pass as many or as little as you want.

3. **Downloads research paper abstract data. Can configure parameters at `abstract/abstract_params.json`:**
```bash
python run.py download_abstracts
```

4. **Extracts drug repurposing candidate data from research paper abstract data:**
```bash
python run.py repurposing
```
