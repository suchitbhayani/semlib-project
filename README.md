# AI-drug-discovery

This project provides scripts to **download and process patient-doctor dialogue and medical paper abstract datasets** related to drug discovery. The code automatically stores CSV files in the `data/` directory and performs ETL (Extract, Transform, Load) operations on them.  

---

## Setup

1. **Clone the repository**:  
```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
cd YOUR_REPO
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
## Usage 

Run the main script with **targets** to specify what to download and process.

### Targets

Approach 1:
- **`dialogue`** – Downloads the patient-doctor dialogue CSV and loads conversations.  
- **`reasons`** – Extracts visit reasons from dialogue data (**requires `dialogue`**).  
- **`illnesses`** – Extracts family illnesses from dialogue data (**requires `dialogue`**).  
- **`symptoms`** – Extracts symptoms from dialogue data (**requires `dialogue`**).

Approach 2:
- **`abstracts`** – Downloads abstracts CSV and extracts drug repurposing candidates.  

### Examples

1. **Download and process dialogues only:**
```bash
python run.py dialogue
```

2. **Download and process dialogues, and extract patient reasons for visits :**
```bash
python run.py dialogue reasons
```
When passing `dialogue`, can also pass `illnesses` and/or `symptoms` to extract family illnesses/symptoms from the dialogue csv data.

3. **Download abstracts and extract repurposing candidates:**
```bash
python run.py abstracts
```
