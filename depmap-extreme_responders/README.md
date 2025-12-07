# DepMap Extreme Responder workflow

- Download missing data files (`CRISPRGeneEffect.csv` and `OmicsCNGeneWGS.csv`) from the corresponding url files under `data/`
- Run the R notebook `Identify extreme responders.ipynb` 
- Fetch mutation and top 10 dependencies:
  - `cd data/mutations && bash download.sh`
  - `cd data/top10deps && bash download.sh`
- Generate prompts: `python generate_prompts.py`
- Set your `OPENAI_API_KEY`: `export OPENAI_API_KEY='...'`
- Submit prompts to `gpt-5-pro`: `python submit_to_gpt.py` (caution, each query costs ~$10)
- Summarize responses with `gpt-5-nano`: `python summarize_responses.py`

All results are under `results/`.
