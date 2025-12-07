#!/usr/bin/env python3
"""
Generate prompt files for GPT-5 Pro queries based on cell line data.

This script reads the results table, mutations data, and top 10 dependencies
to create prompt files for each cell line model.
"""

import csv
import json
import os
from pathlib import Path


def read_tsv(filepath: str) -> list[dict]:
    """Read a TSV file and return a list of dictionaries."""
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return list(reader)


def read_mutations(model_id: str, mutations_dir: str) -> str:
    """
    Read mutation CSV for a model and return a tab-delimited table
    with Gene, Variant Type, Variant Info, and Protein Change columns.
    """
    filepath = Path(mutations_dir) / f"{model_id}.csv"
    if not filepath.exists():
        return "No mutation data available."

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        rows = []
        # Header
        rows.append("Gene\tVariant Type\tVariant Info\tProtein Change")
        for row in reader:
            gene = row.get('Gene', '')
            variant_type = row.get('Variant Type', '')
            variant_info = row.get('Variant Info', '')
            protein_change = row.get('Protein Change', '')
            rows.append(f"{gene}\t{variant_type}\t{variant_info}\t{protein_change}")

    return '\n'.join(rows)


def read_top10_dependencies(model_id: str, top10deps_dir: str) -> str:
    """
    Read top 10 dependencies JSON for a model and return a comma-separated
    string of gene symbols from the 'labels' array.
    """
    filepath = Path(top10deps_dir) / f"{model_id}.json"
    if not filepath.exists():
        return "No dependency data available."

    with open(filepath, 'r') as f:
        data = json.load(f)

    labels = data.get('labels', [])
    return ', '.join(labels)


def generate_prompt(cell_line_name: str, gene: str, top10deps: str, mutation_table: str) -> str:
    """Generate the prompt from the template."""
    template = f"""A cancer cell line, called {cell_line_name}, is extremely sensitive to gene {gene} depletion via CRISPR-KO. Top 10 major dependencies for this cell line are: {top10deps}.

Here are the mutations identified in this cell line:
{mutation_table}

Based on these, can you help me understand why {cell_line_name} might be super sensitive to {gene} depletion? Use markdown when you are answering and make sure you keep the references in."""

    return template


def main():
    # Define paths
    base_dir = Path(__file__).parent
    results_tsv = base_dir / 'results' / 'n1s_tbl.tsv'
    mutations_dir = base_dir / 'data' / 'mutations'
    top10deps_dir = base_dir / 'data' / 'top10deps'
    prompts_dir = base_dir / 'results' / 'prompts'

    # Create prompts directory if it doesn't exist
    prompts_dir.mkdir(parents=True, exist_ok=True)

    # Read the results table
    rows = read_tsv(results_tsv)
    print(f"Found {len(rows)} cell line models to process.")

    # Process each row
    for row in rows:
        model_id = row['ModelID']
        cell_line_name = row['StrippedCellLineName']
        gene = row['gene']

        # Read mutation data
        mutation_table = read_mutations(model_id, mutations_dir)

        # Read top 10 dependencies
        top10deps = read_top10_dependencies(model_id, top10deps_dir)

        # Generate prompt
        prompt = generate_prompt(cell_line_name, gene, top10deps, mutation_table)

        # Write prompt to file
        prompt_file = prompts_dir / f"{model_id}.txt"
        with open(prompt_file, 'w') as f:
            f.write(prompt)

        print(f"Generated prompt for {model_id} ({cell_line_name})")

    print(f"\nDone! Prompts saved to {prompts_dir}")


if __name__ == '__main__':
    main()
