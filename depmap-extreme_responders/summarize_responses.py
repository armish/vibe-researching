#!/usr/bin/env python3
"""
Summarize GPT-5 Pro responses using GPT-5 Nano.

This script reads the GPT-5 Pro responses and submits them to GPT-5 Nano
for summarization into three structured paragraphs.

Usage:
    # First, install dependencies:
    pip install -r requirements.txt

    # Set your OpenAI API key:
    export OPENAI_API_KEY="your-api-key"

    # Dry run to see what would be processed:
    python summarize_responses.py --dry-run

    # Process all responses:
    python summarize_responses.py

    # Process a single model:
    python summarize_responses.py --model-id ACH-000001

    # Reprocess all (including existing):
    python summarize_responses.py --overwrite

    # Adjust concurrency:
    python summarize_responses.py --concurrency 10
"""

import argparse
import asyncio
from pathlib import Path

from openai import AsyncOpenAI

# Configuration
MAX_CONCURRENT_REQUESTS = 10  # Nano is faster, can handle more concurrency
MAX_RETRIES = 3
RETRY_DELAY_SECONDS = 5

# Model configuration
MODEL = "gpt-5-nano"

PROMPT_TEMPLATE = """Below is a set of explanations why a cell line might be super sensitive to a particular gene depletion via CRISPR KO. I want to summarize it in three paragraphs: the first one should explain the most reasonable and lead explanation; the second one should talk about if there are factors that are not as strong as the initial one but those that might stil contribute to the sensitivy; the third should list other interesting facts. Each paragragh can be 5-6 sentences.

Here are the explanations: {gpt5_pro_response}.

Provide your answer in plain text."""


async def summarize_response(
    client: AsyncOpenAI,
    response_file: Path,
    output_dir: Path,
    overwrite: bool = False,
    semaphore: asyncio.Semaphore = None,
) -> tuple[str, bool, str]:
    """
    Summarize a single GPT-5 Pro response using GPT-5 Nano.

    Returns:
        Tuple of (model_id, success, message)
    """
    model_id = response_file.stem
    output_file = output_dir / f"{model_id}.txt"

    # Check if output already exists
    if output_file.exists() and not overwrite:
        return (model_id, True, "Skipped (already exists)")

    # Read the GPT-5 Pro response
    with open(response_file, 'r') as f:
        gpt5_pro_response = f.read()

    # Build the prompt
    prompt = PROMPT_TEMPLATE.format(gpt5_pro_response=gpt5_pro_response)

    async with semaphore:
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                # Submit to GPT-5 Nano (faster model, no background mode needed)
                response = await client.responses.create(
                    model=MODEL,
                    input=prompt,
                )

                # Extract the response text
                response_text = response.output_text

                if not response_text:
                    raise ValueError("Empty response received from API")

                # Save the response
                with open(output_file, 'w') as f:
                    f.write(response_text)

                return (model_id, True, f"Success (attempt {attempt})")

            except Exception as e:
                error_msg = str(e)
                if attempt < MAX_RETRIES:
                    print(f"  [{model_id}] Attempt {attempt} failed: {error_msg}. Retrying in {RETRY_DELAY_SECONDS}s...")
                    await asyncio.sleep(RETRY_DELAY_SECONDS)
                else:
                    return (model_id, False, f"Failed after {MAX_RETRIES} attempts: {error_msg}")


async def main():
    parser = argparse.ArgumentParser(
        description="Summarize GPT-5 Pro responses using GPT-5 Nano."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing summary files",
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=MAX_CONCURRENT_REQUESTS,
        help=f"Maximum concurrent requests (default: {MAX_CONCURRENT_REQUESTS})",
    )
    parser.add_argument(
        "--model-id",
        type=str,
        help="Process only a specific model ID (for testing)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List responses that would be processed without submitting",
    )
    args = parser.parse_args()

    # Define paths
    base_dir = Path(__file__).parent
    responses_dir = base_dir / "results" / "gpt5-pro-response"
    output_dir = base_dir / "results" / "gpt5-nano-summary"

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get list of response files
    if args.model_id:
        response_files = [responses_dir / f"{args.model_id}.md"]
        if not response_files[0].exists():
            print(f"Error: Response file not found: {response_files[0]}")
            return
    else:
        response_files = sorted(responses_dir.glob("*.md"))

    if not response_files:
        print("No response files found. Run submit_to_gpt.py first.")
        return

    print(f"Found {len(response_files)} response(s) to summarize.")
    print(f"Output directory: {output_dir}")
    print(f"Concurrency: {args.concurrency}")
    print(f"Overwrite: {args.overwrite}")
    print()

    # Dry run mode
    if args.dry_run:
        print("Dry run mode - responses that would be processed:")
        for rf in response_files:
            output_file = output_dir / f"{rf.stem}.txt"
            status = "SKIP (exists)" if output_file.exists() and not args.overwrite else "PROCESS"
            print(f"  [{status}] {rf.stem}")
        return

    # Initialize OpenAI client
    # The API key should be set via OPENAI_API_KEY environment variable
    client = AsyncOpenAI()

    # Create semaphore for rate limiting
    semaphore = asyncio.Semaphore(args.concurrency)

    # Submit all responses concurrently
    print("Summarizing responses with GPT-5 Nano...")
    tasks = [
        summarize_response(client, rf, output_dir, args.overwrite, semaphore)
        for rf in response_files
    ]

    results = await asyncio.gather(*tasks)

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    success_count = 0
    skip_count = 0
    fail_count = 0

    for model_id, success, message in results:
        if "Skipped" in message:
            skip_count += 1
            print(f"  [SKIP] {model_id}: {message}")
        elif success:
            success_count += 1
            print(f"  [OK]   {model_id}: {message}")
        else:
            fail_count += 1
            print(f"  [FAIL] {model_id}: {message}")

    print()
    print(f"Total: {len(results)} | Success: {success_count} | Skipped: {skip_count} | Failed: {fail_count}")

    if fail_count > 0:
        print("\nSome summaries failed. Re-run the script to retry failed ones.")
        print("Use --overwrite to reprocess all summaries including successful ones.")


if __name__ == "__main__":
    asyncio.run(main())
