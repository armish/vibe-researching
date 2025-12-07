#!/usr/bin/env python3
"""
Submit prompts to GPT-5 Pro with web search enabled using background mode.

This script reads prompt files from results/prompts/ and submits them to
OpenAI's GPT-5 Pro model with search enabled. Results are saved to
results/gpt5-pro-response/. Supports parallel execution and retry logic.

Uses OpenAI's background mode for reliable long-running queries without
timeout issues.

Usage:
    # First, install dependencies:
    pip install -r requirements.txt

    # Set your OpenAI API key:
    export OPENAI_API_KEY="your-api-key"

    # Dry run to see what would be processed:
    python submit_to_gpt.py --dry-run

    # Process all prompts:
    python submit_to_gpt.py

    # Process a single model:
    python submit_to_gpt.py --model-id ACH-000001

    # Reprocess all (including existing):
    python submit_to_gpt.py --overwrite

    # Adjust concurrency:
    python submit_to_gpt.py --concurrency 3
"""

import argparse
import asyncio
from pathlib import Path
from time import sleep

from openai import AsyncOpenAI

# Configuration
MAX_CONCURRENT_REQUESTS = 5  # Adjust based on your rate limits
MAX_RETRIES = 3
RETRY_DELAY_SECONDS = 10
POLL_INTERVAL_SECONDS = 5  # How often to check background task status

# Model configuration
MODEL = "gpt-5-pro"


async def submit_prompt(
    client: AsyncOpenAI,
    prompt_file: Path,
    output_dir: Path,
    overwrite: bool = False,
    semaphore: asyncio.Semaphore = None,
    poll_interval: int = POLL_INTERVAL_SECONDS,
) -> tuple[str, bool, str]:
    """
    Submit a single prompt to GPT-5 Pro using background mode and save the response.

    Returns:
        Tuple of (model_id, success, message)
    """
    model_id = prompt_file.stem
    output_file = output_dir / f"{model_id}.md"

    # Check if output already exists
    if output_file.exists() and not overwrite:
        return (model_id, True, "Skipped (already exists)")

    # Read the prompt
    with open(prompt_file, 'r') as f:
        prompt_content = f.read()

    async with semaphore:
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                # Submit to GPT-5 Pro with web search enabled using background mode
                # Background mode handles long-running tasks without timeout issues
                response = await client.responses.create(
                    model=MODEL,
                    input=prompt_content,
                    tools=[{"type": "web_search_preview"}],
                    background=True,
                )

                response_id = response.id
                print(f"  [{model_id}] Started background task: {response_id}")

                # Poll for completion
                while response.status in ("queued", "in_progress"):
                    await asyncio.sleep(poll_interval)
                    response = await client.responses.retrieve(response_id)
                    print(f"  [{model_id}] Status: {response.status}")

                # Check final status
                if response.status == "completed":
                    # Extract the response text
                    response_text = response.output_text

                    if not response_text:
                        raise ValueError("Empty response received from API")

                    # Save the response
                    with open(output_file, 'w') as f:
                        f.write(response_text)

                    return (model_id, True, f"Success (attempt {attempt})")

                elif response.status == "failed":
                    error_info = getattr(response, 'error', 'Unknown error')
                    raise Exception(f"Background task failed: {error_info}")

                elif response.status == "cancelled":
                    raise Exception("Background task was cancelled")

                else:
                    raise Exception(f"Unexpected status: {response.status}")

            except Exception as e:
                error_msg = str(e)
                if attempt < MAX_RETRIES:
                    print(f"  [{model_id}] Attempt {attempt} failed: {error_msg}. Retrying in {RETRY_DELAY_SECONDS}s...")
                    await asyncio.sleep(RETRY_DELAY_SECONDS)
                else:
                    return (model_id, False, f"Failed after {MAX_RETRIES} attempts: {error_msg}")


async def main():
    parser = argparse.ArgumentParser(
        description="Submit prompts to GPT-5 Pro with web search enabled."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing response files",
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
        help="List prompts that would be processed without submitting",
    )
    parser.add_argument(
        "--poll-interval",
        type=int,
        default=POLL_INTERVAL_SECONDS,
        help=f"Polling interval in seconds for background tasks (default: {POLL_INTERVAL_SECONDS})",
    )
    args = parser.parse_args()

    # Define paths
    base_dir = Path(__file__).parent
    prompts_dir = base_dir / "results" / "prompts"
    output_dir = base_dir / "results" / "gpt5-pro-response"

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get list of prompt files
    if args.model_id:
        prompt_files = [prompts_dir / f"{args.model_id}.txt"]
        if not prompt_files[0].exists():
            print(f"Error: Prompt file not found: {prompt_files[0]}")
            return
    else:
        prompt_files = sorted(prompts_dir.glob("*.txt"))

    if not prompt_files:
        print("No prompt files found. Run generate_prompts.py first.")
        return

    print(f"Found {len(prompt_files)} prompt(s) to process.")
    print(f"Output directory: {output_dir}")
    print(f"Concurrency: {args.concurrency}")
    print(f"Poll interval: {args.poll_interval}s")
    print(f"Overwrite: {args.overwrite}")
    print()

    # Dry run mode
    if args.dry_run:
        print("Dry run mode - prompts that would be processed:")
        for pf in prompt_files:
            output_file = output_dir / f"{pf.stem}.md"
            status = "SKIP (exists)" if output_file.exists() and not args.overwrite else "PROCESS"
            print(f"  [{status}] {pf.stem}")
        return

    # Initialize OpenAI client
    # The API key should be set via OPENAI_API_KEY environment variable
    client = AsyncOpenAI()

    # Create semaphore for rate limiting
    semaphore = asyncio.Semaphore(args.concurrency)

    # Submit all prompts concurrently
    print("Submitting prompts to GPT-5 Pro (background mode)...")
    tasks = [
        submit_prompt(client, pf, output_dir, args.overwrite, semaphore, args.poll_interval)
        for pf in prompt_files
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
        print("\nSome prompts failed. Re-run the script to retry failed prompts.")
        print("Use --overwrite to reprocess all prompts including successful ones.")


if __name__ == "__main__":
    asyncio.run(main())
