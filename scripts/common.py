"""Shared exact-data reader; accepts plain JSONL or gzip-compressed JSONL."""
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "data" / "exact_results.jsonl.gz"

def read_rows(path=DEFAULT_DATA):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)
