import csv, fcntl
from typing import List, Dict, Union

def append_to_summary(rows: List[Dict[str, Union[int, str, float, None]]], path: str):
    """
    appends a list of dictionaries (rows) to a CSV file in a thread-safe way
    using file locking (via `fcntl`) to prevent parallel write conflicts
    """
    with open(path, "a", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        if f.tell() == 0:
            writer.writeheader()
        writer.writerows(rows)
        fcntl.flock(f, fcntl.LOCK_UN)
