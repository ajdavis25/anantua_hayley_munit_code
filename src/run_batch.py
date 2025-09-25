import os, sys, logging
from config import DUMPS, CSV_SUMMARY
from munit_loader import load_munit_params
from render import process_dump_file
from summary_writer import append_to_summary

logging.basicConfig(level=logging.INFO)

def run():
    munit_params = load_munit_params()
    dump_files = sorted(f for f in os.listdir(DUMPS) if f.endswith(".h5"))

    if "SLURM_ARRAY_TASK_ID" not in os.environ:
        logging.error("SLURM_ARRAY_TASK_ID not set")
        sys.exit(1)

    index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index >= len(dump_files):
        logging.error("Index out of range")
        sys.exit(1)

    dump_file = dump_files[index]
    result_rows = process_dump_file(dump_file, munit_params)

    if result_rows:
        CSV_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
        append_to_summary(result_rows, str(CSV_SUMMARY))
        logging.info(f"Appended results for {dump_file} to summary.")

if __name__ == "__main__":
    run()
