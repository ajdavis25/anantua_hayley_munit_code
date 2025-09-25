import matplotlib
matplotlib.use("Agg")

from pathlib import Path
from .creating_images import plotPositronTestFrame  # Or use relative import if placed inside /src/plot
from tqdm import tqdm  # Optional: for progress bar


def find_all_h5_files(root_dir: Path):
    return list(root_dir.glob("**/*.h5"))


def main():
    root = Path("/work/vmo703/ipole_outputs") # change if needed

    h5_files = find_all_h5_files(root)
    print(f"Found {len(h5_files)} .h5 files")

    for file in tqdm(h5_files, desc="Rendering"):
        png_file = file.with_suffix(".png")
        if png_file.exists():
            continue

        try:
            plotPositronTestFrame(
                str(file),
                cpMax=0.1,
                fractionalCircular=False
            )
        except Exception as e:
            print(f"Error processing {file.name}: {e}")


if __name__ == "__main__":
    main()
