#!/usr/bin/env python3
import sys, os, re, tempfile
from pathlib import Path
from collections import defaultdict

FFMPEG = "/work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts/ffmpeg/ffmpeg"

def folderToMovies(inFolder, framerate=8):
    inFolder = Path(inFolder)
    if not inFolder.is_dir():
        print(f"Error: {inFolder} is not a directory")
        return

    files = sorted(inFolder.glob("*.png"))
    if not files:
        print(f"No PNG files found in {inFolder}")
        return

    # match all model types (case-insensitive)
    pattern = re.compile(
        r'(?i)output_(?P<state>Ma|Sa)[\+\-]\d+\.\d+_\d+_(?P<model>RBETA(?:wJET)?|CRITBETA(?:wJET)?)_(?P<pos>[\d\.]+)'
    )

    groups = defaultdict(list)
    for f in files:
        m = pattern.search(f.name)
        if m:
            state = "MAD" if m.group("state").lower() == "ma" else "SANE"
            model = m.group("model").upper()
            posval = "pos0" if m.group("pos").startswith("0.") else f"pos{m.group('pos')}"
            key = f"{state}_{model}_{posval}"
            groups[key].append(f)
        else:
            print(f"Skipping: {f.name}")

    for key, flist in groups.items():
        if not flist:
            continue

        outMovie = inFolder / f"{key}.avi"
        print(f"\n[+] Making movie for {key} with {len(flist)} frames...")

        # write a temporary file list for ffmpeg
        with tempfile.NamedTemporaryFile("w", delete=False) as tf:
            for f in sorted(flist):
                tf.write(f"file '{f}'\n")
            filelist = tf.name

        # ffmpeg concat input mode
        cmd = (
            f"{FFMPEG} -y -r {int(framerate)} -f concat -safe 0 -i {filelist} "
            f"-s:v 1280x720 -vcodec mjpeg -q:v 3 '{outMovie}'"
        )
        os.system(cmd)
        os.remove(filelist)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python folderToMovies.py <input_folder>")
        sys.exit(1)
    folderToMovies(sys.argv[1])
