import sys
import os

def folderToMovie(inFolder, outMovie, framerate=8, patternMatch="%*.png"):
    if inFolder[-1] != '/':
        inFolder += '/'
    os.system("/work/epc964/ipole_HarvardCannonHPC_copy/aricarte/ipole+e-/ffmpeg/ffmpeg -y -framerate "+str(int(framerate))+" -i "+inFolder+patternMatch+" -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p "+outMovie)#Added path /work/epc964/ipole_HarvardCannonHPC_copy/aricarte/ipole+e-/ffmpeg 11/29/22

if __name__ == '__main__':
    inFolder = sys.argv[1]
    outMovie = sys.argv[2]
    folderToMovie(inFolder, outMovie)
