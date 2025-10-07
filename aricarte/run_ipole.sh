#!/usr/bin/env bash
export LD_LIBRARY_PATH=/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:$LD_LIBRARY_PATH
exec /work/vmo703/aricarte/aricarte-copy/aricarte/ipole+e-/ipole "$@"
