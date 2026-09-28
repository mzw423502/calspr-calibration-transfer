#!/usr/bin/env python3
from make_v24_figures import main
import sys
if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("usage: run_v24_figures.py CONFIG")
    main(sys.argv[1])
