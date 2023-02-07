This repository contains the code for releasing proteinshake datasets.

To release, install `proteinshake`, `TMalign`, and `cd-hit`:

https://github.com/BorgwardtLab/proteinshake
https://zhanggroup.org/TM-align/TMalign.cpp
https://github.com/weizhongli/cdhit

Adjust the `SCRATCH`, `RELEASE`, `NJOBS` variables in `release.py` according to your system.

Then run: `python release.py`