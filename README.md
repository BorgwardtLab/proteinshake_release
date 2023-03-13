This repository contains the code for releasing proteinshake datasets. It also computes the random, sequence, and structure splits.

To release, install `proteinshake`, `foldseek`, and `cd-hit`:

- https://github.com/BorgwardtLab/proteinshake
- https://github.com/steineggerlab/foldseek
- https://github.com/weizhongli/cdhit

Adjust the `TAG`, `SCRATCH`, `DESTINATION`, `NJOBS` variables in `release.py` according to your system.

Then run: `python release.py`

Note that this repository is licensed under GPLv3, due to dependency licenses.
