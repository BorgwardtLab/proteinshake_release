This repository contains the code for releasing proteinshake datasets.

To release, install `proteinshake`, `foldseek`, and `cd-hit`:

https://github.com/BorgwardtLab/proteinshake
https://github.com/steineggerlab/foldseek
https://github.com/weizhongli/cdhit

Adjust the `SCRATCH`, `RELEASE`, `NJOBS` variables in `release.py` according to your system.

Then run: `python release.py`

Note that this repository is licensed under GPLv3, due to dependency licenses.