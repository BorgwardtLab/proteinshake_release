This repository contains the code for releasing proteinshake datasets.

It also outsources packages with incompatible licenses (cd-hit and MADOKA) which are used for sequence and structure clustering.

Clone the repo to Euler, create a venv with proteinshake, activate it, and run:

sbatch --wrap "sh submit.sh" -o logs/wrapper.txt

This will submit a few thousand jobs on Euler and copy the result back to Borg.
