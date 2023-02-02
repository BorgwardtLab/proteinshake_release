This repository contains the code for releasing proteinshake datasets.

It also outsources packages with incompatible licenses (cd-hit and MADOKA) which are used for sequence and structure clustering.

Clone the repo to Euler, create a venv with proteinshake, activate it, and run:

sh submit_parse.sh (few hours)
sh submit_sequence.sh (few minutes)
sh submit_structure.sh (several days)
sh submit_tasks.sh (few minutes)

This will submit a few thousand jobs on Euler and copy the result back to Borg. Wait until each script has finished before submitting the next one.
