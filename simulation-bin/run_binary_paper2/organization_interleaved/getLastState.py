# Script to copy the saved state from the folder with the newest timestamp
# (c) Jannik Luboeinski 2021

import numpy as np
from pathlib import Path
import shutil
import os

rawpaths = Path('.')
folders = sorted([str(x) for x in rawpaths.iterdir() if x.is_dir() and "_TRIPLET" in str(x)])
latest_folder = folders[len(folders)-1]
print("latest_folder: ", latest_folder)

os.remove("saved_state.txt") 

shutil.copyfile(latest_folder + "/saved_state.txt", "saved_state.txt")

