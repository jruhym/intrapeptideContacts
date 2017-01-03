
import sys
import glob
from hydrogen_bonds import *

for currentWildcard in sys.argv[1:]:
    for file_path in glob.glob(currentWildcard):
        reader = PDBATOMFileReader(file_path)
