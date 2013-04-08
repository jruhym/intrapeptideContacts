
import sys
import glob
from hydrogen_bonds import *
from bioinf import PDBATOMFileReader

for currentWildcard in sys.argv[1:]:

	for file_path in glob.glob(currentWildcard):

		reader = PDBATOMFileReader(file_path)
		atomsList = []
		donorList = []
		acceptorList = []
		for atom in reader:
			atomsList.append(atom)
			if atom.is_donor:
				donorList.append(atom)
			if atom.is_acceptor:
				acceptorList.append(atom)
		print len(atomsList)
for donor in donorList:

	for acceptor in acceptorList:

		