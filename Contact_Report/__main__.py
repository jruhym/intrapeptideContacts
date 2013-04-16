
import sys
import glob
from hydrogen_bonds import *
#from bioinf import PDBATOMFileReader

for currentWildcard in sys.argv[1:]:

	for file_path in glob.glob(currentWildcard):

		reader = PDBATOMFileReader(file_path)
		atomDict = {}
		donorDict = {}
		acceptorDict = {}
		for atom in reader:
			atomDict[atom.serial] = atom 
			if atom.is_donor:
				donorDict[atom.serial] = atom
			if atom.is_acceptor:
				acceptorDict[atom.serial] = atom
		print len(atomsList)
for donor in donorList:

	for acceptor in acceptorList:

		