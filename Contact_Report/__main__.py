
import sys
import glob
from .hydrogen_bonds import *

for currentWildcard in sys.argv[1:]:
	for file_path in glob.glob(currentWildcard):
		reader = PDBATOMFileReader(file_path)
		donorDict = {}
		acceptorDict = {}
<<<<<<< HEAD
		for atom in reader:
=======
		atoms = {}
		for atom in reader:
			atoms[atom.serial] = atom
>>>>>>> cad8c4abe604b6fbd52df945a85b31948079b250
			if atom.participant != False:
				if atom.participant.is_donor:
					donorDict[atom.serial] = atom
				if atom.participant.is_acceptor:
					acceptorDict[atom.serial] = atom
for donor in donorDict.itervalues():
	for acceptor in acceptorDict.intervalues():
		


		