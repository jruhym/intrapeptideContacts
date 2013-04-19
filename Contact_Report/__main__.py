
import sys
import glob
from .hydrogen_bonds import *
#from bioinf import PDBATOMFileReader

for currentWildcard in sys.argv[1:]:

	for file_path in glob.glob(currentWildcard):

		reader = PDBATOMFileReader(file_path)
		donorDict = {}
		acceptorDict = {}
		for atom in reader:
			if atom.is_donor:
				donorDict[atom.serial] = atom
			if atom.is_acceptor:
				acceptorDict[atom.serial] = atom
for donor in donorDict.itervalues():

	for acceptor in acceptorDict.intervalues():
		print 'bbq for you'

		