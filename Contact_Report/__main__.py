
import sys
import glob
from hydrogen_bonds import *
#from bioinf import PDBATOMFileReader

for currentWildcard in sys.argv[1:]:

	for file_path in glob.glob(currentWildcard):

		reader = PDBATOMFileReader(file_path)
		#atomDict = {} # dict or ordered dict
		donorDict = {}
		acceptorDict = {}
		for atom in reader:
			#atomDict[atom.serial] = atom 
			if atom.is_donor:
				if atom.valence == 'sp2':
					donorDict[atom.serial] = Sp2DonorIQ(atom)
				else:
					donorDict[atom.serial] = Sp3DonorIQ(atom)
			if atom.is_acceptor:
				if atom.valence = 'sp2':
					acceptorDict[atom.serial] = Sp2AcceptorIQ(atom)
				else:
					acceptorDict[atom.serial] = Sp2AcceptorIQ(atom)
		print len(atomsList)
for donor in donorList:

	for acceptor in acceptorList:
		print 'bbq for you'

		