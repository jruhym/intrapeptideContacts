import unittest
import sure
from .hydrogen_bonds import *
from .constants import *
from StringIO import StringIO


class TestAtomIq(unittest.TestCase):

	def setUp(self):
		self._atom = AtomIQ('ATOM     94  OG  SER A  41     -52.128  11.044  10.561  1.00 82.68           O')

	def test_atom_should_have_correct_serial(self):
		(self._atom.serial).should.equal('94')

	def test_atom_should_have_correct_name(self):
		self._atom.name.should.equal('OG')

	def test_atom_should_have_correct_res_name(self):
		self._atom.res_name.should.equal('SER')

	def test_atom_should_have_correct_uid(self):
		self._atom.uid.should.equal('41')

pdb_contents = """
ATOM     62  N   GLY A  37     -50.959   7.324  15.287  1.00 83.06           N  
ATOM     63  CA  GLY A  37     -51.227   7.713  13.911  1.00 82.58           C  
ATOM     64  C   GLY A  37     -51.710   9.144  13.811  1.00 81.74           C  
ATOM     65  O   GLY A  37     -51.171   9.912  13.024  1.00 83.08           O  
ATOM     66  N   ILE A  38     -52.724   9.494  14.611  1.00 81.34           N  
ATOM     67  CA  ILE A  38     -53.260  10.869  14.659  1.00 80.88           C  
ATOM     68  C   ILE A  38     -52.153  11.874  15.003  1.00 80.13           C  
ATOM     69  O   ILE A  38     -52.121  12.964  14.453  1.00 81.84           O  
ATOM     70  CB  ILE A  38     -54.432  11.032  15.686  1.00 81.07           C  
ATOM     71  CG1 ILE A  38     -55.635  10.132  15.356  1.00 81.29           C  
ATOM     72  CG2 ILE A  38     -54.912  12.479  15.745  1.00 80.45           C  
ATOM     73  CD1 ILE A  38     -56.310  10.413  14.026  1.00 85.05           C  
ATOM     74  N   VAL A  39     -51.255  11.496  15.911  1.00 79.45           N  
ATOM     75  CA  VAL A  39     -50.128  12.346  16.306  1.00 79.67           C  
ATOM     76  C   VAL A  39     -49.166  12.575  15.154  1.00 79.83           C  
ATOM     77  O   VAL A  39     -48.803  13.717  14.863  1.00 81.09           O  
ATOM     78  CB  VAL A  39     -49.352  11.740  17.496  1.00 81.11           C  
ATOM     79  CG1 VAL A  39     -47.996  12.425  17.671  1.00 78.01           C  
ATOM     80  CG2 VAL A  39     -50.185  11.843  18.770  1.00 81.49           C  
ATOM     81  N   MET A  40     -48.754  11.492  14.506  1.00 80.50           N  
ATOM     82  CA  MET A  40     -47.881  11.589  13.320  1.00 80.67           C  
ATOM     83  C   MET A  40     -48.575  12.333  12.155  1.00 80.62           C  
ATOM     84  O   MET A  40     -47.915  13.033  11.390  1.00 78.80           O  
ATOM     85  CB  MET A  40     -47.405  10.197  12.870  1.00 80.79           C  
ATOM     86  CG  MET A  40     -46.294   9.613  13.745  1.00 80.56           C  
ATOM     87  SD  MET A  40     -45.995   7.863  13.419  1.00 83.92           S  
ATOM     88  CE  MET A  40     -44.486   7.588  14.344  1.00 83.67           C  
ATOM     89  N   SER A  41     -49.896  12.177  12.034  1.00 79.41           N  
ATOM     90  CA  SER A  41     -50.672  12.898  11.020  1.00 79.83           C  
ATOM     91  C   SER A  41     -50.651  14.403  11.269  1.00 80.21           C  
ATOM     92  O   SER A  41     -50.521  15.181  10.324  1.00 82.42           O  
ATOM     93  CB  SER A  41     -52.104  12.394  10.976  1.00 79.47           C  
ATOM     94  OG  SER A  41     -52.128  11.044  10.561  1.00 82.68           O  

"""


class TestPdbAtomFileReader(unittest.TestCase):


	def setUp(self):

		self._CAs = []
		pdbfile = StringIO(pdb_contents)
		reader = PDBATOMFileReader(pdbfile)
		self._i = 0
		self._atoms_Dict = {}
		for atom in reader:
			self._i += 1
			self._atoms_Dict[atom.serial] = atom
			if atom.name == 'CA':
				self._CAs.append(atom)

	def test_reader_should_yield_correct_number_of_atoms(self):
		self._i.should.equal(33)
	
	def test_reader_should_yield_correct_number_of_residues(self):
		len(self._CAs).should.equal(5)

	def test_reader_should_construct_residues_and_fill_them_with_atoms(self):
		len(self._CAs[0].residue.atoms).should.equal(4)

	def test_Ser41OG_should_donate_to_Gly37O(self):
		self._atoms_Dict['94'].participant.can_I_bond_to_partner(
			self._atoms_Dict['65']).should.be.ok
	
	def test_Ser41OG_should_be_a_donor(self):
		self._atoms_Dict['94'].participant.is_donor.should.be.ok

	def test_Ser41OG_valence_should_be_sp3(self):
		(self._atoms_Dict['94'].participant.valence).should.equal('sp3')

	def test_Ser41OG_H_bond_donor_radius_should_be_correct(self):
		self._atoms_Dict['94'].participant.H_bond_donor_radius.should.equal(
			1.7
			)

	def test_Ser41OG_max_number_of_donations_should_be_correct(self):
		self._atoms_Dict['94'].participant.max_num_H_donations.should.equal(
			1
			)

	def test_Ser41OG_NN_should_be_CB(self):
		self._atoms_Dict['94'].participant.NN.should.equal('CB')

	def test_Ser41OG_NNN_should_be_CA(self):
		self._atoms_Dict['94'].participant.NNN.should.equal('CA')

	def test_angle_between_SerOG_its_NN_and_GlyO_should_be_correct(self):
		OG = self._atoms_Dict['94']
		b = OG.coordinates
		CB = OG.residue.atoms[OG.participant.NN]
		a = CB.coordinates
		O = self._atoms_Dict['65']
		c = O.coordinates
		ba = a - b
		bc = c - b
		abs(Sp3HBondParticipant.angle_is(ba, bc) - 96.8).should.be.below(0.1)

	def test_torsion_angle_Gly37Ca_C_O_and_Ser41O_should_be_correct(self):
		OG = self._atoms_Dict['94']
		O = self._atoms_Dict['65']
		C = O.residue.atoms[O.participant.NN]
		CA = O.residue.atoms[O.participant.NNN]
		a = CA.coordinates
		b = C.coordinates
		c = O.coordinates
		d = OG.coordinates
		ba = a - b
		bc = c - b
		cd = d - c
		(Sp2HBondParticipant.planarity_is(cd, bc, ba) - 71.2).should.be.below(
			0.1)
