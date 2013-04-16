import bioinf
from numpy import array
from numpy.linalg import norm

class PDBATOMFileReader(object):#FileReader):
    def __init__(self, file_or_path):
        f = file_or_path if not isinstance(file_or_path, basestring) \
        else  open(file_or_path, 'r')
        self._contents = f.read()
        self._atoms = {}
        self._residues = {}
        for line in f:
        	clean_line = line.strip()
        	if clean_line.startswith('ATOM'):
        		atom = AtomIQ(clean_line)
        		self._atoms[atom.serial] = atom
        		try: 
        			self._residues[atom.uid].add_atom(atom)
        		except KeyError:
        			self._residues[atom.uid] = ResidueIQ(atom)
        for atom in self._atoms.itervalues():
        	atom.set_Residue(self._residues[atom.uid])




    def __iter__(self):
        for atom in self._atoms:
                yield atom


class HBondGroup(object):

	valence = property(lambda self: self._valence)
	atoms_str_tupl = property(lambda self: self._atoms_str_tupl)
	residue = property(lambda self: self._residue)
	H_bond_radius = property(lambda self: self._H_bond_radius)
	max_num_H_bonds = property(lambda self: self_.max_num_H_bonds)
	NN = property(lambda self: self_.NN)
	NNN = property(lambda self: self_.NNN)

	def __init__(self, valence=None, atoms_str_tupl=None, residue=None, 
			H_bond_radius=None, max_num_H_bonds=None, NN=None, NNN=None):
	
		HBondGroup.assert_valence(valence)
		HBondGroup.assert_atoms_str_tupl(atoms_str_tupl)
		HBondGroup.assert_residue(residue)
		HBondGroup.assert_H_bond_radius(H_bond_radius)
		HBondGroup.assert_max_num_H_bonds(max_num_H_bonds)
		HBondGroup.assert_NN(NN)
		HBondGroup.assert_NNN(NNN)

		self._valence = valence
		self._atoms_str_tupl = atoms_str_tupl
		self._residue = residue
		self._H_bond_radius = H_bond_radius
		self._max_num_H_bonds = max_num_H_bonds
		self._NN = NN
		self._NNN = NNN

	@staticmethod
	def assert_valence(valence):
		assert isinstance(valence, basestring)

	@staticmethod
	def assert_atoms_str_tupl(atoms_str_tupl):
		assert type(atoms_str_tupl) == tuple
		assert len(atoms_str_tupl) > 0
		for atom_str in atoms_str_tupl:
			assert type(atom_str) == str

	@staticmethod
	def assert_residue(residue):
		assert type(residue) == str
		assert len(residue) == 3 or residue == 'Peptide'

	@staticmethod
	def assert_H_bond_radius(H_bond_radius):
		assert type(H_bond_radius) == float or type(H_bond_radius) == int
		assert H_bond_radius > 0
		assert H_bond_radius <= 2.1

	@staticmethod
	def assert_max_num_H_bonds(max_num_H_bonds):
		assert type(max_num_H_bonds) == int
		assert max_num_H_bonds > 0
		assert max_num_H_bonds < 4

	@staticmethod
	def assert_NN(NN):
		assert type(NN) == str

	@staticmethod
	def assert_NNN(NNN):
		assert type(NNN) == str



class AtomIQ(object):

	def __init__(self, pdbAtomLine):
		assert isinstance(pdbAtomLine, basestring)
		atomPdbProperties = bioinf.PDBAtomLine.parse_string(pdbAtomLine)
		self._res_name = atomPdbProperties.resName
		self._resSeq = atomPdbProperties.resSeq
		self._name = atomPdbProperties.name
		self._serial = atomPdbProperties.serial
		self._valence = None
		self._H_bond_donor_radius = None
		self._H_bond_acceptor_radius = None
		self._residue = None
		self._coordinates = array([float(atomPdbProperties.x),
			float(atomPdbProperties.y),
			float(atomPdbProperties.z)
			])
		for currentDonorGroup in list_of_hbond_donor_groups:
			if self._name in currentDonorGroup.atoms_str_tupl and \
				self._residue in currentDonorGroup.residue:
				self._is_donor = True
				self._valence = currentDonorGroup.valence
				self._H_bond_donor_radius = currentDonorGroup.H_bond_radius
			else:
				self._is_donor = False
		for currentAcceptorGroup in list_of_hbond_acceptor_groups:
			if self._name in currentAcceptorGroup.atoms_str_tupl and \
				self._residue in currentAcceptorGroup.residue:
				self._is_acceptor = True
				self._valence = currentAcceptorGroup.valence
				self._H_bond_acceptor_radius = \
					currentDonorGroup.H_bond_radius
			else:
				self._is_acceptor = False

	def set_Residue(self, residue):
		assert isinstance(residue, ResidueIQ)
		assert residue.uid == self._reqSeq
		self._residue = residue

	res_name = property(lambda self: self._res_name)
	uid = property(lambda self: self._resSeq)
	name = property(lambda self: self._name)
	is_donor = property(lambda self: self._is_donor)
	is_acceptor = property(lambda self: self._is_acceptor)
	coordinates = property(lambda self: self._coordinates)
	serial = property(lambda self: self._serial)
	valence = property(lambda self: self._valence)
	H_bond_donor_radius = property(lambda self: self._H_bond_donor_radius)
	H_bond_acceptor_radius = property(
		lambda self: self._H_bond_acceptor_radius
		)
	residue = property(lambda self: self_.residue, set_Residue)


class ResidueIQ(object):
	def __init__(self, atom):
		assert isinstance(atom, AtomIQ)
		self._atoms = {}
		self._atoms[atom.resSeq] = atom
		self._resSeq = atom.resSeq
		self._abbr = atom.residue

	def add_atom(self, atom):
		assert isinstance(atom, AtomIQ)
		if atom.name not in self._atoms:
			atoms[atom.name] = atom

	atoms = property(lambda self: self._atoms)
	uid = property(lambda self: self._resSeq)
	abbr = property(lambda self: self._abbr)



class Sp3DonorIQ(AtomIQ):

	def __init__(self):
		pass

	def am_I_bonded_to_acceptor(self, acceptor):
		distance = numpy.linalg.norm(self._coordinates - acceptor.coordinates)
		
		if distance < self._H_bond_donor_radius + acceptor.H_bond_acceptor_radius:

			A = acceptor.coordinates
			D = self._coordinates
			DD = 
			if angleA_D_DD > 90. and angleA_D_DD 180.:
				if 


class Sp2DonorIQ(AtomIQ):
	def __init__(self):
		pass


class Sp2AcceptorIQ(AtomIQ):
	pass

class Sp3AcceptorIQ(AtomIQ):
	pass


class BondIQ(object):
	pass


class HBondIQ(BondIQ):
	pass

	#__init__(self, donor, acceptor):
#		__self__.donorValence = donor.getValence()
#		__self__.acceptorValence = acceptor.getValence()



