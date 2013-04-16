import bioinf
from numpy import array
from numpy.linalg import norm

class PDBATOMFileReader(object):#FileReader):
    def __init__(self, file_or_path):
        f = file_or_path if not isinstance(file_or_path, basestring) \
        else  open(file_or_path, 'r')
        self._contents = f.read()

    def __iter__(self):
        for line in self._contents.splitlines():
            clean_line = line.strip()
            if clean_line.startswith('ATOM'):
                yield AtomIQ(clean_line)


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

# These objects come from Table 1 from Stickle et al. 
# DOI: 10.1016/0022-2836(92)91058-W
donor_Nsp2_I_Peptide = HBondGroup(
	valence = 'sp2',
	residue = 'Peptide',
	H_bond_radius = 1.90,
	atoms_str_tupl = ('N',),
	max_num_H_bonds = 1,
	NN = 'CA',
	NNN = 'C'
)

donor_Nsp2_I_Trp = HBondGroup(
	valence = 'sp2',
	residue = 'Trp',
	H_bond_radius = 1.90,
	atoms_str_tupl = ('NE1',),
	max_num_H_bonds = 1,
	NN = 'CD1',
	NNN = 'CE2'
)

donor_Nsp2_II_Asn = HBondGroup(
	valence = 'sp2',
	residue = 'Asn',
	H_bond_radius = 1.90,
	atoms_str_tupl = ("ND2",),
	max_num_H_bonds = 2,
	NN = 'CG',
	NNN = 'CB'
)

donor_Nsp2_II_Gln = HBondGroup(
	valence = "sp2",
	residue = "Gln",
	H_bond_radius = 1.90,
	atoms_str_tupl = ("NE2",),
	max_num_H_bonds = 2,
	NN = "CD",
	NNN = "CG"
)

donor_Nsp2_III_Arg = HBondGroup(
	valence = "sp2",
	residue = "Arg",
	H_bond_radius = 1.90,
	atoms_str_tupl = ("NH1", "NH2"),
	max_num_H_bonds = 2,	
	NN = "CZ",
	NNN = "NE"
)

donor_Nsp2_IV_Arg = HBondGroup(
	valence = "sp2",
	residue = "Arg",
	H_bond_radius = 1.90,
	atoms_str_tupl = ("NE",),
	max_num_H_bonds = 1,
	NN = "CZ",
	NNN = "CD"
)

donor_Nsp2_V_His = HBondGroup(
	valence = "sp2",
	residue = "His",
	H_bond_radius = 1.90,
	atoms_str_tupl = ("NE2",),
	max_num_H_bonds = 1,
	NN = "CE1",
	NNN = "CD2"
)

donor_Nsp3_I_Lys = HBondGroup(
	valence = "sp3",
	residue = "Lys",
	H_bond_radius = 2.10,
	atoms_str_tupl = ("NZ",),
	max_num_H_bonds = 3,
	NN = "CE",
	NNN = "CD"
)

donor_Osp3_I_Ser = HBondGroup(
	valence = "sp3",
	residue = "Ser",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OG",),
	max_num_H_bonds = 1,
	NN = "CB",
	NNN = "CA"
)

donor_Osp3_I_Thr = HBondGroup(
	valence = "sp3",
	residue = "Thr",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OG1",),
	max_num_H_bonds = 1,
	NN = "CB",
	NNN = "CA"
)

donor_Osp2_I_Tyr = HBondGroup(
	valence = "sp2",
	residue = "Tyr",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OH",),
	max_num_H_bonds = 1,
	NN = "CZ",
	NNN = "CE1"
)

acceptor_Nsp2_I_His_ND1 = HBondGroup(
	valence = "sp2",
	residue = "His",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("ND1",),
	max_num_H_bonds = 1,
	NN = "CE1",
	NNN = "CG"
)

acceptor_Nsp2_I_His_NE2 = HBondGroup(
	valence = "sp2",
	residue = "His",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("NE2",),
	max_num_H_bonds = 1,
	NN = "CE1",
	NNN = "CD2"
)

acceptor_Osp3_I_Ser = HBondGroup(
	valence = "sp3",
	residue = "Ser",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OG",),
	max_num_H_bonds = 2,
	NN = "CB",
	NNN = "CA"
)

acceptor_Osp3_I_Thr = HBondGroup(
	valence = "sp3",
	residue = "Thr",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OG1",),
	max_num_H_bonds = 2,
	NN = "CB",
	NNN = "CA"
)

acceptor_Osp2_I_Peptide = HBondGroup(
	valence = "sp2",
	residue = "Peptide",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("O",),
	max_num_H_bonds = 2,
	NN = "C",
	NNN = "CA"
)

acceptor_Osp2_I_Asn = HBondGroup(
	valence = "sp2",
	residue = "Asn",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("OD1",),
	max_num_H_bonds = 2,
	NN = "CG",
	NNN = "CB"
)

acceptor_Osp2_I_Gln = HBondGroup(
	valence = "sp2",
	residue = "Gln",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("OD1",),
	max_num_H_bonds = 2,
	NN = "CD",
	NNN = "CG"
)

acceptor_Osp2_II_Asp = HBondGroup(
	valence = "sp2",
	residue = "Asp",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("OD1", "OD2"),
	max_num_H_bonds = 2,
	NN = "CG",
	NNN = "CB"
)

acceptor_Osp2_II_Glu = HBondGroup(
	valence = "sp2",
	residue = "Glu",
	H_bond_radius = 1.60,
	atoms_str_tupl = ("OE1", "OE2"),
	max_num_H_bonds = 2,
	NN = "CD",
	NNN = "CG"
)

acceptor_Osp2_III_Tyr = HBondGroup(
	valence = "sp2",
	residue = "Tyr",
	H_bond_radius = 1.70,
	atoms_str_tupl = ("OH",),
	max_num_H_bonds = 1,
	NN = "CZ",
	NNN = "CE1"
)

acceptor_Ssp3_I_Met = HBondGroup(
	valence = "sp3",
	residue = "Met",
	H_bond_radius = 1.95,
	atoms_str_tupl = ("SD",),
	max_num_H_bonds = 2,
	NN = "CE",
	NNN = "CG"
)

acceptor_Ssp3_II_Cys = HBondGroup(
	valence = "sp3",
	residue = "Cys",
	H_bond_radius = 2.10,
	atoms_str_tupl = ("SG",),
	max_num_H_bonds = 2,
	NN = "CB",
	NNN = "CA"
)

list_of_hbond_donor_groups = [
	donor_Nsp2_I_Peptide,
	donor_Nsp2_I_Trp,
	donor_Nsp2_II_Asn,
	donor_Nsp2_II_Gln,
	donor_Nsp2_III_Arg,
	donor_Nsp2_IV_Arg,
	donor_Nsp2_V_His,
	donor_Nsp3_I_Lys,
	donor_Osp3_I_Ser,
	donor_Osp3_I_Thr,
	donor_Osp2_I_Tyr
]

list_of_hbond_acceptor_groups = [
	acceptor_Nsp2_I_His_ND1,
	acceptor_Nsp2_I_His_NE2,
	acceptor_Osp3_I_Ser,
	acceptor_Osp3_I_Thr,
	acceptor_Osp2_I_Peptide,
	acceptor_Osp2_I_Asn,
	acceptor_Osp2_I_Gln,
	acceptor_Osp2_II_Asp,
	acceptor_Osp2_II_Glu,
	acceptor_Osp2_III_Tyr,
	acceptor_Ssp3_I_Met,
	acceptor_Ssp3_II_Cys
]



class AtomIQ(object):

	def __init__(self, pdbAtomLine):
		assert isinstance(pdbAtomLine, basestring)
		atomPdbProperties = bioinf.PDBAtomLine.parse_string(pdbAtomLine)
		self._residue = atomPdbProperties.resName
		self._name = atomPdbProperties.name
		self._serial = atomPdbProperties.serial
		self._valence = None
		self._H_bond_donor_radius = None
		self._H_bond_acceptor_radius = None
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

	residue = property(lambda self: self._residue)
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

	def am_I_bonded_to_acceptor(self, acceptor):
		distance = numpy.linalg.norm(self._coordinates - acceptor.coordinates)
		
		if distance < self._H_bond_donor_radius + acceptor.H_bond_acceptor_radius:
			



class DonorIQ(AtomIQ):
	pass


class AcceptorIQ(AtomIQ):
	pass


class BondIQ(object):
	pass


class HBondIQ(BondIQ):
	pass

	#__init__(self, donor, acceptor):
#		__self__.donorValence = donor.getValence()
#		__self__.acceptorValence = acceptor.getValence()



