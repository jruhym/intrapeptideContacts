class HBondGroup():

	def __init__(self):
		self.__valence__ = 'Not yet set'
		self.__acceptor_atoms_str_tupl__ = ['Not yet set']
		self.__donor_atoms_str_tupl__ = ['Not yet set']
		self.__residue__ = 'Not yet set'
		self.__H_bond_radius__ = 'Not yet set'
		self.__max_num_H_bonds__ = 'Not yet set'
		self.__NN__ = 'Not yet set'
		self.__NNN__ = 'Not yet set'

	def get_valance(self):
		return self.__valence__

	def get_atoms_str_tupl(self):
		return self.__donor_atoms_str_tupl__

	def get_residue(self):
		return self.__residue__

	def get_H_bond_radius(self):
		return self.__H_bond_radius__

	def get_max_num_H_bonds(self):
		return self.__max_num_H_bonds__

	def get_NN(self):
		return self.__NN__

	def get_NNN(self):
		return self.__NNN__
	
	def set_valance(self, valence):
		self.__valence__ = valence

	def set_atoms_str_tupl(self, atoms_str_tupl):
		assert type(atoms_str_tupl) == tuple
		assert len(atoms_str_tupl) > 0
		for atom_str in atoms_str_tupl:
			assert type(atom_str) == str
		self.__atoms_str_tupl__ = atoms_str_tupl

	def set_residue(self, residue):
		assert type(residue) == str
		assert len(residue) == 3 or residue == 'Peptide'
		self.__residue__ = residue

	def set_H_bond_radius(self, H_bond_radius):
		assert type(H_bond_radius) == float or type(H_bond_radius) == int
		assert H_bond_radius > 0
		assert H_bond_radius <= 2.1
		self.__H_bond_radius__ = H_bond_radius

	def set_max_num_H_bonds(self, max_num_H_bonds):
		assert type(max_num_H_bonds) == int
		assert max_num_H_bonds > 0
		assert max_num_H_bonds < 4
		self.__max_num_H_bonds__ = max_num_H_bonds

	def set_NN(self, NN):
		assert type(NN) == str
		self.__NN__ = NN

	def set_NNN(self, NNN):
		assert type(NNN) == str
		self.__NNN__ = NNN


# These objects come from Table 1 from Stickle et al. 
# DOI: 10.1016/0022-2836(92)91058-W
donor_Nsp2_I_Peptide = HBondGroup()
donor_Nsp2_I_Peptide.set_valance("sp2")
donor_Nsp2_I_Peptide.set_residue("Peptide")
donor_Nsp2_I_Peptide.set_H_bond_radius(1.90)
donor_Nsp2_I_Peptide.set_atoms_str_tupl(("N",))
donor_Nsp2_I_Peptide.set_max_num_H_bonds(1)
donor_Nsp2_I_Peptide.set_NN("CA")
donor_Nsp2_I_Peptide.set_NNN("C")

donor_Nsp2_I_Trp = HBondGroup()
donor_Nsp2_I_Trp.set_valance("sp2")
donor_Nsp2_I_Trp.set_residue("Trp")
donor_Nsp2_I_Trp.set_H_bond_radius(1.90)
donor_Nsp2_I_Trp.set_atoms_str_tupl(("NE1",))
donor_Nsp2_I_Trp.set_max_num_H_bonds(1)
donor_Nsp2_I_Trp.set_NN("CD1")
donor_Nsp2_I_Trp.set_NNN("CE2")

donor_Nsp2_II_Asn = HBondGroup()
donor_Nsp2_II_Asn.set_valance("sp2")
donor_Nsp2_II_Asn.set_residue("Asn")
donor_Nsp2_II_Asn.set_H_bond_radius(1.90)
donor_Nsp2_II_Asn.set_atoms_str_tupl(("ND2",))
donor_Nsp2_II_Asn.set_max_num_H_bonds(2)
donor_Nsp2_II_Asn.set_NN("CG")
donor_Nsp2_II_Asn.set_NNN("CB")

donor_Nsp2_II_Gln = HBondGroup()
donor_Nsp2_II_Gln.set_valance("sp2")
donor_Nsp2_II_Gln.set_residue("Gln")
donor_Nsp2_II_Gln.set_H_bond_radius(1.90)
donor_Nsp2_II_Gln.set_atoms_str_tupl(("NE2",))
donor_Nsp2_II_Gln.set_max_num_H_bonds(2)
donor_Nsp2_II_Gln.set_NN("CD")
donor_Nsp2_II_Gln.set_NNN("CG")

donor_Nsp2_III_Arg = HBondGroup()
donor_Nsp2_III_Arg.set_valance("sp2")
donor_Nsp2_III_Arg.set_residue("Arg")
donor_Nsp2_III_Arg.set_H_bond_radius(1.90)
donor_Nsp2_III_Arg.set_atoms_str_tupl(("NH1", "NH2"))
donor_Nsp2_III_Arg.set_max_num_H_bonds(2)
donor_Nsp2_III_Arg.set_NN("CZ")
donor_Nsp2_III_Arg.set_NNN("NE")

donor_Nsp2_IV_Arg = HBondGroup()
donor_Nsp2_IV_Arg.set_valance("sp2")
donor_Nsp2_IV_Arg.set_residue("Arg")
donor_Nsp2_IV_Arg.set_H_bond_radius(1.90)
donor_Nsp2_IV_Arg.set_atoms_str_tupl(("NE",))
donor_Nsp2_IV_Arg.set_max_num_H_bonds(1)
donor_Nsp2_IV_Arg.set_NN("CZ")
donor_Nsp2_IV_Arg.set_NNN("CD")

donor_Nsp2_V_His = HBondGroup()
donor_Nsp2_V_His.set_valance("sp2")
donor_Nsp2_V_His.set_residue("His")
donor_Nsp2_V_His.set_H_bond_radius(1.90)
donor_Nsp2_V_His.set_atoms_str_tupl(("NE2",))
donor_Nsp2_V_His.set_max_num_H_bonds(1)
donor_Nsp2_V_His.set_NN("CE1")
donor_Nsp2_V_His.set_NNN("CD2")

donor_Nsp3_I_Lys = HBondGroup()
donor_Nsp3_I_Lys.set_valance("sp3")
donor_Nsp3_I_Lys.set_residue("Lys")
donor_Nsp3_I_Lys.set_H_bond_radius(2.10)
donor_Nsp3_I_Lys.set_atoms_str_tupl(("NZ",))
donor_Nsp3_I_Lys.set_max_num_H_bonds(3)
donor_Nsp3_I_Lys.set_NN("CE")
donor_Nsp3_I_Lys.set_NNN("CD")

donor_Osp3_I_Ser = HBondGroup()
donor_Osp3_I_Ser.set_valance("sp3")
donor_Osp3_I_Ser.set_residue("Ser")
donor_Osp3_I_Ser.set_H_bond_radius(1.70)
donor_Osp3_I_Ser.set_atoms_str_tupl(("OG",))
donor_Osp3_I_Ser.set_max_num_H_bonds(1)
donor_Osp3_I_Ser.set_NN("CB")
donor_Osp3_I_Ser.set_NNN("CA")

donor_Osp3_I_Thr = HBondGroup()
donor_Osp3_I_Thr.set_valance("sp3")
donor_Osp3_I_Thr.set_residue("Thr")
donor_Osp3_I_Thr.set_H_bond_radius(1.70)
donor_Osp3_I_Thr.set_atoms_str_tupl(("OG1",))
donor_Osp3_I_Thr.set_max_num_H_bonds(1)
donor_Osp3_I_Thr.set_NN("CB")
donor_Osp3_I_Thr.set_NNN("CA")

donor_Osp2_I_Tyr = HBondGroup()
donor_Osp2_I_Tyr.set_valance("sp2")
donor_Osp2_I_Tyr.set_residue("Tyr")
donor_Osp2_I_Tyr.set_H_bond_radius(1.70)
donor_Osp2_I_Tyr.set_atoms_str_tupl(("OH",))
donor_Osp2_I_Tyr.set_max_num_H_bonds(1)
donor_Osp2_I_Tyr.set_NN("CZ")
donor_Osp2_I_Tyr.set_NNN("CE1")

acceptor_Nsp2_I_His_ND1 = HBondGroup()
acceptor_Nsp2_I_His_ND1.set_valance("sp2")
acceptor_Nsp2_I_His_ND1.set_residue("His")
acceptor_Nsp2_I_His_ND1.set_H_bond_radius(1.60)
acceptor_Nsp2_I_His_ND1.set_atoms_str_tupl(("ND1",))
acceptor_Nsp2_I_His_ND1.set_max_num_H_bonds(1)
acceptor_Nsp2_I_His_ND1.set_NN("CE1")
acceptor_Nsp2_I_His_ND1.set_NNN("CG")

acceptor_Nsp2_I_His_NE2 = HBondGroup()
acceptor_Nsp2_I_His_NE2.set_valance("sp2")
acceptor_Nsp2_I_His_NE2.set_residue("His")
acceptor_Nsp2_I_His_NE2.set_H_bond_radius(1.60)
acceptor_Nsp2_I_His_NE2.set_atoms_str_tupl(("NE2",))
acceptor_Nsp2_I_His_NE2.set_max_num_H_bonds(1)
acceptor_Nsp2_I_His_NE2.set_NN("CE1")
acceptor_Nsp2_I_His_NE2.set_NNN("CD2")

acceptor_Osp3_I_Ser = HBondGroup()
acceptor_Osp3_I_Ser.set_valance("sp3")
acceptor_Osp3_I_Ser.set_residue("Ser")
acceptor_Osp3_I_Ser.set_H_bond_radius(1.70)
acceptor_Osp3_I_Ser.set_atoms_str_tupl(("OG",))
acceptor_Osp3_I_Ser.set_max_num_H_bonds(2)
acceptor_Osp3_I_Ser.set_NN("CB")
acceptor_Osp3_I_Ser.set_NNN("CA")

acceptor_Osp3_I_Thr = HBondGroup()
acceptor_Osp3_I_Thr.set_valance("sp3")
acceptor_Osp3_I_Thr.set_residue("Thr")
acceptor_Osp3_I_Thr.set_H_bond_radius(1.70)
acceptor_Osp3_I_Thr.set_atoms_str_tupl(("OG1",))
acceptor_Osp3_I_Thr.set_max_num_H_bonds(2)
acceptor_Osp3_I_Thr.set_NN("CB")
acceptor_Osp3_I_Thr.set_NNN("CA")

acceptor_Osp2_I_Peptide = HBondGroup()
acceptor_Osp2_I_Peptide.set_valance("sp2")
acceptor_Osp2_I_Peptide.set_residue("Peptide")
acceptor_Osp2_I_Peptide.set_H_bond_radius(1.60)
acceptor_Osp2_I_Peptide.set_atoms_str_tupl(("O",))
acceptor_Osp2_I_Peptide.set_max_num_H_bonds(2)
acceptor_Osp2_I_Peptide.set_NN("C")
acceptor_Osp2_I_Peptide.set_NNN("CA")

acceptor_Osp2_I_Asn = HBondGroup()
acceptor_Osp2_I_Asn.set_valance("sp2")
acceptor_Osp2_I_Asn.set_residue("Asn")
acceptor_Osp2_I_Asn.set_H_bond_radius(1.60)
acceptor_Osp2_I_Asn.set_atoms_str_tupl(("OD1",))
acceptor_Osp2_I_Asn.set_max_num_H_bonds(2)
acceptor_Osp2_I_Asn.set_NN("CG")
acceptor_Osp2_I_Asn.set_NNN("CB")

acceptor_Osp2_I_Gln = HBondGroup()
acceptor_Osp2_I_Gln.set_valance("sp2")
acceptor_Osp2_I_Gln.set_residue("Gln")
acceptor_Osp2_I_Gln.set_H_bond_radius(1.60)
acceptor_Osp2_I_Gln.set_atoms_str_tupl(("OD1",))
acceptor_Osp2_I_Gln.set_max_num_H_bonds(2)
acceptor_Osp2_I_Gln.set_NN("CD")
acceptor_Osp2_I_Gln.set_NNN("CG")

acceptor_Osp2_II_Asp = HBondGroup()
acceptor_Osp2_II_Asp.set_valance("sp2")
acceptor_Osp2_II_Asp.set_residue("Asp")
acceptor_Osp2_II_Asp.set_H_bond_radius(1.60)
acceptor_Osp2_II_Asp.set_atoms_str_tupl(("OD1", "OD2"))
acceptor_Osp2_II_Asp.set_max_num_H_bonds(2)
acceptor_Osp2_II_Asp.set_NN("CG")
acceptor_Osp2_II_Asp.set_NNN("CB")

acceptor_Osp2_II_Glu = HBondGroup()
acceptor_Osp2_II_Glu.set_valance("sp2")
acceptor_Osp2_II_Glu.set_residue("Glu")
acceptor_Osp2_II_Glu.set_H_bond_radius(1.60)
acceptor_Osp2_II_Glu.set_atoms_str_tupl(("OE1", "OE2"))
acceptor_Osp2_II_Glu.set_max_num_H_bonds(2)
acceptor_Osp2_II_Glu.set_NN("CD")
acceptor_Osp2_II_Glu.set_NNN("CG")

acceptor_Osp2_III_Tyr = HBondGroup()
acceptor_Osp2_III_Tyr.set_valance("sp2")
acceptor_Osp2_III_Tyr.set_residue("Tyr")
acceptor_Osp2_III_Tyr.set_H_bond_radius(1.70)
acceptor_Osp2_III_Tyr.set_atoms_str_tupl(("OH",))
acceptor_Osp2_III_Tyr.set_max_num_H_bonds(1)
acceptor_Osp2_III_Tyr.set_NN("CZ")
acceptor_Osp2_III_Tyr.set_NNN("CE1")

acceptor_Ssp3_I_Met = HBondGroup()
acceptor_Ssp3_I_Met.set_valance("sp3")
acceptor_Ssp3_I_Met.set_residue("Met")
acceptor_Ssp3_I_Met.set_H_bond_radius(1.95)
acceptor_Ssp3_I_Met.set_atoms_str_tupl(("SD",))
acceptor_Ssp3_I_Met.set_max_num_H_bonds(2)
acceptor_Ssp3_I_Met.set_NN("CE")
acceptor_Ssp3_I_Met.set_NNN("CG")

acceptor_Ssp3_II_Cys = HBondGroup()
acceptor_Ssp3_II_Cys.set_valance("sp3")
acceptor_Ssp3_II_Cys.set_residue("Cys")
acceptor_Ssp3_II_Cys.set_H_bond_radius(2.10)
acceptor_Ssp3_II_Cys.set_atoms_str_tupl(("SG",))
acceptor_Ssp3_II_Cys.set_max_num_H_bonds(2)
acceptor_Ssp3_II_Cys.set_NN("CB")
acceptor_Ssp3_II_Cys.set_NNN("CA")



class AtomIQ():

	def __init__(self, pdbAtomLine):
		assert type(pdbAtomLine) == str
		atomPdbProperties = bioinf.PDBAtomLine.parse_string(pdbAtomLine)
		# We may be able to set these at init automatically. 
		self.__residue__ = atomPdbProperties.resName
		self.__name__ = atomPdbProperties.name
		#self.__is_donor__ = -1
		#self.__is_acceptor = -1
		self.__valence__ = 'Not Yet Set'
		self.__coordinates__ = [float(atomPdbProperties.x),
								float(atomPdbProperties.y),
								float(atomPdbProperties.z)]

	def get_residue(self):
		return self.__residue__
	def get_name(self):
		return self.__name__
	def get_is_donor(self):
		return self.__is_donor__
	def get_is_acceptor(self):
		return self.__is_acceptor__
	def get_coordinates(self):
		return self.__coordinates__

	def set_residue(self, residue):
		assert type(residue) == str
		self.__residue__ = residue
	def set_name(self, name):
		assert type(name) == str
		self.__name__ = name
	def set_is_donor(self, is_donor):
		assert type(is_donor) == bool  
		self.__is_donor__ = is_donor
	def set_is_acceptor(self, is_acceptor):
		assert type(is_acceptor) == bool  
		self.__is_acceptor__ = is_acceptor
	def set_coordinates(self, coordinates):
		assert type(coordinates) == list
		for x_i in coordinates:
			assert type(x_i) == float
		self.__coordinates__ = coordinates

class DonorIQ(AtomIQ):
	pass


class AcceptorIQ(AtomIQ):
	pass


class BondIQ(object):
	pass


class HBondIQ(BondIQ):
	pass




