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
		assert len(residue) == 3 or residue == 'peptide'
		self.__residue__ = residue

	def set_H_bond_radius(self, H_bond_radius):
		assert type(H_bond_radius) == float or type(H_bond_radius) == int
		assert H_bond_radius > 0
		assert H_bond_radius <= 2.1
		self.__H_bond_radius__ = H_bond_radius

	def set_max_num_H_bonds(self, max_num_H_bonds):
		assert type(max_num_H_bonds) == int
		assert max_num_H_bonds > 0
		assert max_num_H_bonds < 3
		self.__max_num_H_bonds__ = max_num_H_bonds

	def set_NN(self, NN):
		assert type(NN) == str
		self.__NN__ = NN

	def set_NNN(self, NNN):
		assert type(NNN) == str
		self.__NNN__ = NNN


# These objects come from Table 1 from Stickle et al. 
# DOI: 10.1016/0022-2836(92)91058-W
donor_Nsp2_I_peptide = HBondGroup()
donor_Nsp2_I_peptide.set_valance("sp2")
donor_Nsp2_I_peptide.set_residue("peptide")
donor_Nsp2_I_peptide.set_H_bond_radius(1.90)
donor_Nsp2_I_peptide.set_atoms_str_tupl(("N",))
donor_Nsp2_I_peptide.set_max_num_H_bonds(1)
donor_Nsp2_I_peptide.set_NN("CA")
donor_Nsp2_I_peptide.set_NNN("C")

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
donor_Nsp2_II_Asn.set_max_num_H_bonds(1)
donor_Nsp2_II_Asn.set_NN("CG")
donor_Nsp2_II_Asn.set_NNN("CB")

donor_Nsp2_II_Gln = HBondGroup()
donor_Nsp2_II_Gln.set_valance("sp2")
donor_Nsp2_II_Gln.set_residue("Gln")
donor_Nsp2_II_Gln.set_H_bond_radius(1.90)
donor_Nsp2_II_Gln.set_atoms_str_tupl(("NE2",))
donor_Nsp2_II_Gln.set_max_num_H_bonds(1)
donor_Nsp2_II_Gln.set_NN("CD")
donor_Nsp2_II_Gln.set_NNN("CG")



class AtomIQ():

	def __init__(self):
		# We may be able to set these at init automatically. 
		self.__residue__ = 'Not Yet Set'
		self.__name__ = 'Not Yet Set'
		#self.__is_donor__ = -1
		#self.__is_acceptor = -1
		self.__valence__ = 'Not Yet Set'
		self.__coordinates__ = [0.,0.,0.]

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




