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
		assert type(atom_str_tupl) == tuple
		assert len(atoms_str_tupl) > 0
		for atom_str in atoms_str_tupl:
			assert type(atom_str) == str
		self.__atoms_str_tupl__ = atoms_str_tupl

	def set_residue(self, residue):
		assert type(residue) == str
		assert len(residue) == 3 
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


#class atom_IQ():

#These objects comes from Stickle et al.
donor_Nsp2_I_peptide = HBondGroup()
donor_Nsp2_I_peptide.set_valance("sp2")
donor_Nsp2_I_peptide.set_residue("peptide")
donor_Nsp2_I_peptide.set_H_bondonor_radius(1.90)
donor_Nsp2_I_peptide.set_atoms_str_tupl(["NE1"])
donor_Nsp2_I_peptide.set_max_num_H_bonds(1)
donor_Nsp2_I_peptide.set_NN("CA")
donor_Nsp2_I_peptide.set_NNN("C")

donor_Nsp2_I_Trp = HBondGroup()
donor_Nsp2_I_Trp.set_valance("sp2")
donor_Nsp2_I_Trp.set_residue("Trp")
donor_Nsp2_I_Trp.set_H_bondonor_radius(1.90)
donor_Nsp2_I_Trp.set_atoms_str_tupl(["NE1"])
donor_Nsp2_I_Trp.set_max_num_H_bonds(1)



