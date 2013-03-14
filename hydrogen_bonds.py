import bioinf as bi

class HBondGroup():

	def __init__(self):
		self.__valence__ = 'Not yet set'
		self.__acceptor_atoms_str_list__ = ['Not yet set']
		self.__donor_atoms_str_list__ = ['Not yet set']
		self.__residue__ = 'Not yet set'
		self.__H_bond_radius__ = 'Not yet set'
		self.__max_num_H_bonds__ = 'Not yet set'

	def get_valance(self):
		return self.__valence__

	def get_atoms_str_list(self):
		return self.__donor_atoms_str_list__

	def get_residue(self):
		return self.__residue__

	def get_H_bond_radius(self):
		return self.__H_bond_radius__

	def get_max_num_H_bonds(self):
		return self.__max_num_H_bonds__
	
	def set_valance(self, valence):
		self.__valence__ = valence

	def set__atoms_str_list(self, _atoms_str_list):
		#bi.TestListOfStr.setUp(_atoms_str_list)
		self.___atoms_str_list__ = _atoms_str_list

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


#class atom_IQ():

