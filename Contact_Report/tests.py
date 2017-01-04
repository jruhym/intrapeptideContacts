import unittest
import sure
from .hydrogen_bonds import *
from .constants import *
from StringIO import StringIO
from collections import OrderedDict, namedtuple


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
ATOM     95  N   LEU A  42     -50.778  14.805  12.534  1.00 79.85           N  
ATOM     96  CA  LEU A  42     -50.675  16.219  12.908  1.00 79.83           C  
ATOM     97  C   LEU A  42     -49.256  16.752  12.687  1.00 78.29           C  
ATOM     98  O   LEU A  42     -49.093  17.917  12.355  1.00 80.15           O  
ATOM     99  CB  LEU A  42     -51.106  16.448  14.364  1.00 80.71           C  
ATOM    100  CG  LEU A  42     -52.586  16.242  14.714  1.00 81.76           C  
ATOM    101  CD1 LEU A  42     -52.774  16.400  16.213  1.00 83.86           C  
ATOM    102  CD2 LEU A  42     -53.508  17.191  13.956  1.00 82.80           C  
.
.
.
ATOM    224  N   LYS A  60     -45.271  36.012  -2.184  1.00 76.98           N  
ATOM    225  CA  LYS A  60     -46.336  37.010  -2.323  1.00 76.24           C  
ATOM    226  C   LYS A  60     -45.932  38.438  -1.965  1.00 77.94           C  
ATOM    227  O   LYS A  60     -46.282  39.371  -2.694  1.00 75.40           O  
ATOM    228  CB  LYS A  60     -47.532  36.605  -1.466  1.00 75.25           C  
ATOM    229  CG  LYS A  60     -48.667  37.578  -1.530  1.00 77.75           C  
ATOM    230  CD  LYS A  60     -49.959  36.950  -1.164  1.00 78.52           C  
ATOM    231  CE  LYS A  60     -51.028  37.983  -1.189  1.00 79.74           C  
ATOM    232  NZ  LYS A  60     -52.299  37.333  -1.081  1.00 83.19           N  
.
.
.
ATOM   3505  N   GLU A 338     -51.900  38.271   7.622  1.00 76.38           N  
ATOM   3506  CA  GLU A 338     -52.591  38.694   6.411  1.00 74.47           C  
ATOM   3507  C   GLU A 338     -53.085  37.448   5.670  1.00 73.86           C  
ATOM   3508  O   GLU A 338     -54.258  37.370   5.282  1.00 73.06           O  
ATOM   3509  CB  GLU A 338     -51.659  39.512   5.516  1.00 75.29           C  
ATOM   3510  CG  GLU A 338     -52.343  40.279   4.365  1.00 79.56           C  
ATOM   3511  CD  GLU A 338     -52.695  39.442   3.113  1.00 85.28           C  
ATOM   3512  OE1 GLU A 338     -52.500  38.207   3.087  1.00 91.20           O  
ATOM   3513  OE2 GLU A 338     -53.176  40.048   2.131  1.00 89.12           O  

"""


class TestPdbAtomFileReader(unittest.TestCase):


    def setUp(self):

        self._CAs = []
        pdbfile = StringIO(pdb_contents)
        reader = PDBATOMFileReader(pdbfile)
        self._i = 0
        self._atoms_Dict = OrderedDict()
        for atom in reader:
            self._i += 1
            self._atoms_Dict[atom.serial] = atom
            if atom.name == 'CA':
                self._CAs.append(atom)

    def test_reader_should_yield_correct_number_of_atoms(self):
        self._i.should.equal(59)
    
    def test_reader_should_yield_correct_number_of_residues(self):
        len(self._CAs).should.equal(5)
    
    def test_reader_should_yield_correct_number_of_residues(self):
        len(self._CAs).should.equal(8)

    def test_reader_should_construct_residues_and_fill_them_with_atoms(self):
        len(self._CAs[0].residue.atoms).should.equal(4)

    def test_reader_should_construct_chain_and_fill_it_with_residues(self):
        len(self._CAs[0].chain.residues).should.equal(8)

    def test_Ser41OG_should_donate_to_Gly37O(self):
        self._atoms_Dict['94'].participant.can_bond_to_partner(
            self._atoms_Dict['94'].participant,
            self._atoms_Dict['65'].participant).should.be.ok
    
    def test_donation_between_Ser41OG_and_Gly37O_should_be_mutual(self):
        self._atoms_Dict['94'].participant.is_H_bond_mutual(
            self._atoms_Dict['65'].participant).should.be.ok
    
    def test_Ser41OG_should_be_a_donor(self):
        self._atoms_Dict['94'].participant.is_donor.should.be.ok

    def test_Ser41OG_valence_should_be_sp3(self):
        (self._atoms_Dict['94'].participant.valence).should.equal('sp3')

    def test_Ser41OG_H_bond_donor_radius_should_be_correct(self):
        self._atoms_Dict['94'].participant.H_bond_donor_radius.should.equal(1.7)

    def test_Ser41OG_max_number_of_donations_should_be_correct(self):
        self._atoms_Dict['94'].participant.max_num_H_donations.should.equal(1)

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
        abs(Sp3HBondParticipant.angle(ba, bc) - 96.8).should.be.below(0.1)

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
        (Sp2HBondParticipant.planarity(cd, bc, ba) - 71.2).should.be.below(
            0.1)

    def test_donor_should_have_correct_acceptor_in_list(self):
        self._atoms_Dict['94'].participant.is_H_bond_mutual(self._atoms_Dict['65'].participant)
        self._atoms_Dict['94'].participant.acceptor_list[0].atom.serial.should.equal('65')

    def test_acceptor_should_have_correct_donor_in_list(self):
        self._atoms_Dict['94'].participant.is_H_bond_mutual(self._atoms_Dict['65'].participant)
        self._atoms_Dict['65'].participant.donor_list[0].atom.serial.should.equal('94')
        
    def test_donor_has_too_many_h_bonds(self):
        keys = self._atoms_Dict.keys()
        N = len(keys)
        for i, key_i in enumerate(keys):
            atom_i = self._atoms_Dict[key_i]
            if atom_i.participant and atom_i.participant.is_donor:
                for j in range(i, N - 1):
                    key_j = keys[j]
                    atom_j = self._atoms_Dict[key_j]
                    if atom_j.participant and atom_j.participant.is_acceptor:
                       atom_i.participant.is_H_bond_mutual(atom_j.participant)       
                atom_i.participant.has_excessive_acceptors().shouldnt.be.ok
            
    def test_acceptor_has_too_many_h_bonds(self):
        keys = self._atoms_Dict.keys()
        N = len(keys)
        for i, key_i in enumerate(keys):
            atom_i = self._atoms_Dict[key_i]
            if atom_i.participant and atom_i.participant.is_donor:
                for j in range(i, N - 1):
                    key_j = keys[j]
                    atom_j = self._atoms_Dict[key_j]
                    if atom_j.participant and atom_j.participant.is_acceptor:
                        atom_i.participant.is_H_bond_mutual(atom_j.participant)       
                if atom_i.participant.is_acceptor:
                    atom_i.participant.has_excessive_donors().shouldnt.be.ok


    def test_should_find_hyd_interaction_bn_ILE38CG2_and_LEU42CD1(self):
        self._atoms_Dict['72'].do_I_interact_HYDly_w(
            self._atoms_Dict['102']).should.be.ok

    def test_should_find_ion_bond_bn_LYS60NZ_and_GLU338OE1(self):
        self._atoms_Dict['232'].do_I_bond_IONly_w(
            self._atoms_Dict['3512']).should.be.ok



