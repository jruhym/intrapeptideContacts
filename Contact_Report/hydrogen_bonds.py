import bioinf
from numpy import array, dot, arccos, rad2deg, ndarray, cross
from numpy.linalg import norm
from .constants import *
from collections import OrderedDict



class PDBATOMFileReader(object):#FileReader):
    def __init__(self, file_or_path):

        if isinstance(file_or_path, basestring):
            f = open(file_or_path, 'r')
            print file_or_path
        else:  
            f = file_or_path
        self._atoms = OrderedDict()
        self._residues = OrderedDict()
        #f.next()
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith('ATOM'):
                atom = AtomIQ(clean_line)
                self._atoms[atom.serial] = atom
                try: 
                    self._residues[atom.uid].add_atom(atom)
                except KeyError:
                    self._residues[atom.uid] = ResidueIQ(atom)
        f.close()
        for atom in self._atoms.itervalues():
            atom.set_Residue(self._residues[atom.uid])

    def __iter__(self):
        for atom in self._atoms:
            yield self._atoms[atom]



class AtomIQ(object):
    def __init__(self, pdb_line):
        assert isinstance(pdb_line, basestring) 
        pdb_atom_line = bioinf.PDBAtomLine.parse_string(pdb_line)   
        self._res_name = pdb_atom_line.resName
        self._resSeq = pdb_atom_line.resSeq
        self._name = pdb_atom_line.name
        self._serial = pdb_atom_line.serial
        self._is_donor = False
        self._is_acceptor = False
        self._residue = None
        self._coordinates = array([float(pdb_atom_line.x),
            float(pdb_atom_line.y),
            float(pdb_atom_line.z)
            ])
        self._participant = \
            HBondParticipant.generate_participant_by_valence(self)

    def set_Residue(self, residue):
        assert isinstance(residue, ResidueIQ)
        assert residue.uid == self._resSeq
        self._residue = residue     

    res_name = property(lambda self: self._res_name)
    uid = property(lambda self: self._resSeq)
    name = property(lambda self: self._name)
    is_donor = property(lambda self: self._is_donor)
    is_acceptor = property(lambda self: self._is_acceptor)
    coordinates = property(lambda self: self._coordinates)
    serial = property(lambda self: self._serial)
    residue = property(lambda self: self._residue, set_Residue)
    participant = property(lambda self: self._participant)



class HBondParticipant(object):
    def __init__(self, atom, 
        is_donor=False, H_bond_donor_radius=None, max_num_H_donations=None,
        is_acceptor=False, H_bond_acceptor_radius=None, 
        max_num_H_acceptance=None, NN=None, NNN=None):
        assert isinstance(atom, AtomIQ)
        self._atom = atom
        self._is_acceptor = is_acceptor
        self._is_donor = is_donor
        self._H_bond_acceptor_radius = H_bond_acceptor_radius
        self._H_bond_donor_radius = H_bond_donor_radius
        self._max_num_H_acceptance = max_num_H_acceptance
        self._max_num_H_donations = max_num_H_donations
        # tried to set NN as atoms but residue not yet set by this point so 
        # leave NN as string to index residue.atoms later.
        self._NN = NN
        self._NNN = NNN

    @staticmethod
    def generate_participant_by_valence(atom):
        assert isinstance(atom, AtomIQ)
        is_acceptor = False
        is_donor = False
        valence = None
        H_bond_donor_radius = None
        max_num_H_donations = None
        H_bond_acceptor_radius = None
        max_num_H_acceptance = None
        NN = None
        NNN = None
        for currentDonorGroup in list_of_hbond_donor_groups:
            if (
                atom.name in currentDonorGroup.atoms_str_tupl and \
                atom.res_name == currentDonorGroup.residue.upper()
                ) or (
                    atom.name =='N' and \
                    currentDonorGroup.residue == 'Peptide'
                ):
                is_donor = True
                valence = currentDonorGroup.valence
                H_bond_donor_radius = currentDonorGroup.H_bond_radius
                max_num_H_donations = currentDonorGroup.max_num_H_bonds
                NN = currentDonorGroup.NN
                NNN = currentDonorGroup.NNN

        for currentAcceptorGroup in list_of_hbond_acceptor_groups:
            if (
                atom.name in currentAcceptorGroup.atoms_str_tupl and \
                atom.res_name == currentAcceptorGroup.residue.upper()
                ) or (
                    atom.name =='O' and \
                    currentAcceptorGroup.residue == 'Peptide'
                ):
                is_acceptor = True
                valence = currentAcceptorGroup.valence
                H_bond_acceptor_radius = currentDonorGroup.H_bond_radius
                max_num_H_acceptance = currentDonorGroup.max_num_H_bonds
                NN = currentAcceptorGroup.NN
                NNN = currentAcceptorGroup.NNN
        if is_acceptor or is_donor:
            if valence == 'sp2':
                return Sp2HBondParticipant(atom, 
                    is_donor, H_bond_donor_radius, max_num_H_donations,
                    is_acceptor, H_bond_acceptor_radius, max_num_H_acceptance,
                    NN, NNN
                    )
            else valence == 'sp3':
                return Sp3HBondParticipant(atom,
                    is_donor, H_bond_donor_radius, max_num_H_donations,
                    is_acceptor, H_bond_acceptor_radius, max_num_H_acceptance,
                    NN, NNN
                    )
        else:
            return None

    is_acceptor = property(lambda self: self._is_acceptor)
    is_donor = property(lambda self: self._is_donor)
    H_bond_acceptor_radius = property(
        lambda self: self._H_bond_acceptor_radius)
    H_bond_donor_radius = property(lambda self: self._H_bond_donor_radius)
    max_num_H_acceptance = property(lambda self: self._max_num_H_acceptance)
    max_num_H_donations = property(lambda self: self._max_num_H_donations)
    NN = property(lambda self: self._NN)
    NNN = property(lambda self: self._NNN)



class Sp3HBondParticipant(HBondParticipant):
    def _distance_is_ok(self, M, P, partner):
        distance = norm(M - P)
        if distance < self._H_bond_donor_radius + \
            partner.participant.H_bond_acceptor_radius:
            return distance
        else:
            return False

    @staticmethod        
    def angle_is(ba, bc):
        assert isinstance(ba, ndarray)
        assert isinstance(bc, ndarray)
        return rad2deg(arccos(dot(bc, ba) / (norm(bc) * norm(ba))))

    def _angle_is_ok(self, MtP, MtMM):
        angle = self.angle_is(MtP, MtMM)
        if angle < 180. and angle > self._angle_min:
            return True
        else:
            return False

    def _planarity_is_ok(self, P, M, MM, MMM):
        return True

    def is_H_bond_mutual(self, partner):
<<<<<<< HEAD
        pass

    def can_I_bond_to_partner(self, partner, as_donor=True):
=======
>>>>>>> 3f41ce754fe19262f533dec0fe07e2a2e41ec73b
        M = self._atom.coordinates
        P = partner.coordinates
        distance_or_is_ok = self._distance_is_ok(M, P, partner)
        if distance_or_is_ok:
            return _distance_is_ok and\
                can_I_bond_to_partner(self, partner) and\
                can_I_bond_to_partner(partner, self)

    def can_I_bond_to_partner(self, partner):
        MM = self._atom.residue.atoms[self._NN].coordinates
        MtMM = MM - M
        MtP = P - M
        if self._angle_is_ok(MtP, MtMM):
            MMM = self._atom.residue.atoms[self._NNN].coordinates
            MMtMMM = MMM - MM
            if self._planarity_is_ok(P, M, MM, MMM):
                return True
    
    valence = property(lambda valence:'sp3')



class Sp2HBondParticipant(Sp3HBondParticipant):
    
    @staticmethod
    def planarity_is(ba, bc, cd):
        my_plane_norm = cross(ba, bc)
        perndclr_bc_in_plane = cross(bc, my_plane_norm)
        if dot(cd, perndclr_bc_in_plane) > 0.:
            torsion_angle_center = 0.
        else: 
            torsion_angle_center = 180.
        plane_norm_w_partner = cross(-bc, cd)

        torsion_angle = torsion_angle_center - Sp3HBondParticipant.angle_is(
            my_plane_norm, plane_norm_w_partner)
        print torsion_angle_center
        return torsion_angle


    def _planarity_is_ok(self, MtP, MtMM, MMtMMM):
        MMtM = -MtMM
        my_plane_norm = cross(MMtMMM, MMtM)
        perndclr_MMtM_in_plane = cross(MMtM, my_plane_norm)
        if dot(MtA, perndclr_MMtM_in_plane) > 0.:
            torsion_angle_center = 0.
        else: 
            torsion_angle_center = 180.
        plane_norm_w_partner = cross(MMtM, MtP)
        torsion_angle = _angle_is(my_plane_norm, plane_norm_w_partner)
        if torsion_angle < torsion_angle_center + self._torsion_range and\
            torsion_angle > torsion_angle_center - self._torsion_range:
            return True
        else:
            return False

    valence = property(lambda valence: 'sp2')



class ResidueIQ(object):
    def __init__(self, atom):
        assert isinstance(atom, AtomIQ)
        self._atoms = {
            atom.name: atom
        }
        self._abbr = atom.res_name
        self._uid = atom.uid

    def add_atom(self, atom):
        assert isinstance(atom, AtomIQ)
        assert self.uid == atom.uid
        if atom.name not in self._atoms:
            self._atoms[atom.name] = atom

    atoms = property(lambda self: self._atoms, add_atom)
    uid = property(lambda self: self._uid)
    abbr = property(lambda self: self._abbr)


