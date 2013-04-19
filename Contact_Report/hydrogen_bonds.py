import bioinf
from numpy import array, dot, arccos, rad2deg
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
        f.next()
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
        #self._valence = None
        #self._H_bond_donor_radius = None
        #self._H_bond_acceptor_radius = None
        self._residue = None
        self._coordinates = array([float(pdb_atom_line.x),
            float(pdb_atom_line.y),
            float(pdb_atom_line.z)
            ])
        for currentDonorGroup in list_of_hbond_donor_groups:
            if self._name in currentDonorGroup.atoms_str_tupl and \
                self._res_name == currentDonorGroup.residue.upper():
                self._is_donor = True
                self._valence = currentDonorGroup.valence
                self._H_bond_donor_radius = currentDonorGroup.H_bond_radius
        for currentAcceptorGroup in list_of_hbond_acceptor_groups:
            if pdb_atom_line.name in currentAcceptorGroup.atoms_str_tupl and \
                self._res_name == currentAcceptorGroup.residue.upper():
                self._is_acceptor = True
                self._valence = currentAcceptorGroup.valence
                self._H_bond_acceptor_radius = currentDonorGroup.H_bond_radius
    
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
    #valence = property(lambda self: self._valence)
    #H_bond_donor_radius = property(lambda self: self._H_bond_donor_radius)
    #H_bond_acceptor_radius = property(
        #lambda self: self._H_bond_acceptor_radius
        #)
    residue = property(lambda self: self._residue, set_Residue)


class Hbond_participant(object):

    def __init__(self, atom):
        pass

    @staticmethod
    def generate_participant_valence(self, atom):
        assert isinstance(atom, AtomIQ)
        if atom.valence == 'sp2':
            return Sp2_Hbond_participant(atom)
        elif atom.valence == 'sp3':
            return Sp3_Hbond_participant(atom)


class Sp3_Hbond_participant(Hbond_participant):
    def _distance_is_ok_to(self, M, P, partner):

        distance = norm(M - P)
        if distance < self._H_bond_radius + parner.H_bond_radius:
            return distance
        else:
            return False


    def _angle_is(self, ba, bc):

        return rad2deg(arccos(dot(bc, ba) / (norm(bc) * norm(ba))))

    def _angle_is_ok(self, MtP, MtMM):
        angle = _angle_is(MtP, MtMM)
        if angle < 180. and angle > self._angle_min:
            return True
        else:
            return False

    def _planarity_is_ok(self, P, M, MM, MMM):
        return True

    def am_I_bonded_to_partner(self, parner):
            
        M = self._coordinates
        P = parner.coordinates

        if _distance_is_ok(M,P, partner):

            MM = self._residue[self._NN].coordinates

            MtMM = MM - M
            MtP = P - M
            if _angle_is_ok(MtP, MtMM):

                MMM = self._residue[self._NNN]
                MMtMMM = MMM - MM

                if _planarity_is_ok(P, M, MM, MMM):

                    return True
        
    H_bond_radius = property(lambda self: self._H_bond_radius)




class Sp2_Hbond_participant(Sp3_Hbond_participant):
        
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


