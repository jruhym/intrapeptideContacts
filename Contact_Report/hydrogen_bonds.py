import bioinf
from numpy import array, dot, arccos, rad2deg, ndarray, cross
from numpy.linalg import norm
from .constants import *
from collections import OrderedDict, namedtuple



class PDBATOMFileReader(object):
    def __init__(self, file_or_path):
        self._parse_atom_lines_filling_atoms_residues_and_chains_dicts(file_or_path)
        self._set_residues_and_chains_of_each_atom()
        self._set_chain_of_each_residue_and_add_it_to_itsown_chain()

    def _parse_atom_lines_filling_atoms_residues_and_chains_dicts(self, file_or_path):
        if isinstance(file_or_path, basestring):
            f = open(file_or_path, 'r')
        else:  
            f = file_or_path
        self._atoms = OrderedDict()
        self._residues = OrderedDict()
        self._chains = OrderedDict()
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith('ATOM'):
                atom = AtomIQ(clean_line)
                self._atoms[atom.serial] = atom
                try:
                    self._chains[atom.chainID].add_atom(atom)
                except KeyError:
                    self._chains[atom.chainID] = ChainIQ(atom)
                try: 
                    self._residues[atom.chainID + atom.uid].add_atom(atom)
                except KeyError:
                    self._residues[atom.chainID + atom.uid] = ResidueIQ(atom)
        f.close()

    def _set_residues_and_chains_of_each_atom(self):
        for atom in self._atoms.itervalues():
            atom.set_Residue(self._residues[atom.chainID + atom.uid])
            atom.set_Chain(self._chains[atom.chainID])

    def _set_chain_of_each_residue_and_add_it_to_itsown_chain(self):
        for residue in self._residues.itervalues():
            residue.set_Chain(self._chains[residue.chainID])
            self._chains[residue.chainID].add_residue(residue)

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
        self._residue = None
        self._chain = None
        self._chainID = pdb_atom_line.chainID
        self._coordinates = array([
            float(pdb_atom_line.x),
            float(pdb_atom_line.y),
            float(pdb_atom_line.z)
            ])
        self._participant = \
            HBondParticipant.generate_participant_by_valence(self)

    def set_Residue(self, residue):
        assert isinstance(residue, ResidueIQ)
        assert residue.uid == self._resSeq
        # if None
        self._residue = residue  

    def set_Chain(self, chain):
        assert isinstance(chain, ChainIQ)
        assert chain.chainID == self._chainID
        if self._chain is None:
            self._chain = chain
        else:
            raise TypeError('chain was already set and thus was not None')

    res_name = property(lambda self: self._res_name)
    uid = property(lambda self: self._resSeq)
    name = property(lambda self: self._name)
    chainID = property(lambda self: self._chainID)
    coordinates = property(lambda self: self._coordinates)
    serial = property(lambda self: self._serial)
    residue = property(lambda self: self._residue, set_Residue)
    chain = property(lambda self: self._chain, set_Chain)
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
        self._acceptor_list = []
        self._donor_list = []

    @staticmethod
    def _am_I_when_given(atom, currentGroup, bb_atom_name):
        assert isinstance(atom, AtomIQ)
        assert isinstance(currentGroup, HBondGroup)
        assert bb_atom_name in ('N', 'O')
        return (
            atom.name in currentGroup.atoms_str_tupl and atom.res_name == currentGroup.residue.upper()
            ) or (
            atom.name == bb_atom_name and currentGroup.residue == 'Peptide'
            )

    @staticmethod
    def generate_participant_by_valence(atom):
        assert isinstance(atom, AtomIQ)
        bb = namedtuple('backbone_Hbond_atom_name', ['donor','acceptor'])('N', 'O')
        is_acceptor = False
        is_donor = False
        H_bond_donor_radius = None
        max_num_H_donations = None
        H_bond_acceptor_radius = None
        max_num_H_acceptance = None

        for currentDonorGroup in list_of_hbond_donor_groups:
            if HBondParticipant._am_I_when_given(atom, currentDonorGroup, bb.donor):
                is_donor = True
                valence = currentDonorGroup.valence
                H_bond_donor_radius = currentDonorGroup.H_bond_radius
                max_num_H_donations = currentDonorGroup.max_num_H_bonds
                NN = currentDonorGroup.NN
                NNN = currentDonorGroup.NNN

        for currentAcceptorGroup in list_of_hbond_acceptor_groups:
            if HBondParticipant._am_I_when_given(atom, currentAcceptorGroup, bb.acceptor):
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
 
            elif valence == 'sp3':
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
    atom = property(lambda self: self._atom)
    NN = property(lambda self: self._NN)
    NNN = property(lambda self: self._NNN)



class AngleMinimum(namedtuple('AngleMinimum', ['as_donor', 'as_acceptor'])):
    def is_if(self, donor=True):
        return self.as_donor if donor else self.as_acceptor



class PlaneAngleMaximum(namedtuple('AngleMinimum', ['as_donor', 'as_acceptor'])):
    def is_if(self, donor=True):
        return self.as_donor if donor else self.as_acceptor



class Sp3HBondParticipant(HBondParticipant):
    _angle_min = AngleMinimum(90., 60.)

    def _distance_is_ok(self, partner):
        M = self._atom.coordinates
        P = partner.atom.coordinates
        distance = norm(M - P)
        if distance < self._H_bond_donor_radius + \
            partner.H_bond_acceptor_radius:
            return distance
        else:
            return False

    @staticmethod        
    def angle_is(ba, bc):
        assert isinstance(ba, ndarray)
        assert isinstance(bc, ndarray)
        return rad2deg(arccos(dot(bc, ba) / (norm(bc) * norm(ba))))

    def angle_is_ok(self, MtP, MtMM, as_donor=True):
        angle = self.angle_is(MtP, MtMM)
        return angle < 180. and angle > self._angle_min.is_if(as_donor)

    def planarity_is_ok(self, MtP, MtMM, MMtMMM, as_donor=True):
        return True

    @staticmethod
    def can_I_bond_to_partner(myself, partner, as_donor=True):
        assert isinstance(myself, HBondParticipant)
        assert isinstance(partner, HBondParticipant)        
        M = myself.atom.coordinates
        P = partner.atom.coordinates
        MM = myself.atom.residue.atoms[myself.NN].coordinates
        MtMM = MM - M
        MtP = P - M
        if myself.angle_is_ok(MtP, MtMM, as_donor):
            MMM = myself.atom.residue.atoms[myself.NNN].coordinates
            MMtMMM = MMM - MM
            if myself.planarity_is_ok(MtP, MtMM, MMtMMM, as_donor):
                return True

    def is_H_bond_mutual(self, partner):
        assert isinstance(partner, HBondParticipant)
        distance_or_is_ok = self._distance_is_ok(partner)
        if distance_or_is_ok and \
            self.can_I_bond_to_partner(self, partner) and \
            self.can_I_bond_to_partner(partner, self, as_donor=False):
            partner.append_donor_list(self)
            self.append_acceptor_list(partner)
            return distance_or_is_ok
    
    def append_donor_list(self, potential_h_donor):
        self.donor_list.append(potential_h_donor)

    def append_acceptor_list(self, potential_h_acceptor):
        self.acceptor_list.append(potential_h_acceptor)

    valence = property(lambda valence:'sp3')
    acceptor_list = property(lambda self: self._acceptor_list, append_acceptor_list)
    donor_list = property(lambda self: self._donor_list, append_donor_list)



class Sp2HBondParticipant(Sp3HBondParticipant):
    _angle_min = AngleMinimum(90., 90.)
    _plane_angle_max = PlaneAngleMaximum(60., 90.)

    @staticmethod
    def planarity_is(ba, bc, cd):
        assert isinstance(ba, ndarray)
        assert isinstance(bc, ndarray)
        assert isinstance(cd, ndarray)        
        my_plane_norm = cross(ba, bc)
        perndclr_bc_in_plane = cross(bc, my_plane_norm)
        torsion_angle_center = 0 if dot(cd, perndclr_bc_in_plane) > 0. else 180.
        plane_norm_w_partner = cross(-bc, cd)

        return abs(torsion_angle_center - Sp3HBondParticipant.angle_is(
            my_plane_norm, plane_norm_w_partner))


    def planarity_is_ok(self, MtP, MtMM, MMtMMM, as_donor=True):
        plane_angle = self.planarity_is(MMtMMM, -MtMM, MtP)
        return plane_angle < self._plane_angle_max.is_if(as_donor)

    valence = property(lambda valence: 'sp2')



class ResidueIQ(object):
    def __init__(self, atom):
        assert isinstance(atom, AtomIQ)
        self._atoms = {atom.name: atom}
        self._abbr = atom.res_name
        self._uid = atom.uid
        self._chainID = atom.chainID
        self._chain = None

    def add_atom(self, atom):
        assert isinstance(atom, AtomIQ)
        assert self.uid == atom.uid
        if atom.name not in self._atoms:
            self._atoms[atom.name] = atom

    def set_Chain(self, chain):
        assert isinstance(chain, ChainIQ)
        assert chain.chainID == self._chainID
        if self._chain is None:
            self._chain = chain
        else:
            raise TypeError('chain was already set and thus was not None')

    atoms = property(lambda self: self._atoms, add_atom)
    uid = property(lambda self: self._uid)
    chainID = property(lambda self: self._chainID)
    abbr = property(lambda self: self._abbr)
    chain = property(lambda self: self._chain, set_Chain)



class ChainIQ(object):
    def __init__(self, atom):
        assert isinstance(atom, AtomIQ)
        self._chainID = atom.chainID
        self._atoms = OrderedDict({atom.serial: atom})
        self._residues = OrderedDict({atom.uid: atom.residue})

    def add_atom(self, atom):
        assert isinstance(atom, AtomIQ)
        if atom.serial not in self._atoms:
            self._atoms[atom.serial] = atom
        else:
            raise KeyError('%s already exists in list of atoms for chain %s' %
                (atom.serial, self._chainID))

    def add_residue(self, residue):
        assert isinstance(residue, ResidueIQ)
        if residue.uid not in self._residues:
            self._residues[residue.uid] = residue
        # The following does not work because the first residue is added as 
        # chain is initialized.
        #else:
        #    raise KeyError('%s already exists in list of residues for chain %s' %
        #        (residue.uid, self._chainID))

    atoms = property(lambda self: self._atoms, add_atom)
    residues = property(lambda self: self._residues, add_residue)
    chainID = property(lambda self: self._chainID)





