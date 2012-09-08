# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.           

"""Atom class, used in Structure objects."""

import numpy
import warnings
import copy

from Bio.PDB.Entity import DisorderedEntityWrapper
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Vector import Vector
from Bio.Data import IUPACData

# For sorting atom names
_atom_name_dict={}
_atom_name_dict["N"]=1
_atom_name_dict["CA"]=2
_atom_name_dict["C"]=3
_atom_name_dict["O"]=4

class Atom(object):
    __slots__ = ('level', 'parent', 'name', 'fullname', 'coord', 'bfactor', 'occupancy', 
    'altloc', 'full_id', 'id', 'disordered_flag', 'anisou_array', 'siguij_array', 'sigatm_array', 
    'serial_number', 'xtra', 'element', 'mass')
    def __init__(self, name, coord, bfactor, occupancy, altloc, fullname, serial_number,
                 element=None):
        """
        Atom object.

        The Atom object stores atom name (both with and without spaces), 
        coordinates, B factor, occupancy, alternative location specifier
        and (optionally) anisotropic B factor and standard deviations of 
        B factor and positions.
  
        @param name: atom name (eg. "CA"). Note that spaces are normally stripped.
        @type name: string

        @param coord: atomic coordinates (x,y,z)
        @type coord: Numeric array (Float0, size 3)

        @param bfactor: isotropic B factor
        @type bfactor: number 

        @param occupancy: occupancy (0.0-1.0)
        @type occupancy: number

        @param altloc: alternative location specifier for disordered atoms
        @type altloc: string

        @param fullname: full atom name, including spaces, e.g. " CA ". Normally
        these spaces are stripped from the atom name. 
        @type fullname: string

        @param element: atom element, e.g. "C" for Carbon, "HG" for mercury,
        @type fullname: uppercase string (or None if unknown)
     
        """
        self.level="A"
        # Reference to the residue 
        self.parent=None
        # the atomic data
        self.name=name      # eg. CA, spaces are removed from atom name
        self.fullname=fullname  # e.g. " CA ", spaces included
        self.coord=coord
        self.bfactor=bfactor
        self.occupancy=occupancy
        self.altloc=altloc
        self.full_id=None   # (structure id, model id, chain id, residue id, atom id)
        self.id=name        # id of atom is the atom name (e.g. "CA")
        self.disordered_flag=0
        self.anisou_array=None
        self.siguij_array=None
        self.sigatm_array=None
        self.serial_number=serial_number
        # Dictionary that keeps addictional properties
        self.xtra=None
        # Element assignment routines
        upper = str.upper
        self.element = self._assign_element(upper(element))
        self.mass = self._assign_atom_mass()
        
    def _assign_element(self, element):
        """Tries to guess element from atom name if not recognised."""
        if element not in IUPACData.atom_weights:
            # Inorganic elements have their name shifted left by one position 
            #  (is a convention in PDB, but not part of the standard).
            # isdigit() check on last two characters to avoid mis-assignment of 
            # hydrogens atoms (GLN HE21 for example)

            if self.fullname[0] != " " and not self.fullname[2:].isdigit():
                putative_element = self.name.strip()
            else:
                # Hs may have digit in [0]
                if self.name[0].isdigit():
                    putative_element = self.name[1]
                else:
                    putative_element = self.name[0]
            
            if putative_element in IUPACData.atom_weights:
                # msg = "Used element %r for Atom (name=%s) with given element %r" \
                #       % (putative_element, self.name, element)  
                element = putative_element
            else:
                msg = "Could not assign element %r for Atom (name=%s) with given element %r" \
                      % (putative_element, self.name, element)
                element = " "
                warnings.warn(msg, PDBConstructionWarning)
                
        return element
        
    def _assign_atom_mass(self):
        # Needed for Bio/Struct/Geometry.py C.O.M. function
        return IUPACData.atom_weights.get(self.element, float('NaN'))


    # Special methods   

    def __repr__(self):
        "Print Atom object as <Atom atom_name>."
        return "<Atom %s>" % self.id

    def __sub__(self, other):
        """
        Calculate distance between two atoms.
        
        Example:
            >>> distance=atom1-atom2

        @param other: the other atom
        @type other: L{Atom}
        """
        diff=self.coord-other.coord
        return numpy.sqrt(numpy.dot(diff,diff))
    
    def __cmp__(self, other):
        """
        Compares two atoms based on their atom names.
        N,CA,C,O always come first. Inorganic always come last.
        """                  
        organic = set(['C', 'N', 'O', 'H'])                       
        name1=self.name
        name2=other.name
        # Same name, compare altlocs
        if name1==name2:
            return cmp(self.altloc, other.altloc)
        # Check if priority names
        index1 = _atom_name_dict.get(name1)
        index2 = _atom_name_dict.get(name2)
        if index1 and index2:
            return cmp(index1, index2)
        elif not index1 and index2:
            return 1
        elif not index2 and index1:
            return -1
        else: # If not priority
            # Inorganic come at the end
            if self.element in organic and other.element not in organic:
                return -1
            if self.element not in organic and other.element in organic:
                return 1
            # finally, alphabetatically compare names
            return cmp(name1, name2)        


    # Public methods    

    def flag_disorder(self):
        """Set the disordered flag to 1.

        The disordered flag indicates whether the atom is disordered or not.
        """
        self.disordered_flag=1

    def is_disordered(self):
        "Return the disordered flag (1 if disordered, 0 otherwise)."
        return self.disordered_flag 

    def get_full_id(self):
        """Return the full id of the atom.

        The full id of an atom is the tuple 
        (structure id, model id, chain id, residue id, atom name, altloc).
        """
        return self.parent.get_full_id()+((self.name, self.altloc),)

    def transform(self, rot, tran):
        """
        Apply rotation and translation to the atomic coordinates.

        Example:
                >>> rotation=rotmat(pi, Vector(1,0,0))
                >>> translation=array((0,0,1), 'f')
                >>> atom.transform(rotation, translation)

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        self.coord=numpy.dot(self.coord, rot)+tran
        
    def get_vector(self):
        """
        Return coordinates as Vector.

        @return: coordinates as 3D vector
        @rtype: Vector
        """
        x,y,z=self.coord
        return Vector(x,y,z)

    def copy(self):
        """
        Create a copy of the Atom.
        Parent information is lost.
        """
        # Do a shallow copy then explicitly copy what needs to be deeper.
        shallow = copy.copy(self)
        shallow.parent = None
        shallow.coord = copy.copy(self.coord)
        if isinstance(shallow.xtra, dict):
            shallow.xtra = self.xtra.copy()
        else:
            shallow.xtra = None
        return shallow


class DisorderedAtom(DisorderedEntityWrapper):
    """
    This class contains all Atom objects that represent the same disordered
    atom. One of these atoms is "selected" and all method calls not caught
    by DisorderedAtom are forwarded to the selected Atom object. In that way, a
    DisorderedAtom behaves exactly like a normal Atom. By default, the selected 
    Atom object represents the Atom object with the highest occupancy, but a 
    different Atom object can be selected by using the disordered_select(altloc) 
    method. 
    """
    def __init__(self, id):
        """
        Arguments:
        o id - string, atom name
        """
        self.last_occupancy = None # FIX for very rare negative occupancies.
        DisorderedEntityWrapper.__init__(self, id)

    # Special methods

    def __repr__(self):
        return "<Disordered Atom %s>" % self.id 

    def disordered_add(self, atom):
        "Add a disordered atom."
        # Add atom to dict, use altloc as key   
        atom.flag_disorder()
        # set the residue parent of the added atom
        residue=self.parent
        atom.parent = residue
        altloc=atom.altloc
        occupancy=atom.occupancy
        self[altloc]=atom
        if occupancy>self.last_occupancy:
            self.last_occupancy=occupancy
            self.disordered_select(altloc)
