# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

"""Chain class, used in Structure objects."""

from Bio.PDB.Entity import Entity
import warnings

class Chain(Entity):
    def __init__(self, id):
        self.level="C"
        Entity.__init__(self, id)

    # Private methods

    def _translate_id(self, id):
        """
        A residue id is normally a tuple (hetero flag, sequence identifier, 
        insertion code). Since for most residues the hetero flag and the 
        insertion code are blank (i.e. " "), you can just use the sequence 
        identifier to index a residue in a chain. The _translate_id method
        translates the sequence identifier to the (" ", sequence identifier,
        " ") tuple. 

        Arguments:
        o id - int, residue resseq 
        """
        if isinstance(id, int):
            id=(' ', id, ' ')
        return id
            
    # Special methods   
    def __getitem__(self, id):
        """Return the residue with given id.

        The id of a residue is (hetero flag, sequence identifier, insertion code). 
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method. 

        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.__getitem__(self, id)

    def __contains__(self, id):
        """True if a residue with given id is present in this chain.

        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.__contains__(self, id)

    def __delitem__(self, id):
        """
        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.__delitem__(self, id)

    def __repr__(self):
        return "<Chain id=%s>" % self.id
    
    def __cmp__(self, other):
        id1=self.id
        id2=other.id
        # make sure blank chains come last (often waters)
        if id1==" " and not id2==" ":
            return 1
        elif id2==" " and not id1==" ":
            return -1
        return cmp(id1, id2)

    # Public methods
    def get_unpacked_list(self):
        """Return a list of undisordered residues.

        Some Residue objects hide several disordered residues
        (DisorderedResidue objects). This method unpacks them, 
        ie. it returns a list of simple Residue objects.
        """
        unpacked_list=[]
        for residue in self:
            if residue.is_disordered()==2:
                for dresidue in residue.disordered_get_list():
                    unpacked_list.append(dresidue)
            else:
                unpacked_list.append(residue)
        return unpacked_list
    
    def _renumber_single_residue(self, old_id, new_id):
        """
        To Allow custom renumbering of residues.
        """
        # Identify residue
        residue = self[old_id]

        # Change old references in both dict and list
        del self.child_dict[residue.id]
        residue.id = new_id
        self.child_dict[residue.id] = residue
    
    # Public
    def get_residues(self):
        for r in self:
            yield r

    def get_atoms(self):
        for r in self:
            for a in r:
                yield a
    
    def renumber_residues(self, res_init=1, het_init=0):
        """ Renumbers the residues in a chain. Keeps gaps in the chain if they exist.
            
            Uses SEQRES records to correctly filter HETATMs. In its absence,
            checks residues for alpha carbons (modified residues have it)
            and tries to infer which should be considered as ATOMs.
                        
            Options:
            
            - res_init [int] (default:1)
            First residue number
            - het_init [int] (default:0)
            First HETATM number. Default means first HETATM will be numbered 1000.

        """
        # Derive list of regular residues from SEQRES
        # Safe way to catch MSE and alike.
        h = self.parent.parent.header
        if 'SEQRES' in h and self.id in h['SEQRES']:
            seqres = h['SEQRES'][self.id]
            filter_by_ca = False
        else:
            warnings.warn("WARNING: SEQRES field could not be retrieved for chain %s\n"
                          "HETATM may be oddly renumbered."
                          %self.id)
            filter_by_ca = True
        
        # Calculate displacement value from 1st residue
        residue_list = self.get_list()
        displace = res_init - residue_list[0].id[1]
        
        # To keep track of last residue number
        # Used when renumbering consecutive chains sequentially 
        last_num = 0
        
        for residue in residue_list:
            # ATOMs and Non-HOH / Non-Ions / Non-ModResidues HETATMs
            if  (residue.id[0] == ' ') or \
                (not filter_by_ca and residue.resname in seqres) or \
                (filter_by_ca and 'CA' in residue.child_dict.keys() and residue['CA'].fullname == ' CA '):
                
                last_num = residue.id[1]+displace
                new_id = (residue.id[0], residue.id[1]+displace, residue.id[2])
            else: # HOH, Ions, and other HETATMs
                # het_init is 0 by default and not none or false
                # to allow numbering simple sum if renumbering several chains sequentially
                if het_init == 0:
                    het_init = 1000*(int(last_num/1000)+1) if last_num > 1000 else 1000
                new_id = (residue.id[0], het_init, residue.id[2])
                het_init += 1

            self._renumber_single_residue(residue.id, new_id)

        return (last_num, het_init)
            
