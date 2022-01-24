# -*- coding: utf-8 -*-
"""

PDB File representation

"""
import numpy

from .molecules import Molecule

_resSeq_min = 22
_resSeq_max = 26
_chainId_min = 21
_chainId_max = 22

class PDBFile:
    """Represents a PDB file with a protein-pigment complex structure

    """

    def __init__(self, fname=None):

        self.lines = []
        self.linecount = 0
        
        self.molecules = []

        self._resSeqs = []
        self._uniqueIds = []
        self._res_lines = dict()
        self._unique_lines = dict()

        if fname is not None:
            #load file
            self.linecount = self.load_file(fname)
        else:
            return

    def _reset_helpers(self):
        self._resSeqs = []
        self._uniqueIds = []
        self._res_lines = dict()
        self._unique_lines = dict()       
            
    def load_file(self, fname):
        """Loads a PDB file

        """
        with open(fname) as file:
            k = 0
            for line in file:
                self.lines.append(line)
                k += 1
        return k


    def get_Molecules(self, model=None):    
        """Returns all molecules corresponding to a given model

        """

        if model is None:
            return self.molecules

        else:
            
            molecules = []
            
            #
            #  residue name identifying the molecule in pdb file
            #
            res_name = model.pdbname
            
            #
            # Get all lines with molecules matching residue name
            #
            mollines = self._match_lines(by_recName="HETATM",
                                         by_resName=res_name)
                                      #by_resSeq=378, by_atmName="ND")

            #    return line[_resSeq_min:_resSeq_max]
            # Get resSeq - residue sequence number 
            # (sometimes enough to identify uniquely the molechaincule)
            #)
            for mols in mollines:
                rseq = int(line_resSeq(mols))
                if rseq not in self._resSeqs:
                    # if the sequence numeber not in a list of number
                    self._resSeqs.append(rseq) # append it
                    self._res_lines[rseq] = [] # make a new list of lines 
                self._res_lines[rseq].append(mols) # append line to a list
                                                   # according to a sequence nr
                                                   
            # check if there is one molecule in the list of lines

                
            count = -1
            for mols in mollines:
                rseq = int(line_resSeq(mols))
                chainId = line_chainId(mols)
                comb = chainId+str(rseq)
                if comb not in self._uniqueIds:
                    count += 1
                    self._uniqueIds.append(comb)
                    self._unique_lines[comb] = []
                self._unique_lines[comb].append(mols)

#            # Create a molecule from given lines
#            for rseq in self._resSeqs:
#                m = Molecule(name=str(rseq), elenergies=model.default_energies)
#                r = model.position_of_center(data_type="PDB",
#                                             data=self._res_lines[rseq])
#                m.position = r
#                d = model.transition_dipole(data_type="PDB",
#                                            data=self._res_lines[rseq])
#                m.set_dipole(0,1,d)
#                self.molecules.append(m)

            # Create a molecule from given unique lines
            for rseq in self._uniqueIds:
                m = Molecule(name=str(rseq), elenergies=model.default_energies)
                r = model.position_of_center(data_type="PDB",
                                             data=self._unique_lines[rseq])
                m.position = r
                d = model.transition_dipole(data_type="PDB",
                                            data=self._unique_lines[rseq])
                m.set_dipole(0,1,d)
                m.model = self
                m.data = self._unique_lines[rseq]
                molecules.append(m)
                
        self._reset_helpers()
        self.molecules = molecules
                
        return molecules

    def get_chainId(self,molecule):
        lines = molecule.data
        save = None
        for l in lines:
            cid = line_chainId(l)
            if save is not None:
                if save != cid:
                    raise Exception("No unique chainId")
            else:
                save = cid
        return save
                
    def clear_Molecules(self):
        self.molecules = []
        
    def _match_lines(self, by_recName=None,
                     by_resName=None,
                     by_resSeq=None,
                     by_atmName=None):
        """Matches a line with a given pattern

        """

        matched_lines = []
        k = 0
        for line in self.lines:
            if line_matches(line, by_recName=by_recName,
                            by_resName=by_resName, by_resSeq=by_resSeq,
                            by_atmName=by_atmName):
                matched_lines.append(line)
                k += 1
        #print("matched", k, "lines")
        return matched_lines


def line_resSeq(line):
    """Returns resSeq of the line

    """
    return line[_resSeq_min:_resSeq_max]
    
def line_chainId(line):
    """Returns chainId of a given line
    
    """
    return line[_chainId_min:_chainId_max]
    
def line_xyz(line):
    """Returns coordinates of the line

    """
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    return numpy.array([x,y,z])


def line_matches(line, 
                 by_recName=None,
                 by_resName=None,
                 by_chainId=None,
                 by_resSeq=None,
                 by_atmName=None):
        """Matches a line for several possible patters

        """
        ret = True
        if by_recName is not None:
            rname = line[0:6]
            #print(rname, by_recName)
            if rname.strip() == by_recName.strip():
                ret = ret and True
            else:
                return ret and False
        if by_atmName is not None:
            aname = line[12:16]
            if aname.strip() == by_atmName.strip():
                ret = ret and True
            else:
                return ret and False
        if by_resName is not None:
            resname = line[17:20]
            if resname == by_resName:
                ret = ret and True
            else:
                return False
        if by_chainId is not None:
            chainid = line[21:22]
            if chainid == by_chainId:
                ret = ret and True
            else:
                return False
        if by_resSeq is not None:
            resseq = line[22:26]
            if int(resseq) == by_resSeq:
                ret = ret and True
            else:
                return False
        return ret

