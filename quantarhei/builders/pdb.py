# -*- coding: utf-8 -*-
"""

PDB File representation

"""
import numpy

from .molecules import Molecule

_resSeq_min = 22
_resSeq_max = 26

class PDBFile:
    """Represents a PDB file with a protein-pigment complex structure

    """

    def __init__(self, fname=None):

        self.lines = []
        self.linecount = 0
        self.molecules = []

        self._resSeqs = []
        self._res_lines = dict()

        if fname is not None:
            #load file
            self.linecount = self.load_file(fname)
        else:
            return

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

            res_name = model.pdbname
            #
            # Get BCHLs
            #
            mollines = self._match_lines(by_recName="HETATM",
                                         by_resName=res_name)
                                      #by_resSeq=378, by_atmName="ND")

            #
            # Get resSed
            #
            for mols in mollines:
                rseq = int(line_resSeq(mols))
                if rseq not in self._resSeqs:
                    self._resSeqs.append(rseq)
                    self._res_lines[rseq] = []
                #altloc = mols[16]
                #if altloc.strip() != "":
                #    print(mols)
                self._res_lines[rseq].append(mols)

            # Create a molecule from given lines
            for rseq in self._resSeqs:
                m = Molecule(name=str(rseq), elenergies=model.default_energies)
                r = model.position_of_center(data_type="PDB",
                                             data=self._res_lines[rseq])
                m.position = r
                d = model.transition_dipole(data_type="PDB",
                                            data=self._res_lines[rseq])
                m.set_dipole(0,1,d)
                self.molecules.append(m)

        return self.molecules


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
        print("matched", k, "lines")
        return matched_lines


def line_resSeq(line):
    """Returns resSeq of the line

    """
    return line[_resSeq_min:_resSeq_max]

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
        if by_resSeq is not None:
            resseq = line[22:26]
            if int(resseq) == by_resSeq:
                ret = ret and True
            else:
                return False
        return ret

