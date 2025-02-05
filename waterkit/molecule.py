#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# WaterKit
#
# Class for molecule
#

import os
import re

import numpy as np
import openbabel as ob

from .typer import RotatableBonds
from .typer import HydrogenBonds
from . import utils


class Molecule():

    def __init__(self, OBMol, guess_hydrogen_bonds=True, guess_disordered_hydrogens=True):
        """Initialize a Molecule object.

        Args:
            OBMol (OBMol): OpenBabel molecule object
            guess_hydrogen_bonds (bool): hydrogen bonds typing (Default: True)
            guess_disordered_hydrogens (bool): disordered hydrogens typing (Default: True)

        """
        i = 0
        j = 0
        dtype = [("i", "i4"), ("name", "U4"), ("resname", "U3"), ("resnum", "i4"),
                 ('chain', 'U1'), ("xyz", "f4", (3)), ("q", "f4"), ("t", "U5")]
        self.atoms = np.zeros(OBMol.NumAtoms(), dtype)
        self.hydrogen_bonds = None
        self.rotatable_bonds = None

        # Remove all implicit hydrogens because OpenBabel
        # is doing chemical perception, and we want to read the
        # molecule as is.
        for x in ob.OBMolAtomIter(OBMol):
            if x.GetAtomicNum() != 1 and x.GetImplicitHCount() != 0:
                x.SetImplicitHCount(0)

        for r in ob.OBResidueIter(OBMol):
            for a in ob.OBResidueAtomIter(r):
                xyz = (a.GetX(), a.GetY(), a.GetZ())
                atom = (i+1, r.GetAtomID(a).strip(), r.GetName(), r.GetNum(), r.GetChain(), 
                        xyz, a.GetPartialCharge(), a.GetType())
                self.atoms[i] = atom
                
                i += 1
            j += 1

        # Do the typing for hydrogen bonds and disordered hydrogen atoms
        d = utils.path_module("waterkit")
        
        if guess_disordered_hydrogens:
            dh_file = os.path.join(d, "data/disordered_hydrogens.par")
            dhfield = RotatableBonds(dh_file)
            self._guess_rotatable_bonds(OBMol, dhfield)

        if guess_hydrogen_bonds:
            hb_file = os.path.join(d, "data/waterfield.par")
            hbfield = HydrogenBonds(hb_file)
            self._guess_hydrogen_bond_anchors(OBMol, hbfield)

    @classmethod
    def from_file(cls, fname, guess_hydrogen_bonds=True, guess_disordered_hydrogens=True):
        """Create Molecule object from a PDB file.

        Args:
            fname (str): molecule filename
            guess_hydrogen_bonds (bool): hydrogen bonds typing (Default: True)
            guess_disordered_hydrogens (bool): disordered hydrogens typing (Default: True)

        Returns:
            Molecule

        """
        # Get name and file extension
        name, file_extension = os.path.splitext(fname)
        file_extension = file_extension.split(os.extsep)[-1]

        # Read PDB file
        OBMol = ob.OBMol()
        obconv = ob.OBConversion()

        if file_extension == "pdbqt":
            obconv.SetInFormat("pdb")
        else:
            obconv.SetInFormat(file_extension)
        
         # set error level to avoid warnings about non-standard input
        errlev = ob.obErrorLog.GetOutputLevel()
        ob.obErrorLog.SetOutputLevel(0)

        obconv.ReadFile(OBMol, fname)

        ob.obErrorLog.SetOutputLevel(errlev)

        m = cls(OBMol, guess_hydrogen_bonds, guess_disordered_hydrogens)

        # OpenBabel do chemical perception to define the type
        # So we override the types with AutoDock atom types
        # from the PDBQT file
        if file_extension == "pdbqt":
            qs, ts = m._qt_from_pdbqt_file(fname)
            m.atoms['q'] = qs
            m.atoms['t'] = ts

        return m

    def _qt_from_pdbqt_file(self, fname):
        """Get partial charges and atom types from PDBQT file.

        Args:
            fname (str): molecule filename

        Returns:
            list: partial charges
            list: AutoDock atom types

        """
        atom_types = []
        partial_charges = []

        with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                if re.search("^ATOM", line) or re.search("^HETATM", line):
                    atom_types.append(line[77:79].strip())
                    partial_charges.append(float(line[70:77].strip()))

        return partial_charges, atom_types

    def is_water(self):
        """Tell if it is a water or not."""
        return False

    def atom(self, atom_id):
        """Return the atom i

        Args:
            atom_id (int): atom id
        
        Returns:
            ndarray: 1d ndarray (i, name, resname, resnum, xyz, q, t)

        """
        return self.atoms[i - 1]

    def coordinates(self, atom_ids=None):
        """
        Return coordinates of all atoms or a certain atom
        
        Args:
            atom_ids (int, list): index of one or multiple atoms

        Returns:
            ndarray: 2d ndarray of 3d coordinates

        """
        if atom_ids is not None and self.atoms.size > 1:
            if not isinstance(atom_ids, np.ndarray):
                atom_ids = np.array(atom_ids, dtype=int)
            # -1 because numpy array is 0-based
            atoms = self.atoms[atom_ids - 1]["xyz"]
        else:
            atoms = self.atoms["xyz"]

        return np.atleast_2d(atoms).copy()

    def atom_types(self, atom_ids=None):
        """Return atom types of all atoms or a certain atom.
        
        Args:
            atom_ids (int, list): index of one or multiple atoms

        Returns:
            list: atom types

        """
        if atom_ids is not None and self.atoms.size > 1:
            if not isinstance(atom_ids, np.ndarray):
                atom_ids = np.array(atom_ids, dtype=int)
            # -1 because numpy array is 0-based
            t = self.atoms[atom_ids - 1]['t']
        else:
            t = self.atoms['t']

        return t.tolist()

    def partial_charges(self, atom_ids=None):
        """Get partial charges.

        Args:
            atom_ids (int, list): index of one or multiple atoms

        Returns:
            ndarray: partial charges

        """
        if atom_ids is not None and self.atoms.size > 1:
            if not isinstance(atom_ids, np.ndarray):
                atom_ids = np.array(atom_ids, dtype=int)
            # -1 because numpy array is 0-based
            q = self.atoms[atom_ids - 1]['q']
        else:
            q =  self.atoms['q']

        return q.copy()

    def atom_informations(self, atom_ids=None):
        """Get atom informations (xyz, q, type).
        
        Args:
            atom_ids (int, list): index of one or multiple atoms

        Returns:
            ndarray: atom information (i, xyz, q, t)

        """
        if atom_ids is not None:
            if not isinstance(atom_ids, np.ndarray):
                atom_ids = np.array(atom_ids, dtype=int)
            # -1 because numpy array is 0-based
            atoms = self.atoms[atom_ids - 1][['i', 'xyz', 'q', 't']]
        else:
            atoms = self.atoms[['i', 'xyz', 'q', 't']]

        return atoms.copy()

    def _guess_rotatable_bonds(self, OBMol, rotatable_bond_typer):
        """Guess all the rotatable bonds in the molecule
        based the rotatable forcefield.
        
        Args:
            OBMol (OBMol): OBMolecule object
            rotatable_bond_typer (RotatableBonds): RotatableBonds object

        """
        self.rotatable_bonds = rotatable_bond_typer.match(OBMol)

    def _guess_hydrogen_bond_anchors(self, OBMol, hydrogen_bond_typer):
        """Guess all the hydrogen bonds in the molecule
        based the hydrogen bond forcefield.
        
        Args:
            OBMol (OBMol): OBMolecule object
            hydrogen_bond_typer (HydrogenBonds): HydrogenBonds object

        """
        self.hydrogen_bonds = hydrogen_bond_typer.match(OBMol)

    def update_coordinates(self, xyz, atom_id):
        """Update the coordinates of an atom.

        Args:
            xyz (array_like): 3d coordinates of the new atomic position
            atom_id (int): atom id

        Returns:
            bool: True if successfull, False otherwise

        """
        if atom_id <= self.atoms.size:
            if self.atoms.size > 1:
                self.atoms[atom_id - 1]["xyz"] = xyz
            else:
                self.atoms["xyz"] = xyz
            return True
        else:
            return False

    def to_pdbqt_file(self, fname):
        """Write molecule to PDBQT file

        Args:
            fname (str): output PDBQT filename

        """
        pdbqt_str = "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  1.00    %6.3f %-2s\n"
        output_str = ""

        if self.atoms.size == 1:
            atoms = [self.atoms]
        else:
            atoms = self.atoms

        for atom in atoms:
            x, y, z = atom["xyz"]
            output_str += pdbqt_str % (atom["i"], atom["name"], atom["resname"], atom['chain'], atom["resnum"], 
                                       x, y, z, atom["q"], atom["t"])

        with open(fname, 'w') as w:
            w.write(output_str)


    def to_file(self, fname, fformat, options=None, append=False):
        """Write molecule to file using OpenBabel.
    
        OpenBabel is used to write in format other than PDBQT file because
        it is using the default AutoDock atom types (OA, HD, ...). And so
        it is not working for our special water atom types.

        Args:
            fname (str): output PDBQT filename
            fformat (str): molecule file format
            options (str): OpenBabel wrting options
            append (bool): append to existing PDBQT file (default: False)

        """
        pdbqt_str = "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  1.00    %6.3f %-2s\n"
        output_str = ""

        if self.atoms.size == 1:
            atoms = [self.atoms]
        else:
            atoms = self.atoms

        for atom in atoms:
            x, y, z = atom["xyz"]
            output_str += pdbqt_str % (atom["i"], atom["name"], atom["resname"], atom['chain'], atom["resnum"],
                                       x, y, z, atom["q"], atom["t"])

        if fformat != "pdbqt":
            OBMol = ob.OBMol()
            obconv = ob.OBConversion()
            obconv.SetInAndOutFormats("pdbqt", fformat)
            obconv.ReadString(OBMol, output_str)

            if options is not None:
                for option in options:
                    obconv.AddOption(option)

            output_str = obconv.WriteString(OBMol)

        if append and os.path.isfile(fname):
            with open(fname, "a") as a:
                a.write(output_str)
        else:
            with open(fname, "w") as w:
                w.write(output_str)

    def export_hb_vectors(self, fname):
        """Export all the hb vectors to PDBQT file. 
        
        Args:
            fname (str): filename

        Returns:
            None

        """
        pdbqt_str = "ATOM  %5d  %-3s ANC%2s%4d    %8.3f%8.3f%8.3f%6.2f 1.00    %6.3f %2s\n"

        if self.hydrogen_bonds is not None:
            i = 1
            output_str = ""

            for index, anchor in self.hydrogen_bonds.iterrows():
                x, y, z = anchor.vector_xyz
                atom_type = anchor.anchor_type[0].upper()

                output_str += pdbqt_str % (i, atom_type, "A", index, x, y, z, 1, 1, atom_type)
                i += 1

            with open(fname, "w") as w:
                w.write(output_str)
        else:
            print("Error: There is no hydrogen bond anchors.")
