#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# prepare water cluster
#

import argparse
import copy
import logging
import os
import shutil
import string

import parmed as pmd
from pdb4amber import AmberPDBFixer
from pdb4amber.utils import easy_call

from waterkit import utils

RESSOLV = ('WAT', 'HOH')

AMBER_SUPPORTED_RESNAMES = RESSOLV


def _write_pdb_file(output_name, molecule, overwrite=True, **kwargs):
    '''Write PDB file

    Args:
        output_name (str): pdbqt output filename
        molecule (parmed): parmed molecule object

    '''
    i = 0
    pdb_template = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:>4s}{:>2s}{:2s}\n"
    ter_template = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}\n"
    output_str = ''

    for residue in molecule.residues:
        resname = residue.name
        resid = residue.number
        chain_id = residue.chain

        for atom in residue.atoms:
            if len(atom.name) < 4:
                name = ' %s' % atom.name
            else:
                name = atom.name

            atom_type = atom.name[0]

            if resname in RESSOLV:
                atype = 'HETATM'
            else:
                atype = 'ATOM'

            output_str += pdb_template.format(atype, i + 1, name, " ", resname, chain_id, resid,
                                              " ", atom.xx, atom.xy, atom.xz, 1.0, 1.0, " ",
                                              atom_type, " ")

            i += 1

        if residue.ter:
            output_str += ter_template.format('TER', i + 1, name, " ", resname, chain_id, resid)
            i += 1

    output_str += 'END\n'

    with open(output_name, 'w') as w:
        w.write(output_str)


def _write_pdbqt_file(output_name, molecule):
    '''Write PDBQT file

    Args:
        output_name (str): pdbqt output filename
        molecule (parmed): parmed molecule object

    '''
    i = 0
    pdb_template = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:4s}{:6.3f} {:2s} \n"
    ter_template = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}\n"
    output_str = ''

    for residue in molecule.residues:
        resname = residue.name
        resid = residue.number
        chain_id = residue.chain

        for atom in residue.atoms:
            if len(atom.name) < 4:
                name = ' %s' % atom.name
            else:
                name = atom.name

            # OpenBabel does not like atom types starting with a number
            if atom.type[0].isdigit():
                atom_type = atom.type[::-1]
            else:
                atom_type = atom.type

            # AutoGrid does not accept atom type name of length > 2
            atom_type = atom_type[:2]

            if resname in RESSOLV:
                atype = 'HETATM'
            else:
                atype = 'ATOM'

            output_str += pdb_template.format(atype, i + 1, name, " ", resname, chain_id, resid,
                                              " ", atom.xx, atom.xy, atom.xz, 1.0, 1.0, " ",
                                              atom.charge, atom_type)

            i += 1

        if residue.ter:
            output_str += ter_template.format('TER', i + 1, name, " ", resname, chain_id, resid)
            i += 1

    output_str += 'END\n'

    with open(output_name, 'w') as w:
        w.write(output_str)s


def _convert_amber_to_autodock_types(molecule):
    molecule = copy.deepcopy(molecule)

    amber_autodock_dict = {
        'N3': 'N',
        'H': 'HD',
        'CX': 'C',
        'HP': 'H',
        'CT': 'C',
        'HC': 'H',
        'C': 'C',
        'O': 'OA',
        'N': 'N',
        'H1': 'H',
        'C3': 'C',
        '3C': 'C',
        'C2': 'C',
        '2C': 'C',
        'CO': 'C',
        'O2': 'OA',
        'OH': 'OA',
        'HO': 'HD',
        'SH': 'SA',
        'HS': 'HD',
        'CA': 'A',
        'HA': 'H',
        'S': 'SA',
        'C8': 'C',
        'N2': 'N',
        'CC': 'A',
        'NB': 'NA',
        'CR': 'A',
        'CV': 'A',
        'H5': 'H',
        'NA': 'N',
        'CW': 'A',
        'H4': 'H',
        'C*': 'A',
        'CN': 'A',
        'CB': 'A',
        'Zn2+': 'Zn',
        'Zn': 'Zn',
        'Mn2+': 'Mn',
        'Mn': 'Mn',
        'XC': 'C',
        'br': 'Br',
        'c': 'C',
        'c1': 'C',
        'c2': 'C',
        'c3': 'C',
        'ca': 'A',
        'cc': 'A',
        'cd': 'A',
        'ce': 'C',
        'cf': 'C',
        'cl': 'Cl',
        'cp': 'A',
        'cq': 'A',
        'cu': 'C',
        'cv': 'C',
        'cx': 'C',
        'cy': 'C',
        'c5': 'C',
        'c6': 'C',
        'cz': 'C',
        'cs': 'C',
        'cg': 'C',
        'ch': 'C',
        'f': 'F',
        'h1': 'H',
        'h2': 'H',
        'h3': 'H',
        'h4': 'H',
        'h5': 'H',
        'ha': 'H',
        'hc': 'H',
        'hn': 'HD',
        'ho': 'HD',
        'hp': 'HD',
        'hs': 'HD',
        'hx': 'H',
        'i': 'I',
        'n': 'N',
        'n1': 'NA',
        'n2': 'N',
        'n3': 'N',
        'n4': 'N',
        'n5': 'N',
        'n6': 'N',
        'n7': 'N',
        'n8': 'N',
        'n9': 'N',
        'na': 'N',
        'nb': 'N',
        'nc': 'N',
        'nd': 'N',
        'nh': 'N',
        'ne': 'N',
        'nf': 'N',
        'no': 'N',
        'n+': 'N',
        'nx': 'N',
        'ny': 'N',
        'nz': 'N',
        'ns': 'N',
        'nt': 'N',
        'nu': 'N',
        'nv': 'N',
        'ni': 'N',
        'nj': 'N',
        'nk': 'N',
        'nl': 'N',
        'nm': 'N',
        'nn': 'N',
        'np': 'N',
        'nq': 'N',
        'o': 'OA',
        'oh': 'OA',
        'os': 'OA',
        'op': 'OA',
        'oq': 'OA',
        'p2': 'P',
        'p3': 'P',
        'p4': 'P',
        'p5': 'P',
        'pb': 'P',
        'pc': 'P',
        'pd': 'P',
        'pe': 'P',
        'pf': 'P',
        'px': 'P',
        'py': 'P',
        's': 'S',
        's2': 'SA',
        's4': 'S',
        's6': 'S',
        'sh': 'SA',
        'ss': 'SA',
        'sx': 'S',
        'sy': 'S',
        'sp': 'S',
        'sq': 'S',
        'OW': 'O',  # Add this line
        'HW': 'H'   # Add this line
        # 'Cu':
    }

    for atom in molecule.atoms:
        if atom.residue.name == 'TYR' and atom.name == 'CZ' and atom.type == 'C':
            atom.type = 'A'
        elif atom.residue.name == 'ARG' and atom.name == 'CZ' and atom.type == 'CA':
            atom.type = 'C'
        else:
            atom.type = amber_autodock_dict[atom.type]

    return molecule


def _make_leap_template(parm, input_pdb,
                        prmtop='prmtop', rst7='rst7', lib_files=None, frcmod_files=None):
    default_force_field = ('source leaprc.water.tip3p\n')

    leap_template = ('{force_fields}\n'
                     '{more_force_fields}\n'
                     'x = loadpdb {input_pdb}\n'
                     '{box_info}\n'
                     '{more_leap_cmds}\n'
                     'set default nocenter on\n'
                     'saveAmberParm x {prmtop} {rst7}\n'
                     'quit\n')

    box = parm.box
    if box is not None:
        box_info = 'set x box { %s  %s  %s }' % (box[0], box[1], box[2])
    else:
        box_info = ''

    leap_string = leap_template.format(
        force_fields=default_force_field,
        more_force_fields='',
        box_info=box_info,
        input_pdb=input_pdb,
        prmtop=prmtop,
        rst7=rst7,
        more_leap_cmds='')

    return leap_string


def _generate_resids_chainids_for_tleap(molecule):
    resids = list(range(1, len(molecule.residues) + 1))
    chainids = [" "] * len(resids)
    return resids, chainids


def _ter_flags(molecule):
    return [True if residue.ter else False for residue in molecule.residues]


def _replace_resids_and_chainids(molecule, new_resids, new_chainids):
    n_resids = len(molecule.residues)
    n_new_resids = len(new_resids)

    error_msg = 'Cannot replace resids and chainids.'
    error_msg += ' Number of residues is different (%d != %d).' % (n_resids, n_new_resids)
    assert n_resids == n_new_resids, error_msg

    for i in range(n_resids):
        molecule.residues[i].number = new_resids[i]
        molecule.residues[i].chain = new_chainids[i]


def _add_ter_flags(molecule, ter_flags):
    n_resids = len(molecule.residues)
    n_ter_flags = len(ter_flags)

    error_msg = 'Cannot add TER flags.'
    error_msg += ' Number of TER flags and resids are different (%d != %d).' % (n_resids, n_ter_flags)
    assert n_resids == n_ter_flags, error_msg

    for i in range(n_resids):
        molecule.residues[i].ter = ter_flags[i]


class PrepareWaterCluster:

    def __init__(self, keep_hydrogen=False, keep_altloc=False, ignore_gaps=False, renumbering=True, use_model=1):
        self._keep_hydrogen = keep_hydrogen
        self._keep_altloc = keep_altloc
        self._use_model = use_model
        self._fill_gaps = ignore_gaps
        self._renumbering = renumbering

        self._pdb_filename = None
        self._molecule = None

    def prepare(self, pdb_filename, clean=True):
        '''Prepare water cluster structure

        Args:
            pdb_filename (str): input pdb filename
            clean (bool): remove tleap input and output files (default: True)

        '''
        tleap_input = 'leap.template.in'
        tleap_output = 'leap.template.out'
        tleap_log = 'leap.log'
        pdb_clean_filename = 'tmp.pdb'
        prmtop_filename = 'tmp.prmtop'
        rst7_filename = 'tmp.rst7'

        logger = logging.getLogger('WaterKit water cluster preparation')
        logging.basicConfig(level=os.environ.get('LOGLEVEL', 'INFO'))

        try:
            water_cluster = pmd.load_file(pdb_filename)
        except FileNotFoundError:
            error_msg = 'Water cluster file (%s) cannot be found.' % pdb_filename
            logger.error(error_msg)
            raise

        pdbfixer = AmberPDBFixer(water_cluster)

        # Remove box and symmetry
        pdbfixer.parm.box = None
        pdbfixer.parm.symmetry = None

        # Remove all the hydrogens
        if not self._keep_hydrogen:
            pdbfixer.parm.strip('@/H')
            logger.info('Removed all hydrogen atoms')

        # Keep only water residues
        pdbfixer.parm.strip('!:' + ','.join(RESSOLV))
        logger.info('Removed all non-water residues')

        # Renumbers resids (starting from 1) and renames chainids (starting from A)
        if self._renumbering:
            new_resids, new_chainids = _generate_resids_chainids_for_tleap(pdbfixer.parm)
            _replace_resids_and_chainids(pdbfixer.parm, new_resids, new_chainids)
        else:
            new_resids = [residue.number for residue in pdbfixer.parm.residues]
            new_chainids = [residue.chain for residue in pdbfixer.parm.residues]

        # Take the first model only
        final_coordinates = pdbfixer.parm.get_coordinates()[self._use_model - 1]
        write_kwargs = dict(coordinates=final_coordinates)
        write_kwargs['increase_tercount'] = False  # so CONECT record can work properly
        write_kwargs['altlocs'] = 'occupancy'

        # Create Amber topology and coordinates files
        with utils.temporary_directory(prefix='wk_preparation_', dir='.', clean=clean) as tmp_dir:
            tleap_resids, tleap_chainids = _generate_resids_chainids_for_tleap(pdbfixer.parm)
            _replace_resids_and_chainids(pdbfixer.parm, tleap_resids, tleap_chainids)

            try:
                _write_pdb_file(pdb_clean_filename, pdbfixer.parm, **write_kwargs)
            except:
                error_msg = 'Could not write pdb file %s' % pdb_clean_filename
                logger.error(error_msg)
                raise

            # Generate topology/coordinates files
            with open(tleap_input, 'w') as w:
                content = _make_leap_template(pdbfixer.parm, input_pdb=pdb_clean_filename,
                                              prmtop=prmtop_filename, rst7=rst7_filename)
                w.write(content)

            try:
                easy_call('tleap -s -f %s > %s' % (tleap_input, tleap_output), shell=True)
            except RuntimeError:
                error_msg = 'Could not generate topology/coordinates files with tleap'
                logger.error(error_msg)
                raise RuntimeError(error_msg)

            self._molecule = pmd.load_file(prmtop_filename, rst7_filename)

            # Add back resids, chainids and TER flags to molecule
            _replace_resids_and_chainids(self._molecule, new_resids, new_chainids)
            _add_ter_flags(self._molecule, _ter_flags(pdbfixer.parm))

    def write_pdb_file(self, pdb_filename='water_cluster.pdb'):
        _write_pdb_file(pdb_filename, self._molecule)

    def write_pdbqt_file(self, pdbqt_filename='water_cluster.pdbqt', amber_atom_types=False):
        if amber_atom_types:
            _write_pdbqt_file(pdbqt_filename, self._molecule)
        else:
            molecule = _convert_amber_to_autodock_types(self._molecule)
            _write_pdbqt_file(pdbqt_filename, molecule)


def cmd_lineparser():
    parser = argparse.ArgumentParser(description='prepare water cluster')
    parser.add_argument('-i', '--in', required=True,
                        dest='pdb_filename', help='PDB input file (default: stdin)',
                        default='stdin')
    parser.add_argument('-o', '--out', default='water_cluster',
                        dest='output_prefix', help='output prefix filename (default: water_cluster)')
    parser.add_argument('--keep_hydrogen', action='store_true', default=False,
                        dest='keep_hydrogen', help='keep all hydrogen atoms (default: no)')
    parser.add_argument('--keep_altloc', action='store_true', default=False,
                        dest='keep_altloc', help='keep residue altloc (default is to keep "A")')
    parser.add_argument('--model', type=int, default=1,
                        dest='use_model',
                        help='Model to use from a multi-model pdb file (integer).  (default: use 1st model). '
                             'Use a negative number to keep all models')
    parser.add_argument('--pdb', dest='make_pdb', default=False,
                        action='store_true', help='generate pdb file')
    parser.add_argument('--pdbqt', dest='make_pdbqt', default=False,
                        action='store_true', help='PDBQT file with AutoDock atom types')
    parser.add_argument('--amber_pdbqt', dest='make_amber_pdbqt', default=False,
                        action='store_true', help='PDBQT file with Amber atom types')
    parser.add_argument('--ignore_gaps', dest='ignore_gaps', default=False,
                        action='store_true', help='ignore gaps between residues (automatically add TER records)')
    parser.add_argument('--renumber', dest='renumbering', default=False,
                        action='store_true', help='Residue index will be renumbered (starting from 1).')
    parser.add_argument('--no_clean', dest='no_clean', default=True,
                        action='store_false', help='Does not clean Amber temporay files.')
    return parser.parse_args()


def main():
    args = cmd_lineparser()
    pdb_filename = args.pdb_filename
    output_prefix = args.output_prefix
    keep_hydrogen = args.keep_hydrogen
    keep_altloc = args.keep_altloc
    use_model = args.use_model
    make_pdb = args.make_pdb
    make_pdbqt = args.make_pdbqt
    make_amber_pdbqt = args.make_amber_pdbqt
    ignore_gaps = args.ignore_gaps
    renumbering = args.renumbering
    clean = -args.no_clean

    pr = PrepareWaterCluster(keep_hydrogen, keep_altloc, ignore_gaps, renumbering, use_model)
    pr.prepare(pdb_filename, clean)

    if make_pdb:
        pdb_prepared_filename = '%s.pdb' % output_prefix
        pr.write_pdb_file(pdb_prepared_filename)

    if make_pdbqt:
        pdbqt_prepared_filename = '%s.pdbqt' % output_prefix
        pr.write_pdbqt_file(pdbqt_prepared_filename)

    if make_amber_pdbqt:
        pdbqt_prepared_filename = '%s_amber.pdbqt' % output_prefix
        pr.write_pdbqt_file(pdbqt_prepared_filename, amber_atom_types=True)


if __name__ == '__main__':
    main()
