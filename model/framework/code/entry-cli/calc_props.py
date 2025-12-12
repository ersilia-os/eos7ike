###############################
import argparse
import os
import sys
import csv
import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
###############################

__doc__="""Performs calculation of physiochemical properties of potential antibiotics. SMILES strings are parsed,
conformers are generated using RDKit, and properties calculated. Properties include: chemical formula, molecular weight, rotatable
bonds, globularity, and PBF.
"""

# SMARTS definition for functional groups (using pybel.Smarts)
FUNCTIONAL_GROUP_TO_SMARTS = {
    'primary_amine': pybel.Smarts('[$([N;H2;X3][CX4]),$([N;H3;X4+][CX4])]')
}
FUNCTIONAL_GROUPS = sorted(FUNCTIONAL_GROUP_TO_SMARTS.keys())


def main():
    args = parse_args(sys.argv[1:])
    if(args.smiles):
        # Use new RDKit-based function
        mol = smiles_to_rdkit(args.smiles)
        properties = average_properties(mol)
        properties['smiles'] = args.smiles
        # A file will be written if command line option provide, otherwise write to stdout
        if(args.output):
            mols_to_write = [properties]
            write_csv(mols_to_write, args.output)
        else:
            report_properties(properties)
    elif(args.batch_file):
        mols = parse_batch(args.batch_file)
        mols_to_write = []
        for smiles, name in mols:
            # Use new RDKit-based function
            mol = smiles_to_rdkit(smiles)
            properties = average_properties(mol)
            properties['smiles'] = name
            mols_to_write.append(properties)
        write_csv(mols_to_write, args.output)


def parse_args(arguments):
    """Parse the command line options.
    :return:  All script options
    """
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-s", "--smiles", dest="smiles", metavar="SMILES string", default=None)
    group.add_argument("-b", "--batch", dest="batch_file", metavar="Batch file", default=None)
    parser.add_argument("-o", "--output", dest="output", metavar="Output file", default=None,
                         help="Defaults to csv file with same name as input")

    args = parser.parse_args(arguments)
    if not args.smiles and not args.batch_file:
        parser.error("Input structure is needed")

    # If no output file is specified in batch mode, then replace the file extension of the input with .csv
    if args.batch_file and not args.output:
        args.output = os.path.splitext(args.batch_file)[0] + '.csv'

    return args


def report_properties(properties):
    """
    Write out the results of physiochemical properties to stdout

    :param properties: physiochemical properties to report
    :type properties: dict
    :return: None
    """
    print("Properties for %s" % properties['smiles'])
    print("--------------------------")
    print("Mol. Wt.:\t%f" % properties['molwt'])
    print("Formula:\t%s" % properties['formula'])
    print("RB:\t\t%i" % properties['rb'])
    print("Glob:\t\t%f" % properties['glob'])
    print("PBF:\t\t%f" % properties['pbf'])
    for functional_group in FUNCTIONAL_GROUPS:
        print("%s:\t%s" % (functional_group, properties[functional_group]))


def parse_batch(filename):
    """
    Read a file containing names and SMILES strings

    Expects a file with no header in which each line contains a SMILES string followed by a name for the molecule.
    SMILES and name can be separated by any whitespace.

    :param filename: file to read
    :type filename: str
    :return: List of tuples with names and SMILES
    :rtype: list
    """
    smiles = []
    names = []
    with(open(filename, 'r')) as batch:
        for line in batch:
            (smi, name) = tuple(line.split())
            smiles.append(smi)
            names.append(name)
    return zip(smiles, names)


def write_csv(mols_to_write, filename):
    """
    Write out results of physiochemical properties

    :param mols_to_write: list of molecule properties to write
    :param filename: path to file to write
    :type mols_to_write: list
    :type filename: str
    :return: None
    """
    with(open(filename, 'w')) as out:
        fieldnames = ['smiles', 'formula', 'molwt', 'rb', 'glob', 'pbf'] + FUNCTIONAL_GROUPS
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for mol in mols_to_write:
            writer.writerow(mol)


def average_properties(mol):
    """
    Calculate all relevant properties for a given molecule averaged across conformers.
    The molecule is an RDKit Mol object, which is converted to pybel/OpenBabel
    for some property calculations (molwt, formula, rb, functional groups).

    :param mol: input molecule
    :type mol: rdkit.Chem.Mol
    :return: dictionary of properties
    :rtype dict
    """
    # 1. Generate conformers using RDKit's ETKDG/MMFF
    mol = run_confab(mol)
    num_confs = mol.GetNumConformers()

    if num_confs == 0:
        print("Warning: No conformers generated. Cannot calculate 3D properties.")
        return {
            'formula': 'N/A', 'molwt': 0.0, 'rb': 0, 'glob': -1.0, 'pbf': -1.0,
            **{fg: False for fg in FUNCTIONAL_GROUPS}
        }

    globs = np.empty(num_confs)
    pbfs = np.empty(num_confs)

    # 2. Bridge to pybel/OpenBabel for non-3D properties
    # Convert the RDKit mol with conformers into a MolBlock, then read by OpenBabel
    rdkit_mol_block = Chem.MolToMolBlock(mol, confId=-1) # -1 exports the structure (without conformer info needed here)
    obmol = ob.OBMol()
    obConv = ob.OBConversion()
    obConv.SetInAndOutFormats("mol", "mol")
    obConv.ReadString(obmol, rdkit_mol_block)
    pymol = pybel.Molecule(obmol) # This pybel mol will be used for formula, molwt, rb, and smarts.

    # 3. Iterate over RDKit Conformers for 3D properties
    for i in range(num_confs):
        conf = mol.GetConformer(i)

        # Extract coordinates for the current RDKit conformer
        points = np.empty(shape=(mol.GetNumAtoms(), 3))
        for atom_idx in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            points[atom_idx] = (pos.x, pos.y, pos.z)

        # calculate 3D properties
        globs[i] = calc_glob(points)
        pbfs[i] = calc_pbf(points)

    data = {
        'formula': pymol.formula,
        'molwt': pymol.molwt,
        'rb': rotatable_bonds(pymol),
        'glob': np.mean(globs),
        'pbf': np.mean(pbfs)
    }

    # 4. Functional Groups
    for functional_group, smarts in FUNCTIONAL_GROUP_TO_SMARTS.items():
        data[functional_group] = has_functional_group(pymol, smarts)

    return data


def smiles_to_rdkit(mol_string):
    """
    Reads a SMILES string and creates an RDKit molecule object with added Hydrogens.
    This prepares the molecule for RDKit-based conformer generation.

    :param mol_string: SMILES string
    :type mol_string: str
    :return: RDKit molecule object with Hs
    :rtype: rdkit.Chem.Mol
    """
    m = Chem.MolFromSmiles(mol_string)
    if m is None:
        raise ValueError(f"Could not parse SMILES: {mol_string}")
    m_with_h = Chem.AddHs(m)
    return m_with_h


def run_confab(mol, num_confs=50, max_attempts=500, prune_rms=0.5):
    """
    Generate ensemble of conformers using RDKit's ETKDG method followed by MMFF optimization.

    :param mol: initial molecule (RDKit Mol) to generate conformers from
    :param num_confs: max number of conformers to generate
    :param max_attempts: max attempts for embedding
    :param prune_rms: RMSD cutoff for pruning redundant conformers, default: 0.5
    :type mol: rdkit.Chem.Mol
    :type num_confs: int
    :type prune_rms: float
    :return: molecule object with generated and minimized conformers
    :rtype: rdkit.Chem.Mol
    """
    # 1. Generate conformers using RDKit's ETKDG method
    # 
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, maxAttempts=max_attempts,
                                     pruneRmsThresh=prune_rms, randomSeed=42)

    # 2. Perform energy minimization for each conformer using MMFF
    minimized_mol = Chem.Mol(mol)
    for cid in cids:
        # Use MMFF for minimization
        AllChem.MMFFOptimizeMolecule(minimized_mol, confId=int(cid))

    return minimized_mol


def calc_glob(points):
    """
    Calculates the globularity (glob) of a molecule

    glob varies from 0 to 1 with completely flat molecules like benzene having a
    glob of 0 and spherical molecules like adamantane having a glob of 1

    :param points: numpy array of atom coordinates
    :type points: numpy array
    :return: globularity of molecule
    :rtype: float | int
    """
    if points is None or points.shape[0] < 3:
        return 0
    points = points.T

    # calculate covariance matrix
    cov_mat = np.cov([points[0,:],points[1,:],points[2,:]])

    # calculate eigenvalues of covariance matrix and sort
    vals, vecs = np.linalg.eig(cov_mat)
    vals = np.sort(vals)[::-1]

    # glob is ratio of last eigenvalue and first eigenvalue
    # Globularity = lambda_3 / lambda_1
    if vals[0] > 1e-6: # Check for near-zero to avoid division by zero
        return vals[-1]/vals[0]
    else:
        return 0


def calc_pbf(points):
    """
    Uses SVD to fit atoms in molecule to a plane then calculates the average
    distance to that plane.

    :param points: numpy array of atom coordinates
    :type points: numpy array
    :return: average distance of all atoms to the best fit plane
    :rtype: float
    """
    c, n = svd_fit(points)
    pbf = calc_avg_dist(points, c, n)
    return pbf


def has_functional_group(mol, smarts):
    """
    Determines whether the molecule has the functional group specified by the SMARTS.

    :param mol: pybel molecule object
    :param smarts: pybel SMARTS object
    :return: True if mol has an instance of the functional group, False otherwise
    :rtype: bool
    """
    functional_groups = smarts.findall(mol)

    return len(functional_groups) > 0


def rotatable_bonds(mol):
    """
    Calculates the number of rotatable bonds in a molecules. Rotors are defined
    as any non-terminal bond between heavy atoms, excluding amides

    :param mol: pybel molecule object
    :type mol: pybel.Molecule
    :return rb: number of rotatable bonds
    :rtype int
    """
    rb = 0
    for bond in ob.OBMolBondIter(mol.OBMol):
        if is_rotor(bond):
            rb += 1
    return rb


def is_rotor(bond, include_amides=False):
    """
    Determines if a bond is rotatable

    Rules for rotatable bonds:
    Must be a single or triple bond
    Must include two heavy atoms
    Cannot be terminal
    Cannot be in a ring
    If a single bond to one sp hybridized atom, not rotatable

    :param bond:
    :return: If a bond is rotatable
    :rtype: bool
    """
    # Must be single or triple bond
    if bond.GetBondOrder() == 2: return False

    # Don't count the N-C bond of amides
    if bond.IsAmide() and not include_amides: return False

    # Not in a ring
    if bond.FindSmallestRing() is not None: return False

    # Don't count single bonds adjacent to triple bonds, still want to count the triple bond
    if (bond.GetBeginAtom().GetHyb() == 1) != (bond.GetEndAtom().GetHyb() == 1): return False

    # Cannot be terminal (must be heavy atom degree > 1 on both sides)
    if bond.GetBeginAtom().GetHvyDegree() > 1 and bond.GetEndAtom().GetHvyDegree() > 1: return True
    
    return False


def calc_avg_dist(points, C, N):
    """
    Calculates the average distance a given set of points is from a plane

    :param points: numpy array of points
    :param C: centroid vector of plane
    :param N: normal vector of plane
    :return Average distance of each atom from the best-fit plane
    """
    sum_dist = 0
    for xyz in points:
        sum_dist += abs(distance(xyz, C, N))
    return sum_dist / len(points)


def svd_fit(X):
    """
    Fitting algorithmn was obtained from https://gist.github.com/lambdalisue/7201028
    Find (n - 1) dimensional standard (e.g. line in 2 dimension, plane in 3
    dimension, hyperplane in n dimension) via solving Singular Value
    Decomposition.

    :param X:n x m dimensional matrix which n indicate the number of the dimension and m indicate the number of points
    :return [C, N] where C is a centroid vector and N is a normal vector
    :rtype tuple
    """
    # Find the average of points (centroid) along the columns
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, V = np.linalg.svd(CX)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = V[-1]
    return C, N


def distance(x, C, N):
    """
    Calculate an orthogonal distance between the points and the standard
    Args:
    :param x: n dimensional vector (atom coordinates)
    :param C: n dimensional vector which indicate the centroid of the standard
    :param N: n dimensional vector which indicate the normal vector of the standard
    :return orthogonal distance
    """
    return np.dot(x - C, N)

if __name__ == '__main__':
    main()