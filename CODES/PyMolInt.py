import optparse
from pyrosetta.rosetta import *
from pyrosetta import *
init()
print("Secondary Protein Structure Prediction")
import os
def pose_structure(pose, display_residues = []):
    # store the pose's number of residues, example Python syntax
    nres = pose.total_residue()

    # 1. obtain the pose's sequence
    sequence = pose.sequence()

    # 2. obtain a list of PDB numbering and icode as a single string
    pdb_info = pose.pdb_info()
    PDB_nums = [(str( pdb_info.number(i)) + pdb_info.icode(i)).strip()
        for i in range(1, nres + 1)]
    # 3. obtains a list of the chains organized by residue
    chains = [pdb_info.chain(i) for i in range(1, nres + 1)]
    # 4. extracts a list of the unique chain IDs
    unique_chains = []
    for c in chains:
        if c not in unique_chains:
            unique_chains.append(c)

    # start outputting information to screen
    print('\n' + '='*80)
    print('Loaded from' , pdb_info.name())
    print(nres , 'residues')
    print(len(unique_chains), 'chain(s) ('+ str(unique_chains)[1:-1] + ')')
    print('Sequence:\n' + sequence)

    # this object is contained in PyRosetta v2.0 and above
    # 5. obtain the pose's secondary structure as predicted by PyRosetta's
    #    built-in DSSP algorithm
    DSSP = protocols.moves.DsspMover()
    DSSP.apply(pose)    # populates the pose's Pose.secstruct
    ss = pose.secstruct()
    print( 'Secondary Structure:\n' + ss )
    print( '\t' + str(100. * ss.count('H') / len(ss))[:4] + '% Helical' )
    print( '\t' + str(100. * ss.count('E') / len(ss))[:4] + '% Sheet' )
    print( '\t' + str(100. * ss.count('L') / len(ss))[:4] + '% Loop' )

    # 6. obtain the phi, psi, and omega torsion angles
    phis = [pose.phi(i) for i in range(1, nres + 1)]
    psis = [pose.psi(i) for i in range(1, nres + 1)]
    omegas = [pose.omega(i) for i in range(1, nres + 1)]

    # this object is contained in PyRosetta v2.0 and above
    # create a PyMOLMover for exporting structures directly to PyMOL
    pymover = PyMOLMover()
    pymover.apply(pose)    # export the structure to PyMOL (optional)

    # 7. output information on the requested residues
    # use a simple dictionary to make output nicer
    ss_dict = {'L':'Loop', 'H':'Helix', 'E':'Strand'}
    for i in display_residues:
        print( '='*80 )
        print( 'Pose numbered Residue', i )
        print( 'PDB numbered Residue', PDB_nums[i-1] )
        print( 'Single Letter:', sequence[i-1] )
        print( 'Chain:', chains[i-1] )
        print( 'Secondary Structure:', ss_dict[ss[i-1]] )
        print( 'Phi:', phis[i-1] )
        print( 'Psi:', psis[i-1] )
        print( 'Omega:', omegas[i-1] )
        # extract the chis
        chis = [pose.chi(j + 1, i) for j in range(pose.residue(i).nchi() )]
        for chi_no in range(len(chis)):
            print( 'Chi ' + str(chi_no + 1) + ':', chis[chi_no] )
    print( '='*80 )

parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',default = '5ggf.clean.pdb',  help = 'the PDB file containing the loop to remodel')
parser.add_option('--residues', dest = 'residues',default = '', help = 'the (pose numbered) residues to inspect carefully')
(options, args) = parser.parse_args()

pdb_filename = options.pdb_filename
pose = Pose()
# load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# default to the median residue number
residues = options.residues
if not options.residues:
    residues = [int(pose.total_residue()/2)]
elif options.residues == 'all':
    # accept the word 'all' in place of a residue list
    residues = range(1, pose.total_residue() + 1)
else:
    # please provide the residues of interest as, delimited
    residues = [int(r) for r in options.residues.split(',')]

pose_structure(pose, residues)