import optparse

from pyrosetta.rosetta import *
from pyrosetta import *

init()

def pose_scoring(pose, display_residues = []):
    pymover = PyMOLMover()
    # a. this method returns the hard-coded set of weights for the standard
    #    fullatom ScoreFunction, which is currently called "talaris2013"
    fa_scorefxn = get_fa_scorefxn()
    full_scorefxn = create_score_function('talaris2013')
    ws_patch_scorefxn = create_score_function('talaris2013', 'docking')
    patch_scorefxn = create_score_function('talaris2013')
    patch_scorefxn.apply_patch_from_file('docking')
    # e. here an empty ScoreFunction is created and the weights are set manually
    scorefxn = ScoreFunction()
    scorefxn.set_weight(core.scoring.fa_atr, 0.800)    # full-atom attractive score
    scorefxn.set_weight(core.scoring.fa_rep, 0.440)    # full-atom repulsive score
    scorefxn.set_weight(core.scoring.fa_sol, 0.750)    # full-atom solvation score
    scorefxn.set_weight(core.scoring.fa_intra_rep, 0.004)    # f.a. intraresidue rep. score
    scorefxn.set_weight(core.scoring.fa_elec, 0.700)    # full-atom electronic score
    scorefxn.set_weight(core.scoring.pro_close, 1.000)    # proline closure
    scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.170)    # short-range hbonding
    scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.170)    # long-range hbonding
    scorefxn.set_weight(core.scoring.hbond_bb_sc, 1.170)    # backbone-sidechain hbonding
    scorefxn.set_weight(core.scoring.hbond_sc, 1.100)    # sidechain-sidechain hbonding
    scorefxn.set_weight(core.scoring.dslf_fa13, 1.000)    # disulfide full-atom score
    scorefxn.set_weight(core.scoring.rama, 0.200)    # ramachandran score
    scorefxn.set_weight(core.scoring.omega, 0.500)    # omega torsion score
    scorefxn.set_weight(core.scoring.fa_dun, 0.560)    # fullatom Dunbrack rotamer score
    scorefxn.set_weight(core.scoring.p_aa_pp, 0.320)
    scorefxn.set_weight(core.scoring.ref, 1.000)    # reference identity score

    # ScoreFunction a, b, and e above have the same weights and thus return
    #    the same score for an input pose.  Likewise, c and d should return the
    #    same scores.
    # 2. output the ScoreFunction evaluations
    #ws_patch_scorefxn(pose)    # to prevent verbose output on the next line
    print( '='*80 )
    print( 'ScoreFunction a:', fa_scorefxn(pose) )
    print( 'ScoreFunction b:', full_scorefxn(pose) )
    print( 'ScoreFunction c:', ws_patch_scorefxn(pose) )
    print( 'ScoreFunction d:', patch_scorefxn(pose) )
    print( 'ScoreFunction e:', scorefxn(pose) )
    pose_score = scorefxn(pose)

    # 3. obtain the pose Energies object and all the residue total scores
    energies = pose.energies()
    residue_energies = [energies.residue_total_energy(i)
        for i in range(1, pose.total_residue() + 1)]

    # 4. obtain the non-zero weights of the ScoreFunction, active ScoreTypes
    weights = [core.scoring.ScoreType(s)
        for s in range(1, int(core.scoring.end_of_score_type_enumeration) + 1)
        if scorefxn.weights()[core.scoring.ScoreType(s)]]
    # 5. obtain all the pose energies using the weights list
    # Energies.residue_total_energies returns an EMapVector of the unweighted
    #    score values, here they are multiplied by their weights
    #    remember when performing individual investigation, these are the raw
    #    unweighted score!
    residue_weighted_energies_matrix = [
        [energies.residue_total_energies(i)[w] * scorefxn.weights()[w]
        for i in range(1, pose.total_residue() + 1)]
        for w in weights]

    # Unfortunately, hydrogen bonding scores are NOT stored in the structure
    #    returned by Energies.residue_total_energies
    # 6. hydrogen bonding information must be extracted separately
    pose_hbonds = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set( pose , False , pose_hbonds )

    # 7. create a dictionary with the pose residue numbers as keys and the
    #    residue hydrogen bonding information as values
    # hydrogen bonding information is stored as test in the form:
    # (donor residue) (donor atom) => (acceptor residue) (accecptor atom) |score
    hbond_dictionary = {}
    for residue in range(1, pose.total_residue() + 1):
        hbond_text = ''
        for hbond in range(1, pose_hbonds.nhbonds() + 1):
            hbond = pose_hbonds.hbond(hbond)
            acceptor_residue = hbond.acc_res()
            donor_residue = hbond.don_res()
            if residue == acceptor_residue or residue == donor_residue:
                hbond_text += str(donor_residue).ljust(4) + ' ' + \
                    str(pose.residue(donor_residue).atom_name(\
                        hbond.don_hatm() )).strip().ljust(4) + \
                    ' => ' + str(acceptor_residue).ljust(4) + ' ' + \
                    str(pose.residue(acceptor_residue).atom_name(\
                        hbond.acc_atm() )).strip().ljust(4) + \
                    ' |score: ' + str(hbond.energy()) + '\n'
        hbond_dictionary[residue] = hbond_text

    # 8. approximate the radius of gyration
    # there is an rg ScoreType in PyRosetta for performing this computation so
    #    a ScoreFunction can be made as an Rg calculator, likewise you can
    #    bias a structure towards more or less compact structures using this
    # NOTE: this is NOT the true radius of gyration for a protein, it uses
    #    the Residue.nbr_atom coordinates to save time, this nbr_atom is the
    #    Residue atom closest to the Residue's center of geometry
    RadG = ScoreFunction()
    RadG.set_weight(core.scoring.rg , 1)
    pose_radg = RadG(pose)

    # 9. output the pose information
	# the information is not expressed sequentially as it is produced because
	#    several PyRosetta objects and methods output intermediate information
	#    to screen, this would produce and unattractive output
    print( '='*80 )
    print( 'Loaded from' , pose.pdb_info().name() )
    print( pose.total_residue() , 'residues' )
    print( 'Radius of Gyration ~' , pose_radg )
    print( 'Total Rosetta Score:' , pose_score )
    scorefxn.show(pose)
    # this object is contained in PyRosetta v2.0 and above
    pymover.apply(pose)
    pymover.send_energy(pose)

    # 10. output information on the requested residues
    for i in display_residues:
        print( '='*80 )
        print( 'Pose numbered Residue' , i )
        print( 'Total Residue Score:' , residue_energies[i-1] )
        print( 'Score Breakdown:\n' + '-'*45 )
        # loop over the weights, extract the scores from the matrix
        for w in range(len(weights)):
            print( '\t' + core.scoring.name_from_score_type(weights[w]).ljust(20) + ':\t' ,\
                residue_weighted_energies_matrix[w][i-1] )
        print( '-'*45 )
        # print the hydrogen bond information
        print( 'Hydrogen bonds involving Residue ' + str(i) + ':' )
        print( hbond_dictionary[i][:-1] )
    print( '='*80 )

parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = 'test_in.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel')
parser.add_option('--residues', dest = 'residues',
    default = '',    # default to the median residue number
    help = 'the (pose numbered) residues to inspect carefully')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
# create a pose from the desired PDB file
# -create an empty Pose object
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# default to the median residue number
residues = options.residues
if not options.residues:
    residues = [int(pose.total_residue()/2)]
elif options.residues == 'all':
    # accept the word 'all' in place of a residue list
    residues = range(1, pose.total_residue() + 1)
else:
    # please provide the residues of interest as , delimited
    residues = [int(r) for r in options.residues.split(',')]

pose_scoring(pose, residues)