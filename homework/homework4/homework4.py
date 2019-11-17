# CSCI 5481 Homework 4
# University of Minnesota
# Code by Brian Cooper

# Imports
#   argparse: handle command-line arguments
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

NUM_POSITIONS = 1474
percentages = []

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework4.py',
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--sequence",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence (FASTA) file [required]")

    return parser

# ============================================================================ #
# Read in FASTA sequence file                                                  #
# ============================================================================ #
def read(file):
    # Initialize sequence identifer and sequence lists
    ids = []
    seqs = []

    with open(file) as f:
        # Parse each line individually
        for line in f:
            if line.startswith(">"):
                # Sequence identifier
                ids.append(line.replace(">", "").rstrip())
            else:
                # Sequence
                seqs.append(line.rstrip())

    return ids, seqs

# ============================================================================ #
# Calculate variability at each position in gapped alignment                   #
#   Variability is average identity or fraction of the most common base        #
# ============================================================================ #
def variability(ids, seqs):
    # a dictionary to record the bases at each position
	dic = {}
	for i in range(len(seqs[0])):
		# at first, each position has no A,T,G,or C
		dic[i] = (0,0,0,0)
	for seq in seqs:
		for i in range(len(seq)):
			# update the numbers of A,T,G,C at each positions, if it is a gap, ignore it
			if seq[i] == "A":
				n1 = dic[i][0]+1
				n2 = dic[i][1]
				n3 = dic[i][2]
				n4 = dic[i][3]
				dic[i] = (n1,n2,n3,n4)
            elif seq[i] == "C":
                b1 = dic[i][0]
                b2 = dic[i][1]
                b3 = dic[i][2]
                b4 = dic[i][3]+1
                dic[i] = (b1,b2,b3,b4)
			elif seq[i] == "T":
				m1 = dic[i][0]
				m2 = dic[i][1]+1
				m3 = dic[i][2]
				m4 = dic[i][3]
				dic[i] = (m1,m2,m3,m4)
			elif seq[i] == "G":
				a1 = dic[i][0]
				a2 = dic[i][1]
				a3 = dic[i][2]+1
				a4 = dic[i][3]
				dic[i] = (a1,a2,a3,a4)
			# else:
				# continue

	# write the variability to a new file
	for v1,v2,v3,v4 in dic.values():
		maxVal = max(v1,v2,v3,v4)
		percentage = maxVal * 1.0 / len(ids)
		percentages.append(percentage)

	# write the variability to a text file
	with open("variability.txt","w") as f:
		for j in range(len(percentages)):
			if j < len(percentages) - 1:
				f.write(str(percentages[j]) + "\n")
			else:
				f.write(str(percentages[j]))

	f.close()

def plot(perc):
	x = np.linspace(1, NUM_POSITIONS, NUM_POSITIONS)
	y = perc
	spl = UnivariateSpline(x, y)
	xs = np.linspace(1, NUM_POSITIONS, NUM_POSITIONS)
	spl.set_smoothing_factor(80)
	plt.plot(xs,spl(xs))
	plt.show()

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line input
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    ids, seqs = read(args.sequence)

    # Calculate variability
    variability(ids, seqs)

    # Plot variability
    plot(percentages)
