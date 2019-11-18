# CSCI 5481 Homework 4
# University of Minnesota
# Code by Brian Cooper
# I used SciPy documentation to determine a useful smoothing technique for the
#   variability plot. Reference: https://docs.scipy.org/doc/scipy/reference/ ...
#       ... generated/scipy.interpolate.UnivariateSpline.html

# Imports
#   argparse: handle command-line arguments
#   matplotlib: visual data plotting
#   numpy: useful methods for working with numbers
#   scipy: smoothing method
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
    # Bases at each position
    pos = {}

    # Fill postion dictionary with empty bases
    for i in range(len(seqs[0])):
        pos[i] = (0, 0, 0, 0)

    # Scan through all sequences
    for seq in seqs:
        for i in range(len(seq)):
            # check bases from sequences
            if seq[i] == "A":
                pos[i] = (pos[i][0] + 1, pos[i][1], pos[i][2], pos[i][3])
            elif seq[i] == "T":
                pos[i] = (pos[i][0], pos[i][1] + 1, pos[i][2], pos[i][3])
            elif seq[i] == "G":
                pos[i] = (pos[i][0], pos[i][1], pos[i][2] + 1, pos[i][3])
            elif seq[i] == "C":
                pos[i] = (pos[i][0], pos[i][1], pos[i][2], pos[i][3] + 1)

    # Determine variabilities based on max value
    for v1,v2,v3,v4 in pos.values():
        maxVal = max(v1,v2,v3,v4)
        percentage = maxVal * 1.0 / len(ids)
        percentages.append(percentage)

    # Write variabilities to file
    with open("variability.txt", "w") as f:
        for perc in range(len(percentages)):
            if perc < len(percentages) - 1:
                f.write(str(percentages[perc]) + "\n")
            else:
                f.write(str(percentages[perc]))

    # Close file
    f.close()

# ============================================================================ #
# Plot variability                                                             #
# ============================================================================ #
def plot(perc):
    # Set x and y data
    x = np.linspace(1, NUM_POSITIONS, NUM_POSITIONS)
    y = perc

    # One-dimensional smoothing spline
    spline = UnivariateSpline(x, y)

    # Smooth data
    spline.set_smoothing_factor(10)

    # Scale figure output dimensions
    plt.figure(figsize=(30,8))

    # Add title
    plt.title("Variability", fontsize=24)

    # Add axis labels
    plt.xlabel("Position Index", fontsize=14)
    plt.ylabel("% Sequence Identity", fontsize=14)

    # Plot data
    plt.plot(x, spline(x), color="green")

    # Save figure to file
    plt.savefig("variability.png")

    # Render plot to user's window manager
    plt.show()

# ============================================================================ #
# Find variable regions                                                        #
# ============================================================================ #
def regions(perc):
    # Clustering helper function
    def cluster(data, diff):
        groups = [[data[0]]]
        for x in data[1:]:
            if abs(x - groups[-1][-1]) <= diff:
                groups[-1].append(x)
            else:
                groups.append([x])
        return groups

    numbers = []
    regions = []

    # Store variabilities less than 0.75
    for i in range(len(perc)):
        if(perc[i] < 0.75):
            numbers.append(i)

    # Group numbers by clustering helper function
    groups = cluster(numbers, 9)
    count = 0

    # Select cluster regions whose lengths > 25 base pairs
    for i in groups:
        if len(i) > 25:
            regions.append((i[0], i[-1]))

    # Write variability regions to file
    with open("regions.txt", "w") as f:
        for fst, sec in regions:
            f.write(str(fst) + "\t" + str(sec) + "\n")

    # Close file
    f.close()

    return regions

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

    # Find variable regions
    regions(percentages)
