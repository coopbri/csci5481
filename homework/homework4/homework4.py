# CSCI 5481 Homework 4
# University of Minnesota
# Most code by Brian Cooper, except:
# Credit 1.
# I used SciPy documentation to determine a useful smoothing technique for the
#   variability plot. Reference: https://docs.scipy.org/doc/scipy/reference/ ...
#       ... generated/scipy.interpolate.UnivariateSpline.html
# Credit 2.
# I also used a Stack Overflow answer to determine a number clustering method.
#   The method is by Raymond Hettinger: https://stackoverflow.com/questions/ ...
#       ... 14783947/grouping-clustering-numbers-in-python

# Imports (and uses in this file)
#   argparse: handle command-line arguments
#   matplotlib: visual data plotting
#   numpy: useful methods for working with numbers/data
#   random: random sampling for Question 4
#   scipy: smoothing method
import argparse
import matplotlib.pyplot as plt
import numpy as np
import random
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
    for a, t, g, c in pos.values():
        # Find most frequent base
        maxVal = max(a, t, g, c)

        # Append percentage of most frequent base
        percentages.append(maxVal / len(ids))

    # Write variabilities to file
    with open("variability.txt", "w") as f:
        for p in range(len(percentages)):
            if p < len(percentages) - 1:
                f.write(str(percentages[p]) + "\n")
            else:
                f.write(str(percentages[p]))

    # Close file
    f.close()

# ============================================================================ #
# Plot variability                                                             #
# ============================================================================ #
def plot(perc):
    # Set x data (positions)
    x = np.linspace(0, NUM_POSITIONS, NUM_POSITIONS)

    # Anonymous function to convert y values to true percentages
    convert = lambda n : n * 100

    # Set y data (values mapped to percentages; above lambda applied)
    y = list(map(convert, perc))

    # One-dimensional smoothing spline
    # A univariate smoothing technique was chosen because the only varying data
    #   is the y-axis data, and this technique is designed for one variable
    spline = UnivariateSpline(x, y)

    # Smooth data
    spline.set_smoothing_factor(20)

    # Scale figure output dimensions
    plt.figure(figsize=(30, 8))

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
    # This code is not mine; it is from Raymond Hettinger on Stack Overflow:
    #   https://stackoverflow.com/questions/14783947/grouping-clustering- ...
    #       ... numbers-in-python
    def cluster(data, maxgap):
        data.sort()
        groups = [[data[0]]]
        for x in data[1:]:
            if abs(x - groups[-1][-1]) <= maxgap:
                groups[-1].append(x)
            else:
                groups.append([x])
        return groups

    # Initialize empty lists for regions
    numbers = []
    regions = []

    # Extract and store variabilities less than 70%
    # 70% is a good threshold for cutting out highly variable data; found by
    #   testing of various thresholds between 50% and 90%
    for p in range(len(perc)):
        if(perc[p] < 0.70):
            numbers.append(p)

    # Group numbers using clustering helper function above
    clusters = cluster(numbers, 9)

    # Select and store cluster regions whose length > 25 base pairs
    for c in clusters:
        if len(c) > 25:
            regions.append((c[0], c[-1]))

    # Write variability regions to file
    with open("regions.txt", "w") as f:
        for start, end in regions:
            f.write(str(start) + "\t" + str(end) + "\n")

    # Close file
    f.close()

    return regions

# ============================================================================ #
# Plot variability regions                                                     #
# ============================================================================ #
def plot_regions(regions):
    # Plot positions on x-axis
    plt.plot(NUM_POSITIONS)

    # Fill variability regions in y-dimension
    plt.yticks(np.linspace(0, 0.01, 1, endpoint=True))

    # Plot based on start and end regions from file
    for start, end in regions:
        plt.axvspan(start, end, color="green")

    # Add title
    plt.title("Variability Regions", fontsize=24)

    # Add x-axis label
    plt.xlabel("Position Index", fontsize=14)

    # Save figure to file
    plt.savefig("regions.png")

    # Render plot to user's window manager
    plt.show()

# ============================================================================ #
# Randomly select 100 sequences for analysis                                   #
# ============================================================================ #
def subset(ids, seqs, perc):
    # Initialize dictionaries for whole 16S, region 1, and region 4
    whole = {}
    r1 = {}
    r4 = {}

    # Randomly select 100 sequences
    subset = random.sample(seqs, 100)

    # Append sequences to dictionary
    for item in subset:
        index = seqs.index(item)
        whole[ids[index]] = item

    # Find variable regions among random subset
    reg = regions(perc)

    # set variable regions 1 and 4 for analysis
    v1 = reg[0]
    v4 = reg[3]

    # Write regions based on whole 16S
    for key, val in whole.items():
        r1[key] = whole[key][v1[0]:v1[-1] + 1]
        r4[key] = whole[key][v4[0]:v4[-1] + 1]

    # Write whole 16S variability
    with open("whole.fna", "w") as fw:
        for key1, val1 in whole.items():
            fw.write(">" + str(key1) + "\n" + val1 + "\n")

    # Close whole 16S file
    fw.close()

    # Write variability region 1 sequences
    with open("r1.fna", "w") as fr1:
        for key2, val2 in r1.items():
            fr1.write(">" + str(key2) + "\n" + val2 + "\n")

    # Close variability region 1 file
    fr1.close()

    # Write variability region 4 sequences
    with open("r4.fna", "w") as fr4:
        for key3, val3 in r4.items():
            fr4.write(">" + str(key3) + "\n" + val3 + "\n")

    # Close variability region 4 file
    fr4.close()

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
    reg = regions(percentages)

    # Plot variable regions
    plot_regions(reg)

    # Randomly select 100 sequences for analysis
    subset(ids, seqs, percentages)
