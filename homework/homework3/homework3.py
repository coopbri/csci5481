# CSCI 5481 Homework 3
# University of Minnesota
# Code by Brian Cooper
#
# This file implements the Nei-Saitou neighbor-joining algorithm for phylogeny
#   construction. Bootstrap estimation is supported.

import argparse
import numpy as np
import sys

NUM_SEQS = 62

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework3.py',
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
                    ids.append(line.replace(">","").rstrip())
                else:
                    # Sequence
                    seqs.append(line.rstrip())

        return ids, seqs

# ============================================================================ #
# Calculate genetic distance between two sequences (s1 and s2)                 #
# ============================================================================ #
def distance(s1, s2):
    diff = 0
    length = len(s1)

    for i in range(length):
        if s1[i] != s2[i]:
            diff += 1

    if diff == 0:
        return 0
    else:
        return diff / length

# ============================================================================ #
# Build distance matrix containing pairwise sequence distances                 #
#   Write resultant matrix to file                                             #
# ============================================================================ #
def matrix(ids, seqs, output):
    # Print matrix to file with sequence identifiers
    if output:
        length = len(ids) + 1

        # Generate empty matrix
        matrix = [[(0) for x in range(length)] for y in range(length)]
        matrix[0][0] = ""

        # Fill matrix with sequence identifiers
        for i in range(1, length):
            matrix[0][i] = ids[i-1]
        for j in range(1, length):
            matrix[j][0] = ids[j-1]

        # Fill matrix with pairwise sequence distances
        for m in range(1, length):
            for n in range(1, length):
                matrix[m][n] = distance(seqs[m-1], seqs[n-1])

        # Write matrix to file
        with open("distances.txt", "w") as f:
            for row in matrix:
                f.write("\t".join(map(str, row)) + "\t" + "\n")
    # Generate distance matrix without sequence identifiers
    else:
        length = len(ids)
        matrix = [[(0) for x in range(length)] for y in range(length)]
        for m in range(length):
            for n in range(length):
                matrix[m][n] = distance(seqs[m], seqs[n])

    return matrix

# ============================================================================ #
# Nei-Saitou neighbor-joining algorithm                                        #
# ============================================================================ #
def neighbor(matrix):
    # store nodes in list
    nodes = list(range(len(matrix)))

    # length of matrix
    n = len(matrix)

    # store matrix data as dictionary
    dict = {}
    for i in range(n):
        for j in range(n):
            dict[i, j] = matrix[i][j]

    # start node (at end of matrix)
    x = (NUM_SEQS - 2) * 2

    # perform algorithm while there are more distances to process
    result = []
    while n > 2:
        # build Q-matrix based on distances
        q = {}
        for i, j in dict:
            if i != j:
                q[i, j] = (n-2) * dict[i, j] \
                    - sum([dict[i, node] for node in nodes]) \
                        - sum([dict[j, node] for node in nodes])
            else:
                q[i, j] = 0

        # find minimum value and its corresponding nodes (iterate over Q-matrix)
        min = iMin = jMin = 0
        for i, j in q:
            # Update minimum value if new minimum found
            if q[i, j] < min:
                min = q[i, j]
                iMin = i
                jMin = j

        # calculate distances
        dist = {}
        dist[iMin, x] = (dict[iMin, jMin]) / 2 + (1/(2*(n-2))) \
            * (sum([dict[iMin, node] for node in nodes]) \
                - sum([dict[jMin, node] for node in nodes]))
        dist[jMin, x] = dict[iMin, jMin] - dist[iMin, x]

        # add to result based on node
        if iMin <= NUM_SEQS - 1:
            result.append((x, iMin+1, dist[iMin, x]))
        else:
            result.append((x, iMin, dist[iMin, x]))
        if jMin <= NUM_SEQS - 1:
            result.append((x, jMin+1, dist[jMin, x]))
        else:
            result.append((x, jMin, dist[jMin, x]))

        # update matrix with new distances
        for node in nodes:
            if node != iMin and node != jMin:
                dict[x, node] = 0.5 * (dict[iMin, node] \
                    + dict[jMin, node]-dict[iMin, jMin])

                dict[node, x] = dict[x, node]

        # set first slot to 0
        dict[x, x] = 0

        for i, j in dict.copy():
            # delete unnecessary values
            if i == iMin or i == jMin or j == iMin or j == jMin:
                del dict[i, j]

        # add new node to list
        nodes.append(x)

        # Delete min nodes
        nodes = [node for node in nodes if node not in {iMin, jMin}]

        # move to next node and matrix location
        x -= 1
        n -= 1

    for node in range(len(result)):
        if result[node][0] == NUM_SEQS + 1:
            result[node] = (NUM_SEQS, result[node][1], result[node][2])
        if result[node][0] == NUM_SEQS:
            result[node] = (NUM_SEQS + 1, result[node][1], result[node][2])

    # only two nodes: add to result
    result.append((nodes[1], nodes[0], dict[nodes[0], nodes[1]]))

    # convert tuple list to dictionary
    dic = {}
    for parent, child, dist in result:
        if child not in dic:
            dic[child] = None
        if parent not in dic:
            dic[parent] = [(child, dist)]
        else:
            dic[parent].append((child, dist))

    return dic

# ============================================================================ #
# Determine tree edges; use preorder traversal                                 #
# ============================================================================ #
def edges(result):
    # Preorder traversal of tree
    def preorder(root, d):
        result = []
        if d[root] != None:
            for child, dist in d[root]:
                result.append((root, child, dist))
                # Recursively traverse through tree
                result = result + preorder(child, d)
        return result

    # Write to file (edges)
    with open("edges.txt","w") as f:
        for parent, child, dist in preorder(NUM_SEQS, result):
            f.write(str(parent) + "\t" + str(child) + "\t" + str(dist) + "\n")

# ============================================================================ #
# Newick format of edge distances; use postorder traversal                     #
# ============================================================================ #
def newick(d, root, lst1):
    # Postorder traversal of tree
    def postorder(d, root, lst1):
        lst = []
        if d[root] == None:
            return lst1[root-1]
        for key, value in d[root]:
            # Recursively traverse through tree
            lst.append((postorder(d, key, lst1) + ":" + str(value)))
        result = '(' + ','.join(lst) + ')'
        return result

    # Write to file (Newick format of edges)
    with open('tree.txt', 'w') as f:
        f.write(postorder(d,root,lst1) + ";")

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line parameters
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    ids, seqs = read(args.sequence)

    # Calculate genetic distance between each pair of sequences in FASTA file
    #   Write resultant distance matrix to file `distances.txt`
    m = matrix(ids, seqs, True)

    # Generate similar matrix with no output and no sequence identifiers
    m = matrix(ids, seqs, False)

    # Perform neighbor-joining algorithm with distance matrix
    result = neighbor(m)

    # Write edges to file using preorder traversal
    edges(result)

    # Write edges to file (Newick format) using postorder traversal
    newick(result, NUM_SEQS, ids)
