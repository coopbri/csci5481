# CSCI 5481 Homework 3
# University of Minnesota
# Code by Brian Cooper
#
# This file implements the Nei-Saitou neighbor-joining algorithm for phylogeny
#   construction.

# Algorithm source (for equations, etc.):
#   https://en.wikipedia.org/wiki/Neighbor_joining

# Imports
#   argparse: handle command-line arguments
#   math:     useful math functions
import argparse
import math

# This is actually 1 more than number of 16S sequences (for code convenience)
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
    # Initialize sequence difference as none
    diff = 0

    # Set length as sequence 1's length
    length = len(s1)

    # Iterate over sequence, check base-by-base for differences
    for i in range(length):
        if s1[i] != s2[i]:
            diff += 1

    if diff == 0:
        # Equal sequences; zero difference
        return 0
    else:
        # Difference found
        return diff / length

# ============================================================================ #
# Build distance matrix containing pairwise sequence distances                 #
#   Write resultant matrix to file                                             #
# ============================================================================ #
def matrix(ids, seqs, output=False):
    # Print matrix to file with sequence identifiers
    if output:
        # Length is number of sequence identifiers + 1
        length = len(ids) + 1

        # Generate empty matrix
        matrix = [[(0) for x in range(length)] for y in range(length)]

        # Format beginning of matrix
        matrix[0][0] = ""

        # Fill matrix with sequence identifiers: top row
        for i in range(1, length):
            matrix[0][i] = ids[i-1]

        # Fill matrix with sequence identifiers: left column
        for j in range(1, length):
            matrix[j][0] = ids[j-1]

        # Fill matrix with pairwise sequence distances using `distance` method
        for m in range(1, length):
            for n in range(1, length):
                matrix[m][n] = distance(seqs[m-1], seqs[n-1])

        # Write matrix to file
        with open("distances.txt", "w") as f:
            for row in matrix:
                f.write("\t".join(map(str, row)) + "\t" + "\n")
    # Generate distance matrix without sequence identifiers
    else:
        # Length is number of sequence identifiers
        length = len(ids)

        # Generate empty matrix
        matrix = [[(0) for x in range(length)] for y in range(length)]

        # Fill matrix with pairwise sequence distances using `distance` method
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
                # Use Q-matrix equation from source (Wikipedia) to fill Q-matrix
                q[i, j] = (n-2) * dict[i, j] \
                    - sum([dict[i, node] for node in nodes]) \
                        - sum([dict[j, node] for node in nodes])
            else:
                # i == j, set slot to zero
                q[i, j] = 0

        # find minimum value and corresponding nodes (iterate over Q-matrix)
        min = iMin = jMin = 0
        for i, j in q:
            # Update minimum value if new minimum found
            if q[i, j] < min:
                min = q[i, j]
                iMin = i
                jMin = j

        # calculate distances from pair members (i, j) to new node,
        #   store in dictionary
        dist = {}

        # Use equation from source (Wikipedia) to find distances
        dist[iMin, x] = ((1/2) * dict[iMin, jMin]) + (math.pow(2 * (n-2), -1)) \
            * (sum([dict[iMin, node] for node in nodes]) \
                - sum([dict[jMin, node] for node in nodes]))

        # Use equation from source (Wikipedia) to find difference
        dist[jMin, x] = dict[iMin, jMin] - dist[iMin, x]

        # add tuples to result based on node, checking if location of minimum
        #   exceeds number of sequences
        if iMin >= NUM_SEQS:
            result.append((x, iMin, dist[iMin, x]))
        else:
            result.append((x, iMin + 1, dist[iMin, x]))
        if jMin >= NUM_SEQS:
            result.append((x, jMin, dist[jMin, x]))
        else:
            result.append((x, jMin + 1, dist[jMin, x]))


        # set diagonal slot distance (same sequence identifier) to 0, since
        #   these are the same sequences
        dict[x, x] = 0

        # update matrix with new distances; distance of other taxa from new node
        for node in nodes:
            if node != iMin and node != jMin:
                # Use equation from source (Wikipedia) to find distances
                dict[x, node] = (1/2) * (dict[iMin, node] \
                    + dict[jMin, node] - dict[iMin, jMin])

                dict[node, x] = dict[x, node]

        # copy dictionary to prevent "size-changing during iteration" error
        for i, j in dict.copy():
            # delete nodes that match minimum values
            if i in {iMin, jMin} or j in {iMin, jMin}:
                del dict[i, j]

        # add new node to list
        nodes.append(x)

        # Delete min nodes
        nodes = [node for node in nodes if node not in {iMin, jMin}]

        # move to next node and matrix location, using new node in place of
        #   joined neighbors and updated distances
        x -= 1
        n -= 1

    # Handle edge cases for each component of result
    for node in range(len(result)):
        if result[node][0] == NUM_SEQS + 1:
            result[node] = (NUM_SEQS, result[node][1], result[node][2])
        if result[node][0] == NUM_SEQS:
            result[node] = (NUM_SEQS + 1, result[node][1], result[node][2])

    # only two nodes: add tuple to result list
    result.append((nodes[1], nodes[0], dict[nodes[0], nodes[1]]))

    # use neighbor-joining result to create graph (stored as dictionary)
    graph = {}
    for parent, child, dist in result:
        # Handle child not in tree case
        if child not in graph:
            graph[child] = None
        # Handle parent not in tree case
        if parent not in graph:
            graph[parent] = [(child, dist)]
        # Append if in tree
        else:
            graph[parent].append((child, dist))

    return graph

# ============================================================================ #
# Determine tree edges; use preorder traversal                                 #
# ============================================================================ #
def edges(result):
    # Preorder traversal of tree
    def preorder(root, tree):
        # Initialize resulting traversal (empty)
        result = []

        # Only progress if root is not empty
        if tree[root] != None:
            for child, dist in tree[root]:
                result.append((root, child, dist))
                # Recursively traverse through tree
                result = result + preorder(child, tree)

        return result

    # Write to file (edges)
    with open("edges.txt","w") as f:
        for parent, child, dist in preorder(NUM_SEQS, result):
            f.write(str(parent) + "\t" + str(child) + "\t" + str(dist) + "\n")

# ============================================================================ #
# Newick format of edge distances; use postorder traversal                     #
# ============================================================================ #
def newick(tree, root, ids):
    # Postorder traversal of tree
    def postorder(tree, root, ids):
        # Initialize resulting traversal (empty)
        result = []

        # Handle empty root case
        if tree[root] == None:
            return ids[root - 1]

        for key, value in tree[root]:
            # Recursively traverse through tree
            result.append((postorder(tree, key, ids) + ":" + str(value)))

        # Convert result to Newick format
        return '(' + ','.join(result) + ')'

    # Write to file (Newick format of edges)
    with open('tree.txt', 'w') as f:
        f.write(postorder(tree, root, ids) + ";")

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line input
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    ids, seqs = read(args.sequence)

    # Calculate genetic distance between each pair of sequences in FASTA file,
    #   write resultant distance matrix to file `distances.txt`
    m = matrix(ids, seqs, True)

    # Generate similar matrix with no output and no sequence identifiers
    m = matrix(ids, seqs)

    # Perform Nei-Saitou neighbor-joining algorithm using distance matrix to
    #   generate phylogenetic tree
    tree = neighbor(m)

    # Write edges to file using preorder traversal
    edges(tree)

    # Write edges to file (Newick format) using postorder traversal
    newick(tree, NUM_SEQS, ids)
