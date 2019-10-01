# CSCI 5481 Homework 2
# University of Minnesota
# Code by Brian Cooper
#
# This file implements a modified version (anchored) of the Needleman-Wunsch
#   algorithm for sequence alignment.

import numpy as np
import sys
import argparse

# "constant" score values (not really constant, but functionally)
MATCH = 1
MISMATCH = -3
GAP = -2

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework2.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query sequence (FASTA) [required]")
    parser.add_argument("-r","--reference",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference sequence (FASTA) [required]")
    parser.add_argument("-m","--match",
                      default=None,
                      required=False,
                      help="Path to matches file [optional]")

    return parser

# ============================================================================ #
# Read in FASTA sequence file or match file                                    #
# ============================================================================ #
def read(file, match=False):
    # Match file (for anchored NW)
    if match:
        # Load file (with NumPy for convenience)
        m = np.loadtxt(file)
        # Human sequence indices (first two columns)
        human = m[:, 0:2] - 1
        # Fly sequence indicies (second two columns)
        fly = m[:, 2:4] - 1
        return human, fly

    # FASTA sequence file
    else:
        parsed = ""
        with open(file) as f:
            next(f).rstrip()
            for line in f:
                parsed += line.rstrip()
        return parsed

# ============================================================================ #
# Create similarity matrix for use by Needleman-Wunsch function                #
# ============================================================================ #
def create_matrix(q, r):
    # Initialize similarity matrix (zero matrix)
    V = np.zeros((len(q)+1, len(r)+1))

    # Fill row 1 and column 1 with gap penalty values
    for i in range(len(q)+1):
        V[i][0] = GAP*i
    for j in range(len(r)+1):
        V[0][j] = GAP*j

    # Fill matrix with scores
    for i in range(1, len(q)+1):
        for j in range(1, len(r)+1):
            # Determine match or mismatch
            if (q[i-1] == r[j-1]):
                diag = V[i-1][j-1] + MATCH
            else:
                diag = V[i-1][j-1] + MISMATCH

            left = V[i][j-1] + GAP
            above = V[i-1][j] + GAP

            V[i][j] = max(diag, left, above)

    return V

# ============================================================================ #
# Needleman-Wunsch algorithm                                                   #
#   (param) V: similarity matrix                                               #
#   (param) q: query sequence                                                  #
#   (param) r: reference sequence                                              #
# ============================================================================ #
def needleman_wunsch(V, q, r):
    # Initialize alignments as empty strings
    alignQ = ""
    alignR = ""

    # Length of sequences
    i = len(q)
    j = len(r)

    # Process both sequences entirely
    while i > 0 or j > 0:
        # Determine match or mismatch score
        diff = MATCH if q[i-1] == r[j-1] else MISMATCH

        # Both sequences not fully processed
        if (i > 0 and j > 0):
            # Match or mismatch found
            if V[i-1][j-1] == (V[i][j] - diff):
                alignQ = q[i-1] + alignQ
                alignR = r[j-1] + alignR
                i -= 1
                j -= 1
            # From left
            elif V[i-1][j] == (V[i][j] - GAP):
                alignQ = q[i-1] + alignQ
                alignR = "-" + alignR
                i -= 1
            # From above
            elif V[i][j-1] == (V[i][j] - GAP):
                alignQ = "-" + alignQ
                alignR = r[j-1] + alignR
                j -= 1
        # Query sequence not fully processed
        elif i > 0:
            alignQ = q[i-1] + alignQ
            alignR = "-" + alignR
            i -= 1
        # Reference sequence not fully processed
        elif j > 0:
            alignQ = "-" + alignQ
            alignR = r[j-1] + alignR
            j -= 1

    # Final score, at bottom right of matrix
    finalScore = V[-1][-1]

    return alignQ, alignR, finalScore

# ============================================================================ #
# Anchored Needleman-Wunsch algorithm                                          #
#   (param) q: query sequence                                                  #
#   (param) r: reference sequence                                              #
#   (param) qMatch: query match indices                                        #
#   (param) rMatch: reference match indicies                                   #
# ============================================================================ #
def anchored_nw(q, r, qMatch, rMatch):
    finalScore = 0.0

    new_human = q[:int(qMatch[0,0])]
    new_fly = r[:int(rMatch[0,0])]
    num_matches = np.shape(qMatch)[0]

    for n in range(num_matches):
        qStart = int(qMatch[n,0])
        rStart = int(rMatch[n,0])
        qEnd = int(qMatch[n,1])
        rEnd = int(rMatch[n,1])

        V = create_matrix(q[qStart:qEnd], r[rStart:rEnd])

        # (tmpscore, _, _) = find_alignments(qMatch[qStart:qEnd], rMatch[rStart:rEnd])
        (_, _, tmpscore) = needleman_wunsch(V, q[qStart:qEnd], r[rStart:rEnd])


        finalScore += tmpscore
        new_human += q[qStart:qEnd]
        new_fly += r[rStart:rEnd]

        if n < num_matches-1:
            qNext = int(qMatch[n+1,0])
            rNext = int(rMatch[n+1,0])

            V = create_matrix(q[qEnd:qNext], r[rEnd:rNext])

            (tmphseq, tmpfseq, tmpscore) = needleman_wunsch(V, q[qEnd:qNext], r[rEnd:rNext])

            finalScore += tmpscore
            new_human += tmphseq
            new_fly += tmpfseq

        elif n == num_matches-1:
            new_human += q[qEnd:]
            new_fly += r[rEnd:]

    return (new_human, new_fly, finalScore)

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line parameters
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence files (query and reference)
    q = read(args.query)
    r = read(args.reference)

    # Match file provided
    if args.match:
        # Process provided match file
        qMatch, rMatch = read(args.match, match=True)
        # Perform Needleman-Wunsch algorithm (anchored)
        alignQ, alignR, score = anchored_nw(q, r, qMatch, rMatch)
    # No match file provided
    else:
        # Create and fill similarity matrix for Needleman-Wunsch
        V = create_matrix(q, r)
        # Perform Needleman-Wunsch algorithm (standard)
        alignQ, alignR, score = needleman_wunsch(V, q, r)

    # Output alignments and score
    print("Query alignment:     " + alignQ)
    print("Reference alignment: " + alignR)
    print("\nFinal score: " + str(score))
