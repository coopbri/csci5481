# Homework 1 Question 4
#   Brian Cooper
#   CSCI 5481 at University of Minnesota

# regex library (Question 4c)
import re

if __name__ == '__main__':
    # ==========================================================================
    # Question 4a: fraction of original input query sequences at 97+ percent
    # ==========================================================================

    # read output file, line-by-line
    outfile = open("output.txt").readlines()

    # store all values above threshold (97 percent similarity)
    data = [float(n.split()[2]) for n in outfile if float(n.split()[2]) >= 97.0]

    queryMatch = len(data) / 130727
    print("Question 4a")
    print("===========")
    print(f"Sequences that match database at 97%+ similarity: \
    {round(queryMatch*100, 4)}%\n")

    # ==========================================================================
    # Question 4b: most common bacterial species
    # ==========================================================================

    # regex pattern to find species names
    pattern = "s__(.*)"

    # search each line for species name string, add to list
    species = []
    for line in outfile:
        match = re.search(pattern, line).group()
        species.append(match)

    # Remove indeterminate species from list
    determSpecies = [s for s in species if s != "s__"]

    # Determine most common species in list (highest frequency)
    common = max(set(determSpecies), key=determSpecies.count)
    print("Question 4b")
    print("===========")
    print(f"Most common determinate species:    {common}\n")

    # ==========================================================================
    # Question 4c: average percent similarity of matches
    # ==========================================================================

    avg = sum(data) / len(data)
    print("Question 4c")
    print("===========")
    print(f"Average percent similarity of matches at 97%+ similarity: \
    {round(avg, 4)}%\n")
