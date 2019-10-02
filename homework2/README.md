### CSCI 5481 Homework 2
###### Brian Cooper

## Background Information
The file `homework2.py` implements the <a href="https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm">Needleman-Wunsch</a> global sequence alignment algorithm. Provided a proper match text file, it can implement an anchored version of the algorithm which utilizes information about where matches occur in a set of query and reference genomes. This information is provided as a text file.

## Running the Program
#### Vanilla Needleman-Wunsch
The program can be run as follows:<br>
`python homework2.py -q path_to_query_sequence.fa -r path_to_reference_sequence.fa`

#### Anchored Needleman-Wunsch
To run the anchored version, a match file needs to be provided in addition to the query and reference sequences. It can be run with:<br>
`python homework2.py -q path_to_query_sequence.fa -r path_to_reference_sequence.fa -m path_to_matches.txt`


## Random Permutation Experiment
Another option is to run an experiment that randomly permutes the given sequences and then performs alignment on these permutations. This experiment is performed 10,000 times and can be utilized by providing the `-p` or `--permute` flag when calling the program. After completion, a histogram is displayed which shows the score frequency across each of the permutations.

## Dependencies
Note that the program requires the following Python packages:
- <a href="https://numpy.org/">NumPy</a>
- <a href="https://pypi.org/project/python-string-utils/">string_utils</a>
- <a href="https://matplotlib.org/">matplotlib (pyplot module)</a>
