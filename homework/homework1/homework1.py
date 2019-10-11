# Original code based on Dan Knights' work
#   for CSCI 5481 at University of Minnesota
#
# Modified code by Brian Cooper
#
# Help: python homework1.py -h
import sys, os
import argparse
from subprocess import Popen, PIPE

# setup and establish command line arguments
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework1.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%prog 2.0')
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-t","--taxonomy",
                      default=None,
                      required=True,
                      help="Path to taxonomy file [required]")
    parser.add_argument("-o","--output",
                      default=None,
                      required=True,
                      help="Path to output file [required]")
    parser.add_argument("-c","--command",
                      default='./burst',
                      help="Path to BURST command")
    parser.add_argument("-V","--verbose",
                      action="store_true",
                      help="Verbose output")
    return parser

# runs BURST to search query sequences against reference sequences
def run_burst(query, ref, taxonomy, output, burst_cmd='./burst', verbose=False):
    """thread worker function"""

    # construct command parameters for BURST based on command line input
    cmd = f"{burst_cmd} -r {ref} -q {query} --taxonomy {taxonomy} -o {output}"

    return run_command(cmd, verbose=verbose)

# runs the given command and returns return value and output
def run_command(cmd, verbose=False):
    # verbose output
    if verbose:
        # print entered command
        print(f"\nCommand: {cmd}\n")

    # process command line pipe stream
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    # print return value code
    print(f"Return value: {return_value}\n")

    # print output
    print(stdout)

    # print any errors
    if stderr:
        print(f"Error(s): {stderr}\n")
    else:
        print("No errors!\n")

    return return_value, stdout, stderr

# main program execution
if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    # execute the BURST program with provided command-line arguments
    run_burst(args.query, args.ref, args.taxonomy, args.output, args.command,
              args.verbose)
