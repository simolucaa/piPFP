#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import csv
import numpy as np
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Permute a FASTA file.')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Path to the input FASTA file (for pi calculation) or folder (for pangenome openness).'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='Path to the output folder.'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help="Verbose mode"
    )
    parser.add_argument(
        '-w',
        '--window',
        type=int,
        help="Set the window size for PFP (min 4)"
    )
    parser.add_argument(
        '-p',
        '--modulus',
        type=int,
        help="Set the modulus for PFP (min 10)"
    )
    parser.add_argument(
        '-m',
        '--mode',
        type=str,
        default='pi',
        help="Set the mode for piPFP (pi or openness)"
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        help="Number of threads to use"
    )

    args = parser.parse_args()

    if args.window is None:
        if args.mode == 'pi':
            args.window = 10
        elif args.mode == 'openness':
            args.window = 4
    if args.modulus is None:
        if args.mode == 'pi':
            args.modulus = 100
        elif args.mode == 'openness':
            args.modulus = 10

    return args

def check_paths(input_path, output_path):
    missing_paths = []
    if not os.path.exists(input_path):
        missing_paths.append("Input path")
    if not os.path.exists(output_path):
        missing_paths.append("Output path")

    if missing_paths:
        logging.error(f" {' and '.join(missing_paths)} do not exist")
        exit(1)

def is_fasta(filename):
    try:
        with open(filename, "r") as handle:
            first_line = handle.readline().strip()
            return first_line.startswith(">")
    except Exception as e:
        return False

def write_tsv_pi(output_file, filename, window_size, modulus, dictionary_size, parse_size, pi):
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['Filename', 'Window Size', 'Modulus', 'Dictionary Size', 'Parse Size', 'pi'])
        writer.writerow([filename, window_size, modulus, dictionary_size, parse_size, pi])

def calculate_pi(args):
    logging.info(f"Calculating pi value with w = {args.window} and p = {args.modulus}")
    logging.info(f"Input file: {os.path.abspath(args.input)}")

    filename = os.path.basename(args.input)
    output_file = os.path.join(args.output, f"{filename}_pi.tsv")

    # Handle single-thread or multi-thread
    if args.threads:
        logging.info(f"Using {args.threads} threads")
        cmd = f"./piPFP.x -w {args.window} -p {args.modulus} -t {args.threads} {args.input}" # Update file.x 
        logging.info(f"Running command: {cmd}")
    else:
        logging.info("Running in single-thread mode")
        cmd = f"./piPFP_NT.x -w {args.window} -p {args.modulus} {args.input}" # Update file.x 
        logging.info(f"Running command: {cmd}")

    # Run command
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Parse output
    dict_size = int(result.stdout.split("Sum of lenghts of dictionary words: ")[1].split("\n")[0])
    parse_size = int(result.stdout.split("Total number of words: ")[1].split("\n")[0]) * 4
    pi = dict_size + parse_size
    print(f"\npi: {pi}\n")

    write_tsv_pi(output_file, filename, args.window, args.modulus, dict_size, parse_size, pi)
    logging.info(f"File saved at {output_file}")

def heaps_law (x, y):
    log_x = np.log(x)
    log_y = np.log(y)
    b, log_a = np.polyfit(log_x, log_y, 1)
    a = np.exp(log_a)
    return [a, b]

def load_data(filename):
    with open(filename, "r") as f:
        return [float(x) for x in f.read().split()]
    
def write_tsv_alpha(output_file, filename, window_size, modulus, alpha):
    with open(f"{output_file}_alpha.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['Filename', 'Window Size', 'Modulus', 'alpha'])
        writer.writerow([filename, window_size, modulus, alpha])

def calculate_alpha(args): 
    if os.path.isfile(args.input):
        logging.error("Please provide a folder instead of a file")
        exit(1)
    
    non_fasta_files = []
    for filename in os.listdir(args.input):
        filepath = os.path.join(args.input, filename)
        if os.path.isdir(filepath):
            continue
        if not is_fasta(filepath):
            non_fasta_files.append(filename)

    if non_fasta_files:
        logging.error(f"The following files are not FASTA files:\n{', '.join(non_fasta_files)}\nPlease provide a folder containing only FASTA files.")
        exit(1)

    logging.info(f"Calculating alpha value with w = {args.window} and p = {args.modulus}")
    logging.info(f"Input folder: {os.path.abspath(args.input)}")

    foldername = os.path.basename(args.input)
    output_file = os.path.join(args.output, foldername)

    # Handle single-thread or multi-thread
    if args.threads:
        logging.info(f"Using {args.threads} threads")
        cmd = f"./piPFP_growth.x  -w {args.window} -p {args.modulus} -t {args.threads} -o {output_file} {args.input}" # Update file.x 
        logging.info(f"Running command: {cmd}")
    else:
        logging.info("Running in single-thread mode")
        cmd = f"./piPFP_growth_NT.x -w {args.window} -p {args.modulus} -o {output_file} {args.input}" # Update file.x (single-thread NT not working?)
        logging.info(f"Running command: {cmd}")

    # Run command
    subprocess.run(cmd, shell=True, capture_output=True, text=True)

    logging.info(f"Growth and hist files saved at {output_file}.growth and {output_file}.hist")

    # Load data
    y = load_data(f"{output_file}.growth")
    y_new = y.copy()
    for i in range(1, len(y)):
        y_new[i] = y[i] - y[i - 1]
    y_new = y_new[1:]
    x = np.arange(2, len(y) + 1)
    y = y_new

    # Fit and get alpha
    [a, b] = heaps_law(x, y)
    print(f"\n")
    print(rf"ftot = ({a:.2e}) * N^({b:.2f})")
    print(f"alpha: {abs(b):.2f}\n")

    write_tsv_alpha(output_file, foldername, args.window, args.modulus, abs(b))
    logging.info(f"File saved at {output_file}_alpha.tsv")

###--- MAIN ---###
if __name__ == "__main__":
    args = parse_args()

    # Check paths
    check_paths(args.input, args.output)

    # Activate verbose mode
    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    # Check params
    if args.window < 4 and args.modulus < 10:
        logging.error("Window size must be at least 4 and modulus must be at least 10")
        exit(1)
    if args.window < 4:
        logging.error("Window size must be at least 4")
        exit(1)
    if args.modulus < 10:
        logging.error("Modulus must be at least 10")
        exit(1)

    # Select mode
    logging.info("Verbose mode activated\n")

    if args.mode == 'pi':
        calculate_pi(args)
    elif args.mode == 'openness':
        calculate_alpha(args)
    else:
        logging.error("Invalid mode. Please select 'pi' or 'openness'")
        exit(1)
