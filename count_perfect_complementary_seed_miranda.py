#!/usr/bin/env python3

import sys
import argparse
import re

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Calculate the number of complementary nucleotides in the seed sequence (positions 2-8) for selected transcripts.'
    )
    parser.add_argument('miranda_output', help='Path to miranda output file')
    parser.add_argument('selected_transcripts', help='Path to selected transcripts file')
    return parser.parse_args()

def read_selected_transcripts(selected_transcripts_file):
    """
    Reads the selected transcripts from a file and stores them in a set.

    Args:
        selected_transcripts_file (str): Path to the file containing selected transcript IDs.

    Returns:
        set: A set of selected transcript IDs.
    """
    with open(selected_transcripts_file, 'r') as f:
        # Read each line, strip whitespace, and add to the set if not empty
        selected_transcripts = set(line.strip() for line in f if line.strip())
    return selected_transcripts

def is_perfect_match(base1, base2):
    """
    Determines if two bases form a perfect complementary pair.

    Args:
        base1 (str): Nucleotide from the miRNA.
        base2 (str): Nucleotide from the target mRNA.

    Returns:
        bool: True if the pair is a perfect match, False otherwise.
    """
    # Convert bases to uppercase and replace 'T' with 'U' for RNA
    base1 = base1.upper().replace('T', 'U')
    base2 = base2.upper().replace('T', 'U')
    # Define perfect complementary pairs
    complementary_pairs = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]
    return (base1, base2) in complementary_pairs

def is_wobble_pair(base1, base2):
    """
    Determines if two bases form a wobble pair.

    Args:
        base1 (str): Nucleotide from the miRNA.
        base2 (str): Nucleotide from the target mRNA.

    Returns:
        bool: True if the pair is a wobble pair, False otherwise.
    """
    # Convert bases to uppercase and replace 'T' with 'U' for RNA
    base1 = base1.upper().replace('T', 'U')
    base2 = base2.upper().replace('T', 'U')
    # Define wobble pairs
    wobble_pairs = [('G', 'U'), ('U', 'G')]
    return (base1, base2) in wobble_pairs

def process_alignment(query_line, ref_line):
    """
    Processes an alignment block to count perfect and wobble complementary nucleotides in the seed region.

    Args:
        query_line (str): The line containing the Query sequence from miRNA.
        ref_line (str): The line containing the Ref sequence from mRNA.

    Returns:
        tuple: A tuple containing counts of perfect matches, wobble pairs, and mismatches.
    """
    # Remove the 'Query:' label and extract the sequence
    query_seq = re.sub(r'^Query:\s*', '', query_line).strip()
    # Remove leading '3'' or '5'' and trailing '3'' or '5'' characters
    query_seq = re.sub(r'^[35]\'\s*', '', query_seq)
    query_seq = re.sub(r'\s*[35]\'$', '', query_seq)
    
    # Remove the 'Ref:' label and extract the sequence
    ref_seq = re.sub(r'^Ref:\s*', '', ref_line).strip()
    # Remove leading '3'' or '5'' and trailing '3'' or '5'' characters
    ref_seq = re.sub(r'^[35]\'\s*', '', ref_seq)
    ref_seq = re.sub(r'\s*[35]\'$', '', ref_seq)

    # Reverse the miRNA sequence to align it from 5' to 3'
    query_seq = query_seq[::-1]
    # Reverse the mRNA sequence to align it from 5' to 3'
    ref_seq = ref_seq[::-1]

    # Convert sequences to lists of individual characters for easy iteration
    query_seq = list(query_seq)
    ref_seq = list(ref_seq)

    # Initialize counters for seed region positions
    miRNA_pos = 0          # Position counter for miRNA nucleotides
    seed_start = 2         # Seed region starts at position 2
    seed_end = 8           # Seed region ends at position 8

    # Initialize counters for perfect matches and wobble pairs
    perfect_matches = 0
    wobble_pairs = 0

    # Initialize indices for miRNA and mRNA sequences
    i = j = 0  # i for miRNA, j for mRNA

    # Iterate through both sequences simultaneously
    while i < len(query_seq) and j < len(ref_seq):
        c_query = query_seq[i]  # Current nucleotide in miRNA
        c_ref = ref_seq[j]      # Current nucleotide in mRNA

        # Skip spaces in miRNA sequence
        if c_query == ' ':
            i += 1
            continue

        # Skip spaces in mRNA sequence
        if c_ref == ' ':
            j += 1
            continue

        # If the current miRNA nucleotide is not a gap, increment the position counter
        if c_query != '-':
            miRNA_pos += 1

        # Check if the current position is within the seed region
        if seed_start <= miRNA_pos <= seed_end:
            # Skip if either nucleotide is a gap or space
            if c_query in ['-', ' ']:
                i += 1
                j += 1
                continue
            if c_ref in ['-', ' ']:
                i += 1
                j += 1
                continue
            # Convert bases to uppercase for consistency
            c_query_base = c_query.upper()
            c_ref_base = c_ref.upper()
            # Determine the type of complementarity and update counters accordingly
            if is_perfect_match(c_query_base, c_ref_base):
                perfect_matches += 1
            elif is_wobble_pair(c_query_base, c_ref_base):
                wobble_pairs += 1
            else:
                # Not a complementary pair; currently ignored as per requirements
                pass
        # Move to the next nucleotides in both sequences
        i += 1
        j += 1

    # Return the counts of perfect matches and wobble pairs
    return perfect_matches, wobble_pairs

def main():
    """
    The main function that orchestrates the processing of the Miranda output file
    and calculates complementary nucleotides in the seed regions for selected transcripts.
    """
    # Parse command-line arguments
    args = parse_arguments()
    # Read the set of selected transcripts
    selected_transcripts = read_selected_transcripts(args.selected_transcripts)

    # Open and read the Miranda output file
    with open(args.miranda_output, 'r') as f:
        lines = f.readlines()

    i = 0  # Initialize line counter
    # Iterate through each line in the Miranda output
    while i < len(lines):
        line = lines[i]
        # Check if the line starts with 'Read Sequence:'
        if line.startswith('Read Sequence:'):
            # Extract the transcript name using regex
            match = re.search(r'Read Sequence:(\S+)', line)
            if match:
                transcript = match.group(1)
                # Check if the transcript is among the selected transcripts
                if transcript in selected_transcripts:
                    # Extract gene name and gene length using regex
                    gene_match = re.search(r'gene=(\S+)\((\d+) nt\)', line)
                    if gene_match:
                        gene_name = gene_match.group(1)
                        gene_length = gene_match.group(2)
                    else:
                        gene_name = ''
                        gene_length = ''
                    # Move to the next line to process alignments
                    i += 1
                    # Continue processing until the next 'Read Sequence:' line or end of file
                    while i < len(lines) and not lines[i].startswith('Read Sequence:'):
                        line = lines[i]
                        # Check if the line indicates the start of an alignment block
                        if 'Forward:' in line or 'Reverse:' in line:
                            # Found an alignment block
                            # Move to the 'Query:' line
                            while i < len(lines) and not lines[i].strip().startswith('Query:'):
                                i += 1
                            if i >= len(lines):
                                break
                            query_line = lines[i].strip()
                            i += 1
                            # Skip any empty lines to reach the 'Ref:' line
                            while i < len(lines) and not lines[i].strip().startswith('Ref:'):
                                i += 1
                            if i >= len(lines):
                                break
                            ref_line = lines[i].strip()
                            i += 1
                            # Process the alignment to get perfect matches and wobble pairs
                            perfect_matches, wobble_pairs = process_alignment(query_line, ref_line)
                            # Output the results in the specified format
                            print(f"Read Sequence:{transcript} gene={gene_name}({gene_length} nt) - Complementary nucleotides in Seed: {perfect_matches} (Wobble pairings in Seed: {wobble_pairs})")
                        else:
                            # If not an alignment block, move to the next line
                            i += 1
                else:
                    # If transcript is not selected, skip to the next line
                    i += 1
            else:
                # If transcript name couldn't be extracted, skip to the next line
                i += 1
        else:
            # If the line doesn't start with 'Read Sequence:', move to the next line
            i += 1

if __name__ == '__main__':
    main()
