#!/usr/bin/env python3

import sys
import argparse
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate the number of complementary nucleotides in the seed sequence (positions 2-8) for selected transcripts.')
    parser.add_argument('miranda_output', help='Path to miranda output file')
    parser.add_argument('selected_transcripts', help='Path to selected transcripts file')
    return parser.parse_args()

def read_selected_transcripts(selected_transcripts_file):
    with open(selected_transcripts_file, 'r') as f:
        selected_transcripts = set(line.strip() for line in f if line.strip())
    return selected_transcripts

def is_perfect_match(base1, base2):
    base1 = base1.upper().replace('T', 'U')
    base2 = base2.upper().replace('T', 'U')
    complementary_pairs = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]
    return (base1, base2) in complementary_pairs

def is_wobble_pair(base1, base2):
    base1 = base1.upper().replace('T', 'U')
    base2 = base2.upper().replace('T', 'U')
    wobble_pairs = [('G', 'U'), ('U', 'G')]
    return (base1, base2) in wobble_pairs

def process_alignment(query_line, ref_line):
    # Remove labels and extract sequences
    query_seq = re.sub(r'^Query:\s*', '', query_line).strip()
    query_seq = re.sub(r'^[35]\'\s*', '', query_seq)
    query_seq = re.sub(r'\s*[35]\'$', '', query_seq)
    ref_seq = re.sub(r'^Ref:\s*', '', ref_line).strip()
    ref_seq = re.sub(r'^[35]\'\s*', '', ref_seq)
    ref_seq = re.sub(r'\s*[35]\'$', '', ref_seq)

    # Reverse the miRNA sequence to get 5' to 3'
    query_seq = query_seq[::-1]
    ref_seq = ref_seq[::-1]

    # Convert sequences to lists of characters
    query_seq = list(query_seq)
    ref_seq = list(ref_seq)

    # Initialize positions
    miRNA_pos = 0
    seed_start = 2
    seed_end = 8

    perfect_matches = 0
    wobble_pairs = 0

    i = j = 0  # i for miRNA, j for mRNA

    while i < len(query_seq) and j < len(ref_seq):
        c_query = query_seq[i]
        c_ref = ref_seq[j]

        # Skip spaces in miRNA
        if c_query == ' ':
            i += 1
            continue

        # Skip spaces in mRNA
        if c_ref == ' ':
            j += 1
            continue

        # If c_query is not a gap, increment miRNA_pos
        if c_query != '-':
            miRNA_pos += 1

        # Process positions 2-8 in miRNA
        if seed_start <= miRNA_pos <= seed_end:
            if c_query in ['-', ' ']:
                i += 1
                j += 1
                continue
            if c_ref in ['-', ' ']:
                i += 1
                j += 1
                continue
            c_query_base = c_query.upper()
            c_ref_base = c_ref.upper()
            if is_perfect_match(c_query_base, c_ref_base):
                perfect_matches += 1
            elif is_wobble_pair(c_query_base, c_ref_base):
                wobble_pairs += 1
            else:
                # Not a complementary pair
                pass  # We can ignore mismatches for this calculation
        # Move to next bases
        i += 1
        j += 1

    return perfect_matches, wobble_pairs

def main():
    args = parse_arguments()
    selected_transcripts = read_selected_transcripts(args.selected_transcripts)

    with open(args.miranda_output, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('Read Sequence:'):
            # Extract transcript name
            match = re.search(r'Read Sequence:(\S+)', line)
            if match:
                transcript = match.group(1)
                # Check if transcript is in selected transcripts
                if transcript in selected_transcripts:
                    # Extract gene name and length
                    gene_match = re.search(r'gene=(\S+)\((\d+) nt\)', line)
                    if gene_match:
                        gene_name = gene_match.group(1)
                        gene_length = gene_match.group(2)
                    else:
                        gene_name = ''
                        gene_length = ''
                    # Process alignments for this transcript
                    i += 1
                    while i < len(lines) and not lines[i].startswith('Read Sequence:'):
                        line = lines[i]
                        if 'Forward:' in line or 'Reverse:' in line:
                            # Found an alignment block
                            # Skip until we find 'Query:' line
                            while i < len(lines) and not lines[i].strip().startswith('Query:'):
                                i += 1
                            if i >= len(lines):
                                break
                            query_line = lines[i].strip()
                            i += 1
                            # Skip any lines until we find 'Ref:' line
                            while i < len(lines) and not lines[i].strip().startswith('Ref:'):
                                i += 1
                            if i >= len(lines):
                                break
                            ref_line = lines[i].strip()
                            i += 1
                            # Process the alignment
                            perfect_matches, wobble_pairs = process_alignment(query_line, ref_line)
                            # Output the results
                            print(f"Read Sequence:{transcript} gene={gene_name}({gene_length} nt) - Complementary nucleotides in Seed: {perfect_matches} (Wobble pairings in Seed: {wobble_pairs})")
                        else:
                            i += 1
                else:
                    # Skip to the next block
                    i += 1
            else:
                i += 1
        else:
            i += 1

if __name__ == '__main__':
    main()
