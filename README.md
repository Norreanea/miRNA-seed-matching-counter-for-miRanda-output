# miRNA-seed-matching-counter-for-miRanda-output
miRNA seed matching counter for miRanda output

**miRNA Seed Complementarity Counter** is a Python script that analyzes miRNA-mRNA interactions based on miRanda tool outputs. It specifically focuses on evaluating the complementarity within the seed regions (positions 2-8) of selected miRNAs against their target mRNAs. The script counts the number of perfectly complementary and wobble pairings in these seed regions, providing insights into the strength and specificity of miRNA targeting.

## Features

- **Seed region extraction:** Accurately extracts seed regions (positions 2-8) from both miRNA and mRNA sequences.
- **Complementarity analysis:** Differentiates between perfect and wobble pairings based on standard RNA complementarity rules.
- **Selective processing:** Users can specify a subset of transcripts for targeted analysis.
- **Detailed output:** Provides per-transcript counts of perfect and wobble complementarity in seed regions.
- **Debugging information:** Outputs debug statements to facilitate troubleshooting and verification.

## Requirements

- **Python Version:** Python 3.x
- **Dependencies:** Utilizes only standard Python libraries (`argparse`, `re`, `sys`).

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/miranda-seed-complementarity-counter.git
   cd miranda-seed-complementarity-counter
