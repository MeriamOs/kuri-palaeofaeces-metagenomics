import sys

def filter_kraken_report(input_file, output_file, kmer_threshold=100):
    """
    Filters a Kraken2 report file, keeping rows where the fifth column 
    (distinct k-mers) is greater than a specified threshold.

    Args:
        input_file (str): Path to the input Kraken2 report file.
        output_file (str): Path to save the filtered report.
        kmer_threshold (int): Minimum distinct k-mer count to keep a row. Default is 100.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Write header lines or unclassified lines as-is
            if line.startswith(' ') or line.startswith('unclassified') or line.startswith('\t'):
                columns = line.split()
                if len(columns) >= 5:
                    try:
                        distinct_kmers = int(columns[4])  # Fifth column
                        if distinct_kmers > kmer_threshold:
                            outfile.write(line)
                    except ValueError:
                        # Skip lines where the fifth column isn't an integer
                        continue
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

# Command-line interface
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_kraken_report.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Default threshold is 100 distinct kmers
    filter_kraken_report(input_file, output_file, kmer_threshold=100)

