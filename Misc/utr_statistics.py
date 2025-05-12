import sys

def calculate_stats(utr_lengths):
    if len(utr_lengths) == 0:
        return 0, 0, 0, 0
    return len(utr_lengths), sum(utr_lengths) / len(utr_lengths), min(utr_lengths), max(utr_lengths)

three_prime_lengths = []
five_prime_lengths = []


if len(sys.argv) < 2:
    print("Usage: python script.py <path_to_gtf>")
    sys.exit(1)

# Read the GTF file
with open(sys.argv[1], 'r') as gtf_file:
    for line in gtf_file:
        # Skip comments
        if line.startswith('#'):
            continue

        # Parse the line
        fields = line.strip().split('\t')

        # Check if the line is malformed
        if len(fields) < 8:
            print("Skipping malformed line: ", line.strip())
            continue

        # Extract relevant fields
        utr_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])

        # Check for UTRs and update stats
        if utr_type == 'three_prime_UTR':
            three_prime_lengths.append(end - start + 1)
        elif utr_type == 'five_prime_UTR':
            five_prime_lengths.append(end - start + 1)

three_prime_count, three_prime_avg, three_prime_min, three_prime_max = calculate_stats(three_prime_lengths)
five_prime_count, five_prime_avg, five_prime_min, five_prime_max = calculate_stats(five_prime_lengths)

print(f"Three Prime UTR: Count = {three_prime_count}, Average Length = {three_prime_avg}, Min Length = {three_prime_min}, Max Length = {three_prime_max}")
print(f"Five Prime UTR: Count = {five_prime_count}, Average Length = {five_prime_avg}, Min Length = {five_prime_min}, Max Length = {five_prime_max}")
