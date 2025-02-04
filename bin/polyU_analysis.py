import pysam
from tqdm import tqdm
import operator

samfile_name = "/home/stefan/Synology/Data_nano_ribolyzer/directRNA_004/20231114_RNA004_NP_Nuc/basecalling_output/filtered.bam"
# samfile_name = "/home/stefan/Synology/Data_nano_ribolyzer/directRNA_004/20231114_RNA004_NP_Nuc/filtered_pod5/filtered_rebasecalled_aligned.bam"
samfile = pysam.AlignmentFile(samfile_name, "rb")

number_of_reads = samfile.count()
print(number_of_reads)

last_20 = [{"A": 0, "C": 0, "G": 0, "T": 0} for i in range(20)]

for read in tqdm(samfile.fetch(), total=number_of_reads):
    if read.is_forward:
        read_sequence = read.get_forward_sequence()
    else:
        read_sequence = read.query_sequence
    reversed_read_sequence = read_sequence[::-1]
    # print(read_sequence[-30:])
    # print(reversed_read_sequence[:30])
    position = 0
    first_A_met = False
    nucleotide_counter = 0
    adenin_counter = 0
    for nuc_index, nucleotide in enumerate(reversed_read_sequence):
        if first_A_met:
            if nucleotide == "A":
                nucleotide_counter += 1
                adenin_counter += 1
            else:
                nucleotide_counter += 1
            if (adenin_counter / nucleotide_counter) < 0.7:
                print(position)
                while (
                    reversed_read_sequence[position - 1] != "A"
                    or reversed_read_sequence[position - 2] != "A"
                    or reversed_read_sequence[position - 3] != "A"
                ):
                    if position == 0:
                        break
                    position -= 1
                break
        else:
            if nuc_index + 2 == len(reversed_read_sequence):
                position = 0
                break
            if (
                nucleotide == "A"
                and reversed_read_sequence[nuc_index + 1] == "A"
                and reversed_read_sequence[nuc_index + 2] == "A"
            ):
                nucleotide_counter += 1
                adenin_counter += 1
                first_A_met = True
        position += 1
    cutoff_sequence = read_sequence[-(position + 20) : -(position)]
    print(cutoff_sequence)
    long_cutoff_sequence = read_sequence[-(position + 20) :]
    print(long_cutoff_sequence)
    for index, base in enumerate(cutoff_sequence):
        last_20[index][base] += 1

string = ""
for dictionairy in last_20:
    string += max(dictionairy.items(), key=operator.itemgetter(1))[0]
    print(dictionairy)
print(string)
