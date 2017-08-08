import argparse
import vcf
import csv
import sys

def load_references(paths):
    paths = vars(paths)
    input = paths['input']
    output = paths['output']
    size = paths['sizes']
    return input, output, size


def add_sizes_to_dictionary(dictionary, sizes_file):
    f = open(sizes_file, 'r')
    size_list = f.readlines()
    size_list = map(lambda x : x.strip('\n'), size_list)
    size_list = map(lambda x: x[3:], size_list)
    size_list = filter(lambda x : "_" not in x and "M" not in x, size_list)
    size_list = map(lambda x : x.split('\t'), size_list)
    print size_list
    for item in size_list :
        dictionary[item[0]].append(item[1])


def execute_main(paths):
    input_vcf, output_file, sizes_file = load_references(paths)
    opened_vcf_file_sample = vcf.Reader(open(input_vcf, 'r'))
    dictionary = {}
    for i in range(1,23):
        dictionary[str(i)] = [0,0]
    dictionary['X'] = [0, 0]
    dictionary['Y'] = [0, 0]
    add_sizes_to_dictionary(dictionary, sizes_file)
    for record in opened_vcf_file_sample:
        index = record.CHROM[3:]
        if "_" in index :
            continue
        if "indel" in record.ID :
            dictionary[index][1] += 1
        elif "snp" in record.ID :
            dictionary[index][0] += 1
        else :
            raise Exception("something not found")

    with open(output_file, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for key in dictionary :
            print dictionary[key]
            snp_per_base = float(dictionary[key][0])/float(dictionary[key][2])
            indel_per_base = float(dictionary[key][1]) / float(dictionary[key][2])
            spamwriter.writerow([key]+dictionary[key]+[snp_per_base] + [indel_per_base])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files")
    parser.add_argument('-o', '--output', help="give directories with files")
    parser.add_argument('-s', '--sizes', help="give directories with files")
    paths = parser.parse_args()
    execute_main(paths)