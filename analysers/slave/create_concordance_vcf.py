from CreateConcordance import *

def execute_main(paths):
    concordnumber, inputpaths, output_path = load_references_create(paths)
    vcf_records_object = build_record_object(concordnumber, inputpaths)
    save_file(concordnumber, output_path, vcf_records_object)


def build_record_object(concordnumber, inputpaths):
    full_dictionary, vcf_dictionary = get_dictionary_keys(inputpaths)
    full_dictionary = fill_sample_dictionary(full_dictionary, inputpaths)
    sample_dictionary, list_of_called_samples = filter_dictionary_based_on_concordance(full_dictionary, concordnumber)
    vcf_records_object = []
    for item in list_of_called_samples:
        vcf_records_object.append(vcf_dictionary[item[0]])
    return vcf_records_object


def save_file(concordnumber, output_path, vcf_records_object):
    vcf_reader = vcf.Reader(filename=original_vcf_reader)
    name3 = output_path + "/" + concordnumber + ".vcfoutput.vcf"
    vcf_writer = vcf.Writer(open(name3, 'w'), vcf_reader)
    for record in vcf_records_object:
        vcf_writer.write_record(record)


def load_references_create(paths):
    input_paths = vars(paths)
    inputpaths = input_paths['input']
    os.chdir(inputpaths)
    concordnumber = input_paths['concordance']
    output_path = input_paths['output']
    return concordnumber, inputpaths, output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files")
    parser.add_argument('-C', '--concordance', help="")
    parser.add_argument('-o', '--output', help="")
    paths = parser.parse_args()
    execute_main(paths)