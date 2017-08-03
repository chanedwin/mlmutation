from ANNgeneratematrixes import *
from ANNgenerateresults import *

def vcf_load_references(paths):
    paths = vars(paths)
    sample = paths['sample']
    truth = paths['truth']
    try :
        orig_stdout = sys.stdout
        f = file(paths['output'] + 'scores.txt', 'w')
        sys.stdout = f
    except :
        print "working in normal mode"
    return sample, truth

def execute_main(paths):
    sample, truth = vcf_load_references(paths)
    truth_dictionary = {}
    create_truth_dictionary(truth_dictionary,truth)
    print "length of truth dictionary", len(truth_dictionary)
    sample_dictionary = create_sample_dictionary(sample)
    print "length of sample dictionary", len(sample_dictionary)
    actual_predictions, final_array_of_samples, final_truth_list = match_created_with_truth(truth_dictionary,
                                                                                            sample_dictionary)
    print "lengths of actual_predictions, final_array_of_samples, final_truth_list", len(actual_predictions), len(final_array_of_samples), len(final_truth_list)
    array_of_predicted, array_of_truth, dict_of_samples, list_of_truth = restructure_data(actual_predictions,
                                                                                          final_array_of_samples,
                                                                                          final_truth_list,
                                                                                          truth_dictionary)
    fill_negative_samples(array_of_predicted, array_of_truth, dict_of_samples, list_of_truth)
    print len(list_of_truth)
    print "array of predicted new length is", len(array_of_predicted)
    perf_measure(array_of_truth, array_of_predicted)
    print "final precision score is :", precision_score(array_of_truth, array_of_predicted)
    print "final recall score is :", recall_score(array_of_truth, array_of_predicted)
    print "final F1 score is : ", f1_score(array_of_truth, array_of_predicted)


def fill_negative_samples(array_of_predicted, array_of_truth, dict_of_samples, list_of_truth):
    for item in list_of_truth:
        fillnegative(item, dict_of_samples, array_of_predicted, array_of_truth)


def restructure_data(actual_predictions, final_array_of_samples, final_truth_list, generated_truth_dictionary):
    dict_of_samples = generate_sample_dictionary_compare(final_array_of_samples)
    list_of_truth = generate_list_of_truth(generated_truth_dictionary)
    array_of_predicted = list(actual_predictions)
    array_of_truth = list(final_truth_list)
    return array_of_predicted, array_of_truth, dict_of_samples, list_of_truth


def match_created_with_truth(generated_truth_dictionary, sample_dictionary):
    final_array_of_samples = []
    actual_predictions = []
    final_truth_list = []
    for key in sample_dictionary:
        check_sample_against_truth_dictionary(key, final_truth_list, generated_truth_dictionary)
        final_array_of_samples.append(key)
        actual_predictions.append(1)
    return actual_predictions, final_array_of_samples, final_truth_list


def create_sample_dictionary(sample):
    opened_vcf_file_sample = vcf.Reader(open(sample, 'r'))
    sample_dictionary = create_dictionary_keys(opened_vcf_file_sample, {})
    for record in opened_vcf_file_sample:
        if "GL" in record.CHROM:
            continue
        sample_name = get_sample_name_from_record(record)
        if sample_name not in sample_dictionary :
            sample_dictionary[sample_name] =[1]
        else :
            sample_dictionary[sample_name].append(1)
    return sample_dictionary



def generate_sample_dictionary_compare(list_of_samples):
    dict_of_samples = {}
    for i in range(len(list_of_samples)):
        item = list_of_samples[i]
        new_key = (item[0], item[1], item[2])
        new_value = item[3]
        if new_key not in dict_of_samples:
            dict_of_samples[new_key] = new_value
        else:
            dict_of_samples[new_key] = list(dict_of_samples[new_key])
            dict_of_samples[new_key].extend(new_value)
            dict_of_samples[new_key] = tuple(dict_of_samples[new_key])
            # print dict_of_samples[new_key]
    return dict_of_samples

if __name__ == "__main__":
    np.seterr(divide='raise', invalid='raise')
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-s', '--sample', help="give directories with files")
    parser.add_argument('-t', '--truth', help="give directories with files")
    parser.add_argument('-o', '--output', help="give directories with files")
    paths = parser.parse_args()
    execute_main(paths)