from ANNgeneratematrixes import *
from ANNgenerateresults import *

def vcf_load_references(paths):
    paths = vars(paths)
    sample = paths['sample']
    truth = paths['truth']
    try :
        orig_stdout = sys.stdout
        f = file(paths['output'] + '/interesting_set_of_vcf_scores.txt', 'w')
        sys.stdout = f
        output = paths['output']
    except :
        output = 0
        print "working in normal mode"
    return sample, truth, output

def execute_main(paths):
    sample, truth, output= vcf_load_references(paths)
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
    fill_negative_samples(array_of_predicted, array_of_truth, dict_of_samples, list_of_truth, final_array_of_samples)
    in_ANN_not_in_con_list = []
    in_con_not_in_ANN_list = []
    new_truth_dictionary = create_sample_dictionary(truth)
    sample_dictionary.update(new_truth_dictionary)
    for i in range(len(array_of_predicted)):
        if array_of_predicted[i] == 1 and array_of_truth[i] == 0 :
            if type(final_array_of_samples[i][3]) == list:
                final_array_of_samples[i] = list(final_array_of_samples[i])
                final_array_of_samples[i][3] = tuple(final_array_of_samples[i][3])
                final_array_of_samples[i] = tuple(final_array_of_samples[i])
            in_ANN_not_in_con_list.append(sample_dictionary[tuple(final_array_of_samples[i])])
        if array_of_predicted[i] == 0 and array_of_truth[i] == 1 :
            if type(final_array_of_samples[i][3]) == list:
                final_array_of_samples[i] = list(final_array_of_samples[i])
                final_array_of_samples[i][3] = tuple(final_array_of_samples[i][3])
                final_array_of_samples[i] = tuple(final_array_of_samples[i])
            in_con_not_in_ANN_list.append(sample_dictionary[tuple(final_array_of_samples[i])])
    print "in ANN but not in con list are " , in_ANN_not_in_con_list
    print "in con but not in ANN list are ", in_con_not_in_ANN_list
    print len(list_of_truth)
    print "array of predicted new length is", len(array_of_predicted)
    perf_measure(array_of_truth, array_of_predicted)
    print "final precision score is :", precision_score(array_of_truth, array_of_predicted)
    print "final recall score is :", recall_score(array_of_truth, array_of_predicted)
    print "final F1 score is : ", f1_score(array_of_truth, array_of_predicted)
    if output :
        ann_name = output + "/in_ann_not_in_con.vcf"
        vcf_reader = vcf.Reader(filename=original_vcf_reader)
        vcf_writer = vcf.Writer(open(ann_name, 'w'), vcf_reader)
        for record in in_ANN_not_in_con_list:
            vcf_writer.write_record(record)
        con_name = output + "/in_con_not_in_ann.vcf"
        vcf_reader = vcf.Reader(filename=original_vcf_reader)
        vcf_writer = vcf.Writer(open(con_name, 'w'), vcf_reader)
        for record in in_con_not_in_ANN_list:
            vcf_writer.write_record(record)

def fill_negative_samples(array_of_predicted, array_of_truth, dict_of_samples, list_of_truth, final_array_of_samples):
    for item in list_of_truth:
        fillnegative(item, dict_of_samples, array_of_predicted, array_of_truth, final_array_of_samples)
    if len(final_array_of_samples) !=  len(array_of_predicted) or len(final_array_of_samples) !=  len(array_of_truth):
        raise Exception("length is inconsistent")

def fillnegative(tuple1, sampledict, arrayofsamples, arrayoftruths, samplelist):
    tuple2 = (tuple1[0], tuple1[1], tuple1[2])
    if tuple2 in sampledict:
        for ALT in tuple1[3]:
            if ALT in sampledict[tuple2]:
                return
    arrayofsamples.append(0)
    arrayoftruths.append(1)
    samplelist.append(tuple1)

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
    sample_dictionary = {}
    opened_vcf_file_sample = vcf.Reader(open(sample, 'r'))
    for record in opened_vcf_file_sample:
        if "GL" in record.CHROM:
            continue
        sample_name = get_sample_name_from_record(record)
        sample_dictionary[sample_name] = record
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