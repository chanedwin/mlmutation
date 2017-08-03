import argparse
import os
import sys
import time
import unittest


import numpy as np
from Bio import SeqIO
from sklearn.metrics import *
from ANNgeneratematrixes import *

from methods import *

original_vcf_reader = "/data/backup/metacaller/stage/data/version6.3a/hc.vcf.normalisedtrain.vcf"


# OUTPUT GOAL : This method analyses the concordance score for the respective input sample


# This method take in the user input, and loads it into local variables. It executes the main_analyse_samples file
# It then saves files into a directory determined by the final variables


def load_and_save_data(user_input):
    user_input = vars(user_input)
    inputpaths, concordnumber, output = load_references(user_input)
    my_x_dataset, my_y_dataset, vcf_dictionary = main_analyse_samples_and_truth(inputpaths, concordnumber, output)
    save_files(my_x_dataset, my_y_dataset, output, concordnumber, vcf_dictionary)


# This method is the main method. It first prepares a dictionary of truth to be checked against. It then initialises
# a dictionary of samples with all the keys, each key being a variant call, and then fills it up each key with data
# subsequently, it removes dictionary entries that are the wrong size, and then checks whether
# each entry in the dictionary is actually true or not.
# subsequently it performs array balancing, and returns np.array datastructures of samples, as well as the dictionary
# of truth
# and list of called samples


def main_analyse_samples_and_truth(inputpaths, concordnumber, output):
    os.chdir(inputpaths)
    name3 = output + "/concordance/" + concordnumber + "/finalscores.txt"
    orig_stdout = sys.stdout
    f = file(name3 + '.txt', 'w')
    sys.stdout = f
    truthdict = generate_truth_list(inputpaths)
    full_dictionary, vcf_dictionary = get_dictionary_keys(inputpaths)
    full_dictionary = fill_sample_dictionary(full_dictionary, inputpaths)
    sample_dictionary, list_of_called_samples = filter_dictionary_based_on_concordance(full_dictionary, concordnumber)
    cleaned_sample_array, clean_truth_array, vcf_records_object = check_predicted_with_truth(list_of_called_samples,
                                                                                             truthdict, vcf_dictionary)
    cleaned_sample_array, clean_truth_array = fill_arrays_with_false_negatives(truthdict, sample_dictionary,
                                                                               clean_truth_array, cleaned_sample_array)
    calculate_results(cleaned_sample_array, clean_truth_array, concordnumber, output)
    return cleaned_sample_array, clean_truth_array, vcf_records_object


def calculate_results(scikittestlist, finaltruthlist, concordnumber, output_location):
    false_positive_final, true_negative_final = perf_measure(finaltruthlist, scikittestlist)
    print "final false positive is :", false_positive_final
    print "final true negative is :", true_negative_final
    print "final precision score is :", precision_score(finaltruthlist, scikittestlist)
    print "final recall score is :", recall_score(finaltruthlist, scikittestlist)
    print "final F1 score is : ", f1_score(finaltruthlist, scikittestlist)

def generate_list_of_truth(dict_of_truth):
    list_of_truth = []
    for key in dict_of_truth:
        mytuple = dict_of_truth[key]
        temptuple =[]
        for item in mytuple:
            temptuple.append(item)
        list_of_truth.append([key[0], key[1], key[2], temptuple])
    return list_of_truth

def fill_arrays_with_false_negatives(truthdict, sampledictionary, finaltruthlist, scikittestlist):
    list_of_truth = generate_list_of_truth(truthdict)
    for key in list_of_truth:
        check_truth_value_with_positive_called_samples(key, sampledictionary, finaltruthlist, scikittestlist)
    return scikittestlist, finaltruthlist


# This method goes through all the training variant calling files and extracts unique calls into a sample dictionary

def get_dictionary_keys(path):
    sample_dictionary = {}
    vcf_dictionary = {}
    for vcf_file in os.listdir(path):
        if ignore_file(vcf_file):
            continue
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        vcf_reader_2 = vcf.Reader(open(vcf_file, 'r'))
        sample_dictionary = create_dictionary_keys(vcf_reader, sample_dictionary)
        vcf_dictionary = create_vcf_dictionary(vcf_reader_2, vcf_dictionary)
    return sample_dictionary, vcf_dictionary


def create_vcf_dictionary(vcf_reader, vcf_dictionary):
    for record in vcf_reader:
        if "GL" in record.CHROM:
            continue
        sample_name = get_sample_name_from_record(record)
        vcf_dictionary[sample_name] = record  # fullname has become a key in fulldictionary
    return vcf_dictionary


# This method goes through all the training variant calling files and fills each entry in a sample dictionary
# with data. If it is empty, it returns an array of length n, where n is the number of variables
# that same caller would have provided.
# Each caller has a different amount of variables because it contains different datasets

def fill_sample_dictionary(sample_dictionary, path):
    index = 0
    for vcf_file in os.listdir(path):
        if ignore_file(vcf_file):
            continue
        index += 1
        opened_vcf_file = vcf.Reader(open(vcf_file, 'r'))
        iterate_over_file_to_extract_data(sample_dictionary, opened_vcf_file)
        refill_dictionary_with_zero_arrays_for_each_file(sample_dictionary, index)
    return sample_dictionary


# this method fills the dictionary with empty arrays with the same length as the ones that were supposed to be added

def refill_dictionary_with_zero_arrays_for_each_file(full_dictionary, index):
    empty_set = [0, 0]
    for item in full_dictionary:
        checksum = len(full_dictionary[item])
        if checksum < index:
            full_dictionary[item].append(empty_set)


# this method iterates through all the files to extract data from each sample. It uses methods from the
# methods.py function, which parses each record for data.

def iterate_over_file_to_extract_data(sample_dictionary, vcf_reader1):
    for record in vcf_reader1:
        if "GL" in record.CHROM:
            continue
        sample_name = get_sample_name_from_record(record)
        sample_data = [1, 0]
        sample_dictionary[sample_name].append(sample_data)


# this method prepares

def get_reference_dictionary_for_entropy(reference_path):
    record_dictionary = SeqIO.to_dict(SeqIO.parse(reference_path, "fasta"), key_function=get_chr)
    return record_dictionary


def ignore_file(vcf_file):
    if "vcf" not in vcf_file or "truth" in vcf_file or "breakseq" in vcf_file:
        return True
    return False


def create_dictionary_keys(vcf_reader, sample_dictionary):
    for record in vcf_reader:
        if "GL" in record.CHROM:
            continue
        sample_name = get_sample_name_from_record(record)
        sample_dictionary[sample_name] = []  # fullname has become a key in fulldictionary
    return sample_dictionary


def get_sample_name_from_record(record):
    templist = []
    for item in record.ALT:
        templist.append(str(item).upper())
    sample_name = (str(record.CHROM), str(record.POS), str(record.REF).upper(), tuple(templist))
    return sample_name


def check_sample_against_truth_dictionary(tuple_name, final_sample_list, final_truth_list, truth_dictionary):
    temp_tuple = (tuple_name[0], tuple_name[1], tuple_name[2])
    if temp_tuple in truth_dictionary:
        for alternate in tuple_name[3]:
            if alternate in truth_dictionary[temp_tuple]:
                final_sample_list.append(1)
                final_truth_list.append(1)
                return
    final_sample_list.append(1)
    final_truth_list.append(0)
    return


def load_references(input_paths):
    inputpaths = input_paths['input'][0]
    concordnumber = input_paths['concordance']
    output = input_paths['output']
    return inputpaths, concordnumber, output


def save_files(scikittestlist, finaltruthlist, output_location, concordance, vcf_reader_object):
    name1 = output_location + "/concordance/" + concordance + "/predicted.onlyprec"
    name2 = output_location + "/concordance/" + concordance + "/truth.onlyprec"
    name3 = output_location + "/concordance/" + concordance + "/vcfoutput.vcf"
    np.save(name1, scikittestlist)
    np.save(name2, finaltruthlist)
    vcf_reader = vcf.Reader(filename=original_vcf_reader)
    vcf_writer = vcf.Writer(open(name3, 'w'), vcf_reader)
    for record in vcf_reader_object:
        vcf_writer.write_record(record)


def check_predicted_with_truth(passed_list_of_samples, dictionary_of_truth, vcf_dictionary):
    final_array_of_samples = []
    final_truth_list = []
    vcf_records_object = []
    for item in passed_list_of_samples:
        check_sample_against_truth_dictionary(item[0], final_array_of_samples, final_truth_list,
                                                       dictionary_of_truth)
        vcf_records_object.append(vcf_dictionary[item[0]])
    return final_array_of_samples, final_truth_list, vcf_records_object


def filter_dictionary_based_on_concordance(fulldictionary, concordance):
    sampledictionary = {}
    samplelist = []
    for key in fulldictionary:
        count = 0
        for row in fulldictionary[key]:
            count += int(row[0])
        if count >= int(concordance):
            sampledictionary[(key[0], key[1], key[2])] = key[3]
            samplelist.append([key, fulldictionary[key]])
    return sampledictionary, samplelist


def add_mode_values_into_list_of_samples(full_dictionary, mode_value):
    removed = 0
    list_of_passed_samples = []
    for key in full_dictionary:
        second_count = 0
        for item in full_dictionary[key]:
            second_count += len(item)
        if second_count != mode_value:
            removed += 1
            continue
        list_of_passed_samples.append([key, full_dictionary[key]])
    return list_of_passed_samples


def perf_measure(y_actual, y_hat):
    true_positive = 0
    false_positive = 0
    false_negative = 0
    true_negative = 0

    for i in range(len(y_hat)):
        if y_actual[i] == 1 and y_hat[i] == 1:
            true_positive += 1
    for i in range(len(y_hat)):
        if y_hat[i] == 1 and y_actual[i] == 0:
            false_positive += 1
    for i in range(len(y_hat)):
        if y_actual[i] == 1 and y_hat[i] == 0:
            false_negative += 1
    for i in range(len(y_hat)):
        if y_hat[i] == 0 and y_actual[i] == 0:
            true_negative += 1

    print "true positives :", true_positive
    print "false positives :", false_positive
    print "false negatives :", false_negative
    print "true negatives :", true_negative

    true_positive = float(true_positive)
    false_positive = float(false_positive)
    false_negative = float(false_negative)
    if false_positive == 0 and true_positive == 0:
        false_positive_rate = 0
    else:
        false_positive_rate = false_positive / (false_positive + true_positive)
    if false_negative == 0 and true_positive == 0:
        true_negative_rate = 0
    else:
        true_negative_rate = false_negative / (false_negative + true_positive)

    return false_positive_rate, true_negative_rate


def check_truth_value_with_positive_called_samples(key, dictionary_of_samples, truth_list,
                                                   scikittestlist):
    tuplekey = (key[0],key[1],key[2])
    if tuplekey in dictionary_of_samples:
        for alternate in key[3]:
            if alternate in dictionary_of_samples[tuplekey]:
                return
    scikittestlist.append(0)
    truth_list.append(1)

class MyTest(unittest.TestCase):
    def test(self):
        self.assertEqual(fun(3), 4)

if __name__ == "__main__":
    start_time = time.time()
    np.seterr(divide='raise', invalid='raise')
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files", nargs='+')
    parser.add_argument('-C', '--concordance', help="")
    parser.add_argument('-o', '--output', help="")
    paths = parser.parse_args()
    load_and_save_data(paths)
