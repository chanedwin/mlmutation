import argparse
import csv
import logging
import re
import sys
import time

import matplotlib
import numpy as np
import sklearn.manifold
import sklearn.preprocessing
import vcf

matplotlib.use('Agg')

from pomegranate import *
from buildcdpmethods import *
from multiprocessing import Pool

COMMON_SNP_SCORE = 0.4
UNCOMMON_SNP_SCORE = 0.6
TOTAL_NUMBER_OF_BAYESIAN_FEATURES = 5
NEURAL_NETWORK_PREDICTION_NULL_VALUE = -1
CLINVAR_NON_DELETERIOUS_SCORE = 0.2
CLINVAR_DELETERIOUS_SCORE = 0.8
LOAD_VCF_LOG = "loading annovar scores"
VCF_FILE_NOT_FOUND = "Could not find vcf file, exiting programme"
NAMES_OF_RELEVANT_SCORE_FIELDS = ['SIFT_score', 'LRT_score', 'MutationAssessor_score',
                                  'Polyphen2_HVAR_score', 'FATHMM_score']
INDEX_SIFT_SCORE = 0
INDEX_LRT_SCORE = 1
INDEX_MUTATIONASSESSOR_SCORE = 2
INDEX_POLYPHEN_SCORE = 3
INDEX_FATHMM_SCORE = 4
NAME_OF_GENE_FIELD = 'Gene.refGene'


def execute_main(paths):
    vcf_object = parse_vcf_data_from_vcf_file(paths)
    #perform_tsne(vcf_object)
    normalised_vcf_object = normalise_scores_in_list_of_scores(vcf_object)
    generate_bayesian_network_and_inference(normalised_vcf_object)
    print_data_of_full_list_of_scores(normalised_vcf_object)
    # end_time = time.time()
    # print "time taken", (end_time - start_time)


def fast_parse_vcf(vcf_file):
    """lazy vcf reader"""
    file_object = open(vcf_file, "rb")
    while True:
        data = file_object.readline()
        if not data:
            break
        yield data


def parse_vcf_data_from_vcf_file(paths):
    paths = vars(paths)  # convert paths to a dictionary object
    input = paths['input']
    full_dataset = []
    pool = Pool(processes=20)
    full_dataset.extend(pool.map(async_parse_data, fast_parse_vcf(input), 10000))
    full_dataset = list(filter(lambda x: x, full_dataset))  # throw away none returns
    return full_dataset


def async_parse_data(chunk):
    """
    SIFT_score : Deleterious (sift<=0.05), Tolerated (sift>0.05) ; Domain = [0,1]
    'LRT_score' : Deleterious (lrt<0.001), Tolerated (lrt >= 0.001) ; Domain = [0,1]
    'MutationAssessor_score' : High effect on Function (>3), Medium (3>x>1.5), Low (1.5>x>0), Neutral (0>x)
    Polyphen2_HVAR_score : Probably damaging (>=0.909), Possibly damaging (0.447<=pp2_hdiv<=0.909)
    , benign (pp2_hdiv<=0.446) ; Domain  [0,1]
    'FATHMM_score' D :  < -1.5 ,T : > -1.5
    'MutationTaster_pred' : Use Discrete Variable, "A" ("disease_causing_automatic"), "D" ("disease_causing"),
                        "N" ("polymorphism") or "P" ("polymorphism_automatic")
                        A -> 0.8
                        D -> 0.6
                        N -> 0.4
                        P -> 0.2

    """
    if re.match("^#.*", chunk):
        return None  # throw away headers
    split_data = re.split(";", chunk)
    filtered_split_data = list(filter(lambda x: "=." not in x, split_data))  # filter empty fields
    filtered_split_data = list(filter(lambda x: "=" in x, filtered_split_data))  # filter non-score fields
    dict_split_data = dict(map(lambda x: x.split("="), filtered_split_data))  # each field now a key:value
    relevant_data = []
    iterate_through_relevant_scores(relevant_data, dict_split_data, NAMES_OF_RELEVANT_SCORE_FIELDS)
    if not list(filter(lambda x: x, relevant_data)):
        return None  # throw away if no entries in relevant scores
    relevant_data.insert(0, dict_split_data[NAME_OF_GENE_FIELD])
    append_mutation_taster_score(relevant_data, dict_split_data)
    append_clinvar_scores_if_present(relevant_data, dict_split_data)
    append_snp_score_if_present(relevant_data, dict_split_data)
    append_neural_network_score_if_present(relevant_data, dict_split_data)
    return relevant_data


def append_mutation_taster_score(my_data, record):
    if 'MutationTaster_pred' in record:
        if record['MutationTaster_pred'] == 'A':
            my_data.append(0.8)
        elif record['MutationTaster_pred'] == 'D':
            my_data.append(0.6)
        elif record['MutationTaster_pred'] == 'N':
            my_data.append(0.4)
        elif record['MutationTaster_pred'] == 'P':
            my_data.append(0.2)
        else:
            my_data.append(0.2)
    else:
        my_data.append(0.2)


def iterate_through_relevant_scores(my_data, record, fields):
    for field in fields:
        add_if_record_present_else_add_zero(my_data, record, field)


def add_if_record_present_else_add_zero(my_data, record, key):
    if key in record:
        if record[key] != ".":
            my_data.append(float(record[key]))
    else:
        my_data.append(0.0)


def append_snp_score_if_present(my_data, record):
    if 'snp138' in record:
        if record['snp138'] != ".":
            snp_present = COMMON_SNP_SCORE
    else:
        snp_present = UNCOMMON_SNP_SCORE
    my_data.append(snp_present)


def append_clinvar_scores_if_present(list_of_scores, record):
    if 'clinvar_20150629' in record:
        if record['clinvar_20150629'] != ".":
            clinvar_score = CLINVAR_DELETERIOUS_SCORE
    else:
        clinvar_score = CLINVAR_NON_DELETERIOUS_SCORE
    list_of_scores.append(clinvar_score)


def append_neural_network_score_if_present(list_of_important_mutations, record):
    if 'NN_prediction' in record:
        NN_prediction = float(record['NN_prediction'])
    else:
        NN_prediction = NEURAL_NETWORK_PREDICTION_NULL_VALUE
    list_of_important_mutations.append(NN_prediction)


# MUST CHANGE USE BINARY MAPPING
def normalise_scores_in_list_of_scores(full_list_of_scores):
    bayesian_scores = list(map(lambda x: normalise_inputs(x), full_list_of_scores))
    floored_bayesian_scores = list(map(lambda x: floor_function(x), bayesian_scores))
    ceiling_bayesian_scores = list(map(lambda x: ceiling_function(x), floored_bayesian_scores))
    return ceiling_bayesian_scores


def floor_function(my_data):
    return list(map(lambda x: 0.2 if (type(x) == float and x < 0.2) else x, my_data))


def ceiling_function(my_data):
    return list(map(lambda x: 0.8 if (type(x) == float and x > 0.8) else x, my_data))


def normalise_inputs(my_data):
    """
    Normalise inputs from LRT,SIFT,MutationTaster,PolyPhen and FATHMM_score
    => Normalising strategy :
    Final Domain : [0,1]  where boundaries indicate
    0 - Neutral Mutation
    1 - Functionally Important Mutation
    Uses sigmoid function

    Input Dataset
    SIFT_score : Deleterious (sift<=0.05), Tolerated (sift>0.05)
    'LRT_score' : Deleterious (lrt<0.001), Tolerated (lrt >= 0.001)
    'MutationAssessor_score' : High(ly likely deleterious) (>3), Medium (3>x>1.5), Low (1.5>x>0), Neutral (0>x)
    Polyphen2_HVAR_score : Probably damaging (>=0.909), Possibly damaging (0.447<=pp2_hdiv<=0.909)
    , benign (pp2_hdiv<=0.446) ; Domain  [0,1]
    'FATHMM_score' D(eleterous) :  < -1.5 ,T(olerated) : > -1.5
    """
    data_without_header = my_data[1:]
    SIFT_score = data_without_header[INDEX_SIFT_SCORE]
    LRT_score = data_without_header[INDEX_LRT_SCORE]
    MutationAssessor_score = data_without_header[INDEX_MUTATIONASSESSOR_SCORE]
    PolyPhen_score = data_without_header[INDEX_POLYPHEN_SCORE]
    FATHMM_score = data_without_header[INDEX_FATHMM_SCORE]
    data_without_header[INDEX_SIFT_SCORE] = 1 - SIFT_score  # from [1,0] to [0,1]
    data_without_header[INDEX_LRT_SCORE] = 1 - LRT_score  # from [1,0] to [0,1]
    data_without_header[INDEX_MUTATIONASSESSOR_SCORE] = sigmoid(MutationAssessor_score)  # from [0,inf) to [0,1]
    data_without_header[INDEX_POLYPHEN_SCORE] = PolyPhen_score  # already [0,1]
    data_without_header[INDEX_FATHMM_SCORE] = sigmoid(-FATHMM_score)  # from (inf,-inf) to [0,1]
    return [my_data[0], ] + data_without_header


def mapped_sigmoid(input_list):
    return list(map(lambda x: sigmoid(x), input_list))


def sigmoid(x):
    return 1 / (1 + np.exp(1) ** (-x))


def perform_tsne(full_list_of_scores):
    array_of_scores = map(lambda x: x[:-3], full_list_of_scores)
    array_of_mutations = map(lambda x: x[1:], array_of_scores)
    array_of_names = map(lambda x: x[0], array_of_scores)
    processed_array_of_scores = count_and_impute_nan_values(array_of_mutations)
    results = sklearn.manifold.TSNE().fit_transform(processed_array_of_scores)
    assert len(results) == len(full_list_of_scores)
    my_list_1, my_list_2 = zip(*results)
    zipped = zip(array_of_names, my_list_1, my_list_2)
    write_to_csv(zipped, "tsne_results.csv")
    sys.exit()


def count_and_impute_nan_values(array_of_scores):
    num_nan = np.isnan(array_of_scores).sum()
    if num_nan > 0:
        logging.warning("found " + str(num_nan) + " nan values in dataset, imputing via mean strategy")
        imp = sklearn.preprocessing.Imputer()
        processed_array_of_scores = imp.fit_transform(array_of_scores)
        return processed_array_of_scores
    else:
        return array_of_scores


def generate_bayesian_network_and_inference(full_list_of_scores):
    for record in full_list_of_scores:
        ClinVar_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        PolyPhen2_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        LRT_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        MutationTaster_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        MutationAssessor_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        SIFT_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        FATHMM_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        rs_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        import_cdp = get_cdp(3, [record[-1], record[-2], 0.8])
        functional_cdp = get_cdp(7, record[1:8])
        functional_gene = ConditionalProbabilityTable(functional_cdp, [SIFT_gene, LRT_gene, MutationAssessor_gene,
                                                                       PolyPhen2_gene, FATHMM_gene, MutationTaster_gene
            , ClinVar_gene])
        real_gene = DiscreteDistribution({'True': 0.5, 'False': 0.5})
        importgene = ConditionalProbabilityTable(import_cdp, [real_gene, rs_gene, functional_gene])

        # set up states
        s1 = State(real_gene, name="Real Gene")
        s2 = State(functional_gene, name="Functional Gene")
        s3 = State(importgene, name="Important Gene")
        s4 = State(ClinVar_gene, name="ClinVar")
        s5 = State(PolyPhen2_gene, name="PolyPhen")
        s6 = State(LRT_gene, name="LRT")
        s7 = State(MutationTaster_gene, name="MutationTaster")
        s8 = State(MutationAssessor_gene, name="MutationAssessor")
        s9 = State(SIFT_gene, name="SIFT")
        s10 = State(FATHMM_gene, name="FATHMM_gene")
        s11 = State(rs_gene, name="rs_gene")

        # set up network
        network = BayesianNetwork("Gene Prediction")
        network.add_states(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11)
        network.add_edge(s1, s3)
        network.add_edge(s2, s3)
        network.add_edge(s4, s2)
        network.add_edge(s5, s2)
        network.add_edge(s6, s2)
        network.add_edge(s7, s2)
        network.add_edge(s8, s2)
        network.add_edge(s9, s2)
        network.add_edge(s10, s2)
        network.add_edge(s11, s3)
        network.bake()
        # network.plot()
        # savefig("path.png")
        # predict probability that gene is important given other probabilities
        beliefs = network.predict_proba({'Real Gene': 'True', 'ClinVar': 'True', 'PolyPhen': 'True', 'LRT': 'True',
                                         'MutationTaster': 'True', 'MutationAssessor': 'True', 'SIFT': 'True',
                                         'FATHMM_gene': 'True', 'rs_gene': 'True'})
        # print "\n".join("{}\t{}".format(state.name, belief) for state, belief in zip(network.states, beliefs))
        # get the probability that the gene is important
        prob_gene_important = beliefs[2].values()[1]
        beliefs = map(str, beliefs)
        record.append(prob_gene_important)
    full_list_of_scores.sort(key=lambda x: x[-1], reverse=True)


def print_data_of_full_list_of_scores(full_list_of_scores):
    for item in full_list_of_scores:
        print item


def write_to_csv(data_matrix, file_name):
    with open(file_name, 'wb') as csvfile:
        file_writer = csv.writer(csvfile, delimiter=' ',
                                 quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for line in data_matrix:
            file_writer.writerow(line)


"""
DEPRECATED METHODS
"""


# Deprecated
def get_scores_of_each_record(record):
    raw_list_of_important_mutations = [record.INFO['SIFT_score'], record.INFO['LRT_score'],
                                       record.INFO['MutationAssessor_score'],
                                       record.INFO['MutationTaster_score'], record.INFO['Polyphen2_HVAR_score'],
                                       record.INFO['FATHMM_score']]
    list_of_important_mutations = map(lambda x: 0.0 if x[0] == None else float(x[0]), raw_list_of_important_mutations)
    if not list(filter(lambda x: x, list_of_important_mutations)):  # if it is an empty list, throw it away
        return None
    list_of_important_mutations = map(lambda x: 0.0 if x == None else float(x), list_of_important_mutations)
    record_name = record.INFO['Gene.refGene'][0] if len(record.INFO['Gene.refGene']) == 1 else str(
        reduce(lambda x, y: x + "," + y, record.INFO['Gene.refGene']))
    list_of_important_mutations.insert(0, record_name)
    append_clinvar_scores(list_of_important_mutations, record)  # TO-FIX make append and generate the same
    append_snp_score(list_of_important_mutations, record)
    append_neural_network_score(list_of_important_mutations, record)
    return list_of_important_mutations


# Deprecated
def append_neural_network_score(list_of_important_mutations, record):
    if 'NN_prediction' in record.INFO:
        NN_prediction = float(record.INFO['NN_prediction'][0])
    else:
        NN_prediction = NEURAL_NETWORK_PREDICTION_NULL_VALUE
    list_of_important_mutations.append(NN_prediction)


# Deprecated
def append_snp_score(list_of_scores, record):
    if record.INFO['snp138'][0] != None:
        snp_present = COMMON_SNP_SCORE
    else:
        snp_present = UNCOMMON_SNP_SCORE
    list_of_scores.append(snp_present)


# Deprecated
def append_clinvar_scores(list_of_scores, record):
    if record.INFO['clinvar_20150629'][0] != None:
        clinvar_score = CLINVAR_DELETERIOUS_SCORE
    else:
        clinvar_score = CLINVAR_NON_DELETERIOUS_SCORE
    list_of_scores.append(clinvar_score)


# Deprecated
def load_vcf_object_from_input_path(paths):
    logging.info('Loading VCF from file')
    paths = vars(paths)  # convert paths to a dictionary object
    input = paths['input']
    logging.info(LOAD_VCF_LOG)
    try:
        opened_vcf_file = vcf.Reader(open(input, 'r'))
    except:
        logging.info(VCF_FILE_NOT_FOUND)
        exit()
    logging.info('Done loading VCF from file')
    new_list_of_scores = []
    for record in opened_vcf_file:
        new_list_of_scores.append(get_scores_of_each_record(record))
    new_list_of_scores = list(filter(lambda x: x, new_list_of_scores))
    return new_list_of_scores


# Deprecated
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


# Deprecated
def obtain_full_list_of_scores(vcf_object):
    logging.info('Extracting features from each vcf record')
    pool = Pool(processes=100)
    print "starting pool processes"
    start = time.time()
    full_list_of_scores = []
    num = 0
    for chunk in pool.imap_unordered(get_scores_of_each_record, vcf_object, 100):
        full_list_of_scores.extend(chunk)
        num += 1
        print "finished", num, "chunk"
    full_list_of_scores = filter(lambda x: x, full_list_of_scores)
    end = time.time()
    print "pooled time is", (end - start)
    return full_list_of_scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files")
    paths = parser.parse_args()
    logging.info('Starting Prediction Pipeline')
    execute_main(paths)
