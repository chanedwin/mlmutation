import argparse
import logging
import sys
import csv
import time

import matplotlib
import vcf
import sklearn.manifold
import sklearn.preprocessing

matplotlib.use('Agg')

from pomegranate import *
from buildcdpmethods import *
from multiprocessing import Pool


#this script generates the bayesian network that ranks the mutations obtained from vcf files. Firstly, it processes
#the file by extract all the relevant features from the VCF records (stored using pyvcf). Subsequently, it
#generates a bayesian network, and updates this network using the features obtained from the VCF records.


#declare constants and scores for probability network

COMMON_SNP_SCORE = 0.4

UNCOMMON_SNP_SCORE = 0.6

TOTAL_NUMBER_OF_BAYESIAN_FEATURES = 6

NEURAL_NETWORK_PREDICTION_NULL_VALUE = -1

CLINVAR_NON_DELETERIOUS_SCORE = 0.2

CLINVAR_DELETERIOUS_SCORE = 0.8

LOAD_VCF_LOG = "loading annovar scores"
VCF_FILE_NOT_FOUND = "Could not find vcf file, exiting programme"


#main control function that loads list of vcf records, parses vcf records into a matrix of scores, performs TSNE,
#normalises the score and then generates the bayesian inference network. Finally, it prints the full list of scores

def execute_main(paths):
    vcf_object = load_vcf_object_from_input_path(paths)
    full_list_of_scores = obtain_full_list_of_scores(vcf_object)
    perform_tsne(full_list_of_scores)
    normalise_scores_in_list_of_scores(full_list_of_scores)
    generate_bayesian_network_and_inference(full_list_of_scores)
    print_data_of_full_list_of_scores(full_list_of_scores)

def load_vcf_object_from_input_path(paths):
    logging.info('Loading VCF from file')
    paths = vars(paths)  #convert paths to a dictionary object
    input = paths['input']
    logging.info(LOAD_VCF_LOG)
    try:
        opened_vcf_file = vcf.Reader(open(input, 'r'))
    except :
        logging.info(VCF_FILE_NOT_FOUND)
        exit()
    return opened_vcf_file



def obtain_full_list_of_scores(vcf_object):
    logging.info('Extracting features from each vcf record')
    pool = Pool(processes=20)
    number_of_counts = 0
    full_list_of_scores = []
    start = time.time()
    full_list_of_scores.append(pool.map(get_scores_of_each_record,vcf_object))
    end = time.time()
    print end - start
    return full_list_of_scores


def append_snp_score(list_of_scores, record):
    if record.INFO['snp138'][0] != None:
        snp_present = COMMON_SNP_SCORE
    else :
        snp_present = UNCOMMON_SNP_SCORE
    list_of_scores.append(snp_present)


def append_clinvar_scores(list_of_scores, record):
    if record.INFO['clinvar_20150629'][0] != None:
        clinvar_score = CLINVAR_DELETERIOUS_SCORE
    else:
        clinvar_score = CLINVAR_NON_DELETERIOUS_SCORE
    list_of_scores.append(clinvar_score)


    #SIFT_score [1,0]
    #LRT_score [1,0]
    #FATHMM_score[1,0]
    #MutationTaster[0,1]
    #MutationAssessor[0,1] = ???
    #Polyphen2_HVAR_score[0,1]

def get_scores_of_each_record(record):
    list_of_important_mutations = [record.INFO['SIFT_score'], record.INFO['LRT_score'], record.INFO['MutationAssessor_score'],
                          record.INFO['MutationTaster_score'], record.INFO['Polyphen2_HVAR_score'],
                          record.INFO['FATHMM_score']]
    list_of_important_mutations = map(lambda x: None if x[0] == None else float(x[0]), list_of_important_mutations)
    if not list(filter(lambda x: x != None, list_of_important_mutations)):  # if it is an empty list, throw it away
        return None
    append_clinvar_scores(list_of_important_mutations, record)  # TO-FIX make append and generate the same
    append_snp_score(list_of_important_mutations, record)
    append_neural_network_score(list_of_important_mutations, record)
    return list_of_important_mutations


def append_neural_network_score(list_of_important_mutations, record):
    if 'NN_prediction' in record.INFO:
        NN_prediction = record.INFO['NN_prediction'][0]
    else:
        NN_prediction = NEURAL_NETWORK_PREDICTION_NULL_VALUE
    list_of_important_mutations.append(NN_prediction)


def normalise_scores_in_list_of_scores(full_list_of_scores):
    for i in range(TOTAL_NUMBER_OF_BAYESIAN_FEATURES):
        min_num = 1000000
        max_num = -1000000
        for item in full_list_of_scores:
            if item[1][i] != None:
                min_num = min(min_num, item[1][i])
                max_num = max(max_num, item[1][i])
        for item in full_list_of_scores:
            if item[1][i] != None:
                value = ((item[1][i] - min_num) / (max_num - min_num) + 0.2) / 1.3
                item[1][i] = value
            else:
                item[1][i] = 0.5


def perform_tsne(full_list_of_scores):
    array_of_scores = map(lambda x : x[1][:-1], full_list_of_scores)
    print full_list_of_scores[:10]
    print array_of_scores[:10]
    imp = sklearn.preprocessing.Imputer()
    processed_array_of_scores = imp.fit_transform(array_of_scores)
    results = sklearn.manifold.TSNE().fit_transform(processed_array_of_scores)
    write_to_csv(results,"tsne_results.csv")
    sys.exit()

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
        import_cdp = get_cdp(3, [(record[0] + 0.2) / 1.3, record[3], 0.8])
        functional_cdp = get_cdp(7, record[1])
        functional_gene = ConditionalProbabilityTable(functional_cdp, [ClinVar_gene, PolyPhen2_gene, LRT_gene,
                                                                       MutationTaster_gene, MutationAssessor_gene,
                                                                       SIFT_gene, FATHMM_gene])
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
        record.append(record[2].INFO['snp138'])
        record.append(record[3])
    full_list_of_scores.sort(key=lambda x: x[4], reverse=True)


def print_data_of_full_list_of_scores(full_list_of_scores):
    for item in full_list_of_scores:
        print item[2], item, item[2].INFO['Gene.refGene']

def write_to_csv(data_matrix, file_name):
    with open(file_name, 'wb') as csvfile:
        file_writer = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for line in data_matrix :
            file_writer.writerow(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files")
    paths = parser.parse_args()
    logging.info('Starting Prediction Pipeline')
    execute_main(paths)
