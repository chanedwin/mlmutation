import vcf
from Bio import SeqIO
import scipy.stats

import warnings
import functools
import sklearn

#File containing methods used to extract features from pyvcf record entries.


ENTROPY_CONSTANT_RANGE = 5


def getallvalues(record, reference_dictionary, base_entropy, file_name):
    is_snp = has_snp(record)
    is_indel = has_indel(record)
    record.ALT = list(filter(lambda x : ">" or "<" not in x, record.ALT)) #remove illegal entries
    temp_list = []
    for item in record.ALT:
        temp_list.append(str(item).upper())
    sample_entropy, kl_entropy, gc_content, homo_run, alt_div = extract_information_from_dna_string(record, reference_dictionary, base_entropy)
    # get indels
    indel_list = []
    if "breakseq" in file_name:
        indel_list = bs2_parse_indel(record)
    elif "pind" in file_name:
        indel_list = pindel_parse_indel(record)
    elif "st" in file_name:
        indel_list = st_parse_indel(record)
    elif "ug" in file_name:
        indel_list = ug_parse_indel(record)
    elif "hc" in file_name:
        indel_list = hc_parse_indel(record)
    elif "fb" in file_name:
        indel_list = fb_parse_indel(record)
    fullinfo = [1, kl_entropy, gc_content, homo_run, alt_div, is_snp, is_indel]
    fullinfo.extend(indel_list)
    fullinfo = [x[0] if isinstance(x, list) else x for x in fullinfo]
    fullinfo = [0 if x is None else x for x in fullinfo]
    return fullinfo


def extract_information_from_dna_string(record, reference_dictionary, base_entropy):
    index = record.POS - 1

    #use the index to get the relevant DNA strings
    base_dna_sequence = get_base_dna_sequence(index, record, reference_dictionary)
    adjusted_dna_sequence = dna_sequence_adjustment(base_dna_sequence, record, reference_dictionary, index)

    #get information about the adjusted_dna_sequence
    gc_content = get_gc_content(adjusted_dna_sequence)
    sample_probability = get_entropy_prob(adjusted_dna_sequence)
    homo_run = get_homo_run(adjusted_dna_sequence)
    sample_entropy = scipy.stats.entropy(sample_probability)
    kl_entropy = scipy.stats.entropy(sample_entropy, qk=base_entropy)

    # alt div measures the informational difference between reference and alternate sequences
    alt_ref_info_divergence = extract_alt_div_variable(base_entropy, kl_entropy, record)

    return sample_entropy, kl_entropy, gc_content, homo_run, alt_ref_info_divergence


def get_base_dna_sequence(index, record, reference_dictionary):
    base_string = reference_dictionary[str(record.CHROM)].seq[
                  index - ENTROPY_CONSTANT_RANGE:index + ENTROPY_CONSTANT_RANGE]
    return base_string


def extract_alt_div_variable(base_entropy, kl_entropy, record):
    for alternative in record.ALT:
        if not str(alternative):
            alt_kl_entropy = 0
        else:
            alt_kl_entropy = scipy.stats.entropy(get_entropy_prob(str(alternative)), qk=base_entropy)
        temp_alt_div = abs(kl_entropy - alt_kl_entropy)
        alt_div = max(alt_div, temp_alt_div)
    return alt_div


def get_entropy_prob(sample_string):
    a_count, c_count, t_count, g_count = 0 , 0 , 0 , 0
    a_count += sample_string.count('A') + sample_string.count('a')
    c_count += sample_string.count('C') + sample_string.count('c')
    t_count += sample_string.count('T') + sample_string.count('t')
    g_count += sample_string.count('G') + sample_string.count('g')

    all__values = a_count + c_count + t_count + g_count
    all__values = float(all__values)
    if not all__values:
        return (0.25, 0.25, 0.25, 0.25)
    a_prob = a_count / all__values
    c_prob = c_count / all__values
    t_prob = t_count / all__values
    g_prob = g_count / all__values

    return (a_prob, c_prob, t_prob, g_prob)


def get_gc_content(sample_string):
    a_count, c_count, t_count, g_count = 0 , 0 , 0 , 0
    a_count += sample_string.count('A') + sample_string.count('a')
    c_count += sample_string.count('C') + sample_string.count('c')
    t_count += sample_string.count('T') + sample_string.count('t')
    g_count += sample_string.count('G') + sample_string.count('g')

    all__values = a_count + c_count + t_count + g_count
    all__values = float(all__values)
    if not all__values:
        return 0
    return float(g_count + c_count) / all__values


def get_homo_run(sample_string):
    if not sample_string:
        return 0
    assert sample_string.isupper()
    curr_letter = sample_string[0]
    max_count = 0
    count = 0
    for letter in sample_string:
        if curr_letter == letter:
            count += 1
        else:
            curr_letter = letter
            count = 1
        max_count = max(max_count, count)
    return max_count


def dna_sequence_adjustment(string, record, reference_dictionary, index):
    maxlength = 0
    for alternate in record.ALT:
        maxlength = max(maxlength, len(alternate))
    if maxlength < len(string):
        return string
    else:
        return reference_dictionary[str(record.CHROM)].seq[index - (maxlength / 2):index + (maxlength / 2)]


def get_ref_entropy(path):
    a_count, c_count, t_count, g_count = 0 , 0 , 0 , 0

    for record in SeqIO.parse(path, "fasta"):
        a_count += sample_string.count('A') + sample_string.count('a')
        c_count += sample_string.count('C') + sample_string.count('c')
        t_count += sample_string.count('T') + sample_string.count('t')
        g_count += sample_string.count('G') + sample_string.count('g')

    all__values = a_count + c_count + t_count + g_count

    all__values = float(all__values)

    a_prob = a_count / all__values
    c_prob = c_count / all__values
    t_prob = t_count / all__values
    g_prob = g_count / all__values

    return (a_prob, c_prob, t_prob, g_prob)

def has_snp(record):
    if len(record.REF) == 1:
        for alternate in record.ALT:
            if "<" in str(alternate) or ">" in str(alternate) :
                return 0
            elif len(alternate) == 1:
                return 1
    return 0


def has_indel(record):
    for alternate in record.ALT:
        if (len(alternate) > len(record.REF)) or (len(alternate) < len(record.REF)):
            return 1
    return 0

def get_chr(record):
    return record.name

#deprecated
def bs2_parse_indel(record):
    if True:
        return record
    # checked
    fullinfo = []
    #if "<DEL>" in record.ALT:
    #    fullinfo.extend([0, 0, 1, 0])
    #elif "INS" in record.ALT:
    #    fullinfo.extend([0, 1, 0, 0])
    #else:
    #    fullinfo.extend([0, 0, 0, 0])
    return fullinfo


def fb_parse_indel(record):
    fullinfo = []
    fullinfo.append(record.INFO['DP'])
    fullinfo.append(record.INFO['DPB'])
    fullinfo.append(record.samples[0].data.DPR[0])
    fullinfo.append(record.samples[0].data.DPR[1])
    fullinfo.append(record.INFO['MQM'])
    if ['MQMR'] in record.INFO.keys():
        fullinfo.append(record.INFO['MQMR'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.QUAL)
    fullinfo.append(record.INFO['QA'])
    fullinfo.append(record.INFO['QR'])
    fullinfo.append(record.INFO['AB'])
    fullinfo.append(record.INFO['ABP'])
    GL_score = 0.5*record.samples[0].data.GL[1] + record.samples[0].data.GL[2]
    fullinfo.append(GL_score)
    fullinfo = [x[0] if (type(x) is list) else x for x in fullinfo]
    return fullinfo


def hc_parse_indel(record):
    fullinfo = []
    fullinfo.append(record.INFO['DP'])
    if hasattr(record.samples[0].data, 'AD'):
        fullinfo.extend(record.samples[0].data.AD)
    else:
        fullinfo.extend([0, 0])
    fullinfo.append(record.INFO['MQ'])
    if 'MQRankSum' in record.INFO.keys():
        fullinfo.append(record.INFO['MQRankSum'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.QUAL)
    fullinfo.append(record.samples[0].data.GQ)
    if 'BaseQRankSum' in record.INFO.keys():
        fullinfo.append(record.INFO['BaseQRankSum'])
    else:
        fullinfo.append(0)
    if 'QD' in record.INFO.keys():
        fullinfo.append(record.INFO['QD'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.INFO['AC'])
    fullinfo.append(record.INFO['AF'])
    PL_score = 0.5*record.samples[0].data.PL[1] + record.samples[0].data.PL[2]
    fullinfo.append(PL_score)
    fullinfo = [x[0] if (type(x) is list) else x for x in fullinfo]
    return fullinfo


def ug_parse_indel(record):
    fullinfo = []
    fullinfo.append(record.INFO['DP'])
    fullinfo.extend(record.samples[0].data.AD)
    fullinfo.append(record.INFO['MQ'])
    fullinfo.append(record.INFO['MQ0'])
    if 'MQRankSum' in record.INFO.keys():
        fullinfo.append(record.INFO['MQRankSum'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.QUAL)
    fullinfo.append(record.samples[0].data.GQ)
    if 'BaseQRankSum' in record.INFO.keys():
        fullinfo.append(record.INFO['BaseQRankSum'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.INFO['QD'])
    fullinfo.append(record.INFO['AC'])
    fullinfo.append(record.INFO['AF'])
    PL_score = 0.5 * record.samples[0].data.PL[1] + record.samples[0].data.PL[2]
    fullinfo.append(PL_score)
    fullinfo = [x[0] if (type(x) is list) else x for x in fullinfo]
    return fullinfo


def pindel_parse_indel(record):
    fullinfo = []
    fullinfo.extend(record.samples[0].data.AD)
    return fullinfo


def st_parse_indel(record):
    fullinfo = []
    fullinfo.append(record.INFO['DP'])
    if ['IDV'] in record.INFO.keys():
        fullinfo.append(record.INFO['IDV'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.INFO['MQ'])
    if ['MQB'] in record.INFO.keys():
        fullinfo.append(record.INFO['MQB'])
    else:
        fullinfo.append(0)
    fullinfo.append(record.QUAL)
    PL_score = 0.5*record.samples[0].data.PL[1] + record.samples[0].data.PL[2]
    fullinfo.append(PL_score)
    DP4 = record.INFO['DP4'][0] + record.INFO['DP4'][1] + record.INFO['DP4'][2] + record.INFO['DP4'][3]
    fullinfo.append(DP4)
    return fullinfo

