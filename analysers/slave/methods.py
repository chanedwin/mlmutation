import vcf
from Bio import SeqIO
from scipy.stats import entropy


def getallvalues(record, record_list, base_entropy, file_name):
    is_snp = has_snp(record)
    is_indel = has_indel(record)
    for alternate in record.ALT:
       if "<" in str(alternate) or ">" in str(alternate) :
            alternate = ""
    temp_list = []
    for item in record.ALT:
        temp_list.append(str(item).upper())
    sample_entropy, kl_entropy, gc_content, homo_run, alt_div = cal_entropy(record, record_list, base_entropy)
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


def cal_entropy(record, record_list, base_entropy):
    index = record.POS - 1
    string = record_list[str(record.CHROM)].seq[index - 5:index + 5]
    string = string_adjustment(string, record, record_list, index)

    gc_content = get_gc_content(string)

    sample_entropy = get_entropy_prob(string)

    homo_run = get_homo_run(string)
    # alt div measures the informational difference between reference and alternate sequences
    alt_div = 0
    kl_entropy = entropy(sample_entropy, qk=base_entropy)
    for alternative in record.ALT:
        if str(alternative) == "":
            alt_kl_entropy = 0
        else:
            alt_kl_entropy = entropy(get_entropy_prob(str(alternative)), qk=base_entropy)
        temp_alt_div = abs(kl_entropy - alt_kl_entropy)
        alt_div = max(alt_div, temp_alt_div)

    return entropy(sample_entropy), kl_entropy, gc_content, homo_run, alt_div


def get_entropy_prob(string):
    a_count = 0
    c_count = 0
    t_count = 0
    g_count = 0
    a_count += string.count('A')
    a_count += string.count('a')
    c_count += string.count('C')
    c_count += string.count('c')
    t_count += string.count('T')
    t_count += string.count('t')
    g_count += string.count('G')
    g_count += string.count('g')

    all__values = a_count + c_count + t_count + g_count
    all__values = float(all__values)
    if all__values == 0:
        return (0.25, 0.25, 0.25, 0.25)

    a_prob = a_count / all__values
    c_prob = c_count / all__values
    t_prob = t_count / all__values
    g_prob = g_count / all__values

    return (a_prob, c_prob, t_prob, g_prob)


def get_gc_content(string):
    a_count = 0
    c_count = 0
    t_count = 0
    g_count = 0
    a_count += string.count('A')
    a_count += string.count('a')
    c_count += string.count('C')
    c_count += string.count('c')
    t_count += string.count('T')
    t_count += string.count('t')
    g_count += string.count('G')
    g_count += string.count('g')

    all__values = a_count + c_count + t_count + g_count
    all__values = float(all__values)
    if all__values == 0:
        return 0
    return (g_count + c_count) / all__values


def get_homo_run(string):
    if string == "" :
        return 0
    curr_letter = string[0]
    maxcount = 0
    count = 0
    for letter in string:
        if curr_letter == letter:
            count += 1
        else:
            curr_letter = letter
            count = 1
        maxcount = max(maxcount, count)
    return maxcount


def string_adjustment(string, record, record_list, index):
    maxlength = 0
    for alternate in record.ALT:
        maxlength = max(maxlength, len(alternate))
    if maxlength < len(string):
        return string
    else:
        return record_list[str(record.CHROM)].seq[index - (maxlength / 2):index + (maxlength / 2)]


def get_ref_entropy(path):
    a_count = 0
    c_count = 0
    t_count = 0
    g_count = 0

    for record in SeqIO.parse(path, "fasta"):
        a_count += record.seq.count('A')
        c_count += record.seq.count('C')
        t_count += record.seq.count('T')
        g_count += record.seq.count('G')
        a_count += record.seq.count('a')
        c_count += record.seq.count('c')
        t_count += record.seq.count('t')
        g_count += record.seq.count('g')

    all__values = a_count + c_count + t_count + g_count

    all__values = float(all__values)

    a_prob = a_count / all__values
    c_prob = c_count / all__values
    t_prob = t_count / all__values
    g_prob = g_count / all__values

    return (a_prob, c_prob, t_prob, g_prob)


def get_chr(record):
    return record.name


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
    #temp = record.INFO["TYPE"][0]
    #if "np" in temp:
    #    fullinfo.extend([1, 0, 0, 0])
    #elif "ins" in temp:
    #    fullinfo.extend([0, 1, 0, 0])
    #elif "del" in temp:
    #    fullinfo.extend([0, 0, 1, 0])
    #elif "complex" in temp:
    #    fullinfo.extend([0, 0, 0, 1])
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
    #if len(record.REF) == len(record.ALT[0]):
    #    fullinfo.extend([1, 0, 0, 0])
    #elif len(record.REF) < len(record.ALT[0]):
    #    fullinfo.extend([0, 1, 0, 0])
    #elif len(record.REF) > len(record.ALT[0]):
    #    fullinfo.extend([0, 0, 1, 0])
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
    #if len(record.REF) == len(record.ALT[0]):
    #    fullinfo.extend([1, 0, 0, 0])
    #elif len(record.REF) < len(record.ALT[0]):
    #    fullinfo.extend([0, 1, 0, 0])
    #elif len(record.REF) > len(record.ALT[0]):
    #    fullinfo.extend([0, 0, 1, 0])
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
    # checked
    # print record.INFO["SVTYPE"]
    #if "RPL" in record.INFO["SVTYPE"]:
        #   print "RPL was returned"
    #    fullinfo.extend([1, 0, 0, 0])
    #elif "INS" in record.INFO["SVTYPE"]:
    #    fullinfo.extend([0, 1, 0, 0])
        # print "INS was returned"
    #elif "DEL" in record.INFO["SVTYPE"]:
    #    fullinfo.extend([0, 0, 1, 0])
    #elif "INV" in record.INFO["SVTYPE"]:
    #    fullinfo.extend([0, 0, 0, 1])
    #elif "DUP" in record.INFO["SVTYPE"]:
    #    fullinfo.extend([0, 0, 0, 1])
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
    #if len(record.REF) == len(record.ALT[0]):
    #    fullinfo.extend([1, 0, 0, 0])
    #elif len(record.REF) < len(record.ALT[0]):
    #    fullinfo.extend([0, 1, 0, 0])
    #elif len(record.REF) > len(record.ALT[0]):
    #    fullinfo.extend([0, 0, 1, 0])
    #fullinfo.append(record.INFO['SGB'])
    #if ['INDEL'] in record.INFO.keys():
    #    fullinfo.append(record.INFO['INDEL'])
    #else:
    #    fullinfo.append(0)
    PL_score = 0.5*record.samples[0].data.PL[1] + record.samples[0].data.PL[2]
    fullinfo.append(PL_score)
    DP4 = record.INFO['DP4'][0] + record.INFO['DP4'][1] + record.INFO['DP4'][2] + record.INFO['DP4'][3]
    fullinfo.append(DP4)
    return fullinfo


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
