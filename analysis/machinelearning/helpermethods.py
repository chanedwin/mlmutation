def is_vcf_record_in_dictionary(sample_tuple, reference_dictionary):
    assert type(sample_tuple) is tuple
    assert len(sample_tuple) == 4
    temp_tuple = (sample_tuple[0], sample_tuple[1], sample_tuple[2])
    if temp_tuple in reference_dictionary:
        for alternate in sample_tuple[3]:
            if alternate in reference_dictionary[temp_tuple]:
                return True
    return False
