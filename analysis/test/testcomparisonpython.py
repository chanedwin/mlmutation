import unittest
import machinelearning.helpermethods as helper

class TestFeatureMethods(unittest.TestCase):


    def test_is_vcf_in_dictionary(self):
        test_dictionary_case_one = {(1,1,"A"):("A","C"),(1,1,"C"):("A","G")}
        test_dictionary_case_two = {(3,1,"A"):("A","C"),(2,1,"A"):("A","C")}
        tuple_record_case_one = (1,1,"A",("A","C"))
        tuple_record_case_two = (2,1,"A",("A","C"))
        tuple_record_case_three = (3,1,"A",("A","C"))
        tuple_record_case_four = (1,1,"A",("A",))
        tuple_record_case_five = (1,1,"A",("G",))
        tuple_record_case_six = (3,1,"A",("G","C"))
        tuple_record_case_seven = (3,2,"A",("G","C"))
        tuple_record_case_eight = (3,1,"C",("G","C"))
        self.assertTrue(helper.is_vcf_record_in_dictionary(tuple_record_case_one,test_dictionary_case_one))
        self.assertTrue(helper.is_vcf_record_in_dictionary(tuple_record_case_four,test_dictionary_case_one))
        self.assertTrue(helper.is_vcf_record_in_dictionary(tuple_record_case_six,test_dictionary_case_two))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_one,test_dictionary_case_two))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_two,test_dictionary_case_one))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_three,test_dictionary_case_one))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_five,test_dictionary_case_one))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_five,test_dictionary_case_two))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_four,test_dictionary_case_two))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_seven,test_dictionary_case_two))
        self.assertFalse(helper.is_vcf_record_in_dictionary(tuple_record_case_eight,test_dictionary_case_two))


if __name__ == '__main__':
    unittest.main()