import machinelearning.extractfeaturesfromvcf as mlextract
import unittest


class RecordInjection:

    def __init__(self, record, alt):
        assert (type(record) is str)
        assert (type(alt) is list or type(alt) is tuple)
        self.REF = record
        self.ALT = alt

class TestFeatureMethods(unittest.TestCase):


    def test_get_entropy_prob(self):
        self.assertEqual(mlextract.get_entropy_prob("AAA"),(1, 0, 0, 0))
        self.assertEqual(mlextract.get_entropy_prob("ACTG"),(0.25, 0.25, 0.25, 0.25))
        self.assertEqual(mlextract.get_entropy_prob("G"),(0, 0, 0, 1))
        self.assertEqual(mlextract.get_entropy_prob("c"),(0, 1, 0, 0))
        self.assertEqual(mlextract.get_entropy_prob("t"),(0, 0, 1, 0))
        self.assertEqual(mlextract.get_entropy_prob("actgATCG"),(0.25, 0.25, 0.25, 0.25))

    def test_get_gc_content(self):
        self.assertEqual(mlextract.get_gc_content("AAA"),0.0)
        self.assertEqual(mlextract.get_gc_content("ACTG"),0.5)
        self.assertEqual(mlextract.get_gc_content("G"),1.0)
        self.assertEqual(mlextract.get_gc_content("c"),1.0)
        self.assertEqual(mlextract.get_gc_content("t"),0.0)
        self.assertEqual(mlextract.get_gc_content("actgATCG"),0.5)

    def test_get_homo_run(self):
        self.assertEqual(mlextract.get_homo_run("AAA"),3)
        self.assertEqual(mlextract.get_homo_run("ACTG"),1)
        self.assertEqual(mlextract.get_homo_run("G"),1)
        self.assertEqual(mlextract.get_homo_run("C"),1)
        with self.assertRaises(AssertionError):
            mlextract.get_homo_run("t")
        self.assertEqual(mlextract.get_homo_run("ACTGATCG"),1)
        self.assertEqual(mlextract.get_homo_run("ACTCGGGGGGGC"),7)
        self.assertEqual(mlextract.get_homo_run("ACTCGGGCGGGGC"),4)

    def test_has_snp(self):
        test_case_snp_one = RecordInjection("A",["C","G"])
        test_case_snp_two = RecordInjection("T",["G"])
        test_case_snp_three = RecordInjection("G",["GC","G"])
        test_case_not_snp_one = RecordInjection("AA",["C","G"])
        test_case_not_snp_two = RecordInjection("A",["CG","GC"])
        test_case_not_snp_three = RecordInjection("A",["",])
        self.assertEqual(mlextract.has_snp(test_case_snp_one),1)
        self.assertEqual(mlextract.has_snp(test_case_snp_two),1)
        self.assertEqual(mlextract.has_snp(test_case_snp_three),1)
        self.assertEqual(mlextract.has_snp(test_case_not_snp_one),0)
        self.assertEqual(mlextract.has_snp(test_case_not_snp_two),0)
        self.assertEqual(mlextract.has_snp(test_case_not_snp_three),0)

    def test_has_indel(self):
        test_case_not_indel_one = RecordInjection("A",["C","G"])
        test_case_not_indel_two = RecordInjection("T",["G"])
        test_case_indel_one = RecordInjection("AA",["C","G"])
        test_case_indel_two = RecordInjection("A",["CG","GC"])
        test_case_indel_three = RecordInjection("T",["GC","C"])
        test_case_indel_four = RecordInjection("",["C"])
        test_case_indel_five = RecordInjection("G",[""])
        self.assertEqual(mlextract.has_indel(test_case_not_indel_one),0)
        self.assertEqual(mlextract.has_indel(test_case_not_indel_one),0)
        self.assertEqual(mlextract.has_indel(test_case_indel_one),1)
        self.assertEqual(mlextract.has_indel(test_case_indel_two),1)
        self.assertEqual(mlextract.has_indel(test_case_indel_three),1)
        self.assertEqual(mlextract.has_indel(test_case_indel_four),1)
        self.assertEqual(mlextract.has_indel(test_case_indel_five),1)



if __name__ == '__main__':
    unittest.main()