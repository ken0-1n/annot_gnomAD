#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import annot_gnomad


class TestFilter(unittest.TestCase):

    def setUp(self):
        self.parser = annot_gnomad.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV1.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test1"])
        args.func(args)

        output_file = tmp_dir+"/genomonsv_test1.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test1.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV1.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test2", "--vcf_filter"])
        args.func(args)

        output_file = tmp_dir+"/genomonsv_test2.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test2.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV2.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test3"])
        args.func(args)

        output_file = tmp_dir+"/genomonsv_test3.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test3.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV3.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test4"])
        args.func(args)

        output_file = tmp_dir+"/genomonsv_test4.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test4.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
        
    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV4.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test5", "--margin", "1"])
        args.func(args)

        output_file = tmp_dir+"/genomonsv_test5.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test5.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
        
    def test6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/genomonSV5.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf.gz"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/genomonsv_test6", "--vcf_filter"])
        args.func(args)
        
        output_file = tmp_dir+"/genomonsv_test6.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/genomonsv_test6.genomonSV.result.gnomad.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
        
    def sample_test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_file = cur_dir + "/../data/HG002_NA24385_son.genomonSV.result.realign3.txt"
        gnomad_vcf = cur_dir + "/../data/gnomad_v2_sv.sites.vcf"
        output_file = tmp_dir + "/HG002_NA24385_son.genomonSV.result.gnomad.txt"
        answer_file = cur_dir + "/../data/HG002_NA24385_son.genomonSV.result.gnomad.answer.txt"
        args = self.parser.parse_args(["comp", "--input_sv", input_file, "--gnomad_vcf", gnomad_vcf, "--out_pref", tmp_dir+"/HG002_NA24385_son"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()
    