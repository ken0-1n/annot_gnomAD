from cyvcf2 import VCF
import sys, os
import subprocess


def get_strand(sv_type):
    
    taple1 = ('+','-')
    taple2 = ('-','+')
    taple3 = ('+','+')
    taple4 = ('-','-')
    
    ret = []
    if sv_type == "BND":
        ret = [taple1,taple2,taple3,taple4]
    elif sv_type == "CPX":
        ret = []
    elif sv_type == "CTX":
        ret = []
    elif sv_type == "DEL":
        ret = [taple1]
    elif sv_type == "DUP":
        ret = [taple2]
    elif sv_type == "INS":
        ret = []
    elif sv_type == "INV":
        ret = [taple3,taple4]
    elif sv_type == "MCNV":
        ret = []
    return ret


def convert_vcf_to_bedpe(input_file, output_file, margin_major, margin_minor, method, is_filter):

    hOUT = open(output_file, 'w')
    vcf = VCF(input_file, "r")
    for variant in vcf:
        if variant.FILTER != None and is_filter: continue # a value of PASS in the VCF will give None 

        chr1, POS = variant.CHROM, variant.POS
        chr2, END  = variant.INFO.get("CHR2"), variant.INFO.get("END")
        svtype, svlen = variant.INFO.get("SVTYPE"), variant.INFO.get("SVLEN")
        ac = variant.INFO.get("AC") # Number of non-reference alleles observed
        af = variant.INFO.get("AF") # Allele frequency
        an = variant.INFO.get("AN") # Total number of alleles genotyped 
        qual = variant.QUAL
        start1, end1, start2, end2 = POS, POS, END, END
        chr1 = chr1.replace("chr", "")
        chr2 = chr2.replace("chr", "")

        for strand_pair in get_strand(svtype):
            dir1 = strand_pair[0]
            dir2 = strand_pair[1]
            if dir1 == '+':
                start1 = str(int(start1) - int(margin_minor))
                end1 = str(int(end1) + int(margin_major))
            else:
                start1 = str(int(start1) - int(margin_major))
                end1 = str(int(end1) + int(margin_minor))
            if dir2 == '+':
                start2 = str(int(start2) - int(margin_minor))
                end2 = str(int(end2) + int(margin_major))
            else:
                start2 = str(int(start2) - int(margin_major))
                end2 = str(int(end2) + int(margin_minor))

            ID = chr1 + ':' + dir1 + start1 + '-' + chr2 + ':' + dir2 + start2
            
            # Chr_1    Start_1    End_1
            # Chr_2    Start_2    End_2
            # Name  Score
            # Strand1   Strand2
            # AN, AC, AF sv-type
            print("\t".join([chr1, str(start1), str(end1), chr2, str(start2), str(end2),
                ID, str("{:.0f}".format(qual)), dir1, dir2, str(an), str(ac), str("{:.5f}".format(af))]), file = hOUT)
            
    vcf.close()
    hOUT.close()


def convert_to_bedpe(input_file, output_file, margin_major, margin_minor, method):

    hOUT = open(output_file, 'w')
    with open(input_file, 'r') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Chr_1" and method == "genomonSV": continue

            chr1, pos1, dir1, chr2, pos2, dir2 = F[0], F[1], F[2], F[3], F[4], F[5]
            start1, end1, start2, end2 = pos1, pos1, pos2, pos2
            ID = chr1 + ':' + dir1 + start1 + '-' + chr2 + ':' + dir2 + start2

            if dir1 == '+':
                start1 = str(int(start1) - int(margin_minor))
                end1 = str(int(end1) + int(margin_major))
            else:
                start1 = str(int(start1) - int(margin_major))
                end1 = str(int(end1) + int(margin_minor))
            if dir2 == '+':
                start2 = str(int(start2) - int(margin_minor))
                end2 = str(int(end2) + int(margin_major))
            else:
                start2 = str(int(start2) - int(margin_major))
                end2 = str(int(end2) + int(margin_minor))

            print("\t".join([chr1, start1, end1, chr2, start2, end2, ID, '.', dir1, dir2]), file = hOUT)

    hOUT.close()


def comp_main(args):
    
    convert_to_bedpe(args.input_sv, args.out_pref + ".sv1.bedpe", args.margin, args.margin, "genomonSV")
    convert_vcf_to_bedpe(args.gnomad_vcf, args.out_pref + ".sv2.bedpe", args.margin, args.margin, "gnomAD", args.vcf_filter)
    
    hOUT = open(args.out_pref + ".sv_comp.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.out_pref + ".sv1.bedpe", "-b", args.out_pref + ".sv2.bedpe"], stdout = hOUT)
    hOUT.close()

    # create dictionary
    sv_comp = {}
    with open(args.out_pref + ".sv_comp.bedpe", 'r') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            sv_comp[F[6]] = "\t".join([F[20],F[21],F[22]])
            # print('\t'.join(F) + "\t"+ "\t".join([F[19],F[20],F[21]]))

    # add SV annotation to sv
    hOUT = open(args.out_pref + ".genomonSV.result.gnomad.txt", 'w')
    with open(args.input_sv, 'r') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            
            if F[0] == "Chr_1":
                SV_info = "gnomAD_AN\tgnomAD_AC\tgnomAD_AF"
            else:
                chr1, pos1, dir1, chr2, pos2, dir2 = F[0], F[1], F[2], F[3], F[4], F[5]
                start1, end1, start2, end2 = pos1, pos1, pos2, pos2
                ID = chr1 + ':' + dir1 + start1 + '-' + chr2 + ':' + dir2 + start2
                SV_info = sv_comp[ID] if ID in sv_comp else "---\t---\t---"
            print('\t'.join(F) + '\t' + SV_info, file=hOUT)
    hOUT.close()

    # remove intermediate files
    os.remove(args.out_pref + ".sv1.bedpe")
    os.remove(args.out_pref + ".sv2.bedpe")
    os.remove(args.out_pref + ".sv_comp.bedpe")
    
