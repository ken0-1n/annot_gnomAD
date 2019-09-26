# annot_gnomad
Annotate structural variants in gnomAD

## Dependency

### Python
Python (>= 3.6), [cyvcf2](https://github.com/brentp/cyvcf2)

### Software
[bedtools](http://bedtools.readthedocs.org/en/latest/])

### Database
[gnomAD gnomad_v2_sv.sites.vcf.gz](https://gnomad.broadinstitute.org/downloads)

## Install

```
python setup.py install
```

## Commands

```
annot_gnomad comp [-h] --input_sv INPUT_SV --gnomad_vcf GNOMAD_VCF
                   --out_pref OUT_PREF [--margin MARGIN] [--vcf_filter]
```

You can check the manual by typing
```
annot_gnomad comp -h
```

## Results

The primary result is ${out_pref}.genomonSV.result.gnomad.txt
