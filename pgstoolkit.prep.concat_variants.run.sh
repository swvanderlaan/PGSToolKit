#!/bin/bash
# Array job: list variants
STUDYDIR=$1
# 8 bit
echo "variantid alt_ids rsid_aegs chromosome position number_of_alleles first_allele alt_alleles" > ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.variantlist.txt
for CHR in $(seq 1 22); do 
	cat ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.chr${CHR}.variantlist.txt | \
	awk '{ print $3":"$4, $1, $2, $3, $4, $5, $6, $7}' | tail -n +2 >> ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.variantlist.txt 
done
gzip -v ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.8bit.variantlist.txt

# 16 bit
echo "variantid alt_ids rsid_aegs chromosome position number_of_alleles first_allele alt_alleles" > ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.variantlist.txt
for CHR in $(seq 1 22); do 
	cat ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.chr${CHR}.variantlist.txt | \
	awk '{ print $3":"$4, $1, $2, $3, $4, $5, $6, $7}' | tail -n +2 >> ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.variantlist.txt 
done
gzip -v ${STUDYDIR}/aegscombo.topmed_r3_f10_b38.split_norm_af_filter.16bit.variantlist.txt