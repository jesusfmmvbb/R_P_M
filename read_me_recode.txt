If genotype for your data are numerical values (fragment length, number of repeats, etc), the program "recode" can recode your data to fit with requirements of the programme.

Name of the input file should be ***_raw.txt, where *** is a user defined name the programme ask for. This file should include a first row with the number of the genotyped individuals, the number of loci and the maximum number of alleles you expect in any of the loci (if not sure put a large number below 1,000).
Next rows should contain the genotype for each individual in the typical format of first allele first locus, second allele first locus, first allele second locus, second allele second locus, ... (see the following example). Missing values must be coded as 0 (zero).

402 6 100
122	126	222	274	278	278	128	130	191	196	212	216
148	150	230	235	274	322	134	136	0	0	212	224
126	161	222	250	266	274	0	0	0	0	0	0
122	126	230	274	322	322	110	130	191	0	224	232
126	126	222	248	286	322	130	130	0	0	224	237
126	127	248	252	274	286	130	136	0	0	220	228
148	148	222	252	294	298	110	136	182	182	212	228
0	0	236	244	274	298	0	0	0	0	0	0
126	147	235	243	281	285	0	0	0	0	0	0
0	0	235	243	274	322	110	130	0	0	216	220

Output file will be called ***_recod.txt and will include the number of individuals, the number of loci, the actual number of alleles per locus and the recoded genotypes in the same format as above.

         402
           6
          30          39          52          24          20          22
  3  5  4 27 15 15  8  9 13 15  6  7
 19 20  6  9 13 43 12 14  0  0  6 12
  5 24  4 19  6 13  0  0  0  0  0  0
  3  5  6 27 43 43  1  9 13  0 12 17
  5  5  4 18 22 43  9  9  0  0 12 21
  5  6 18 21 13 22  9 14  0  0 10 15
 19 19  4 21 27 30  1 14  9  9  6 15
  0  0 10 15 13 30  0  0  0  0  0  0
  5 18  9 14 18 21  0  0  0  0  0  0
  0  0  9 14 13 43  1  9  0  0  7 10
