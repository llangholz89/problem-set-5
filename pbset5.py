#! usr/env/bin python

# Problem Set 5

#-------------------
#Problem 1:

#files:
lamina='/Users/loganlangholz/Documents/class_files/Graduate/Spring_2016/Genome_analysis/Projects/data-sets/bed/lamina.bed'

#----------
#Problem 1
#On what chr is region with largest start position in 'lamina.bed' file? 

maxstart = 0
maxend = 0
maxchrom = ''

for line in open(lamina):
    if line.startswith('#'): continue

    fields = line.strip().split('\t')

    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])

    if start > maxstart:
        maxstart = start
        maxchrom = chrom
        maxend = end

print 'answer1- ' + maxchrom
#-----------
#What is the region with the largest end position on chrY in 'lamina.bed'?
#Report as: chrom start end value region length

maxchrom = ''
maxend = 0
maxstart = 0

for line in open(lamina):
    if line.startswith('#'): continue

    fields = line.strip().split('\t')

    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])

    if chrom == 'chrY' : 
        if end > maxend:
            maxend = end
            maxchrom = chrom
            maxstart = start

print 'answer2- ' + 'Chrom: %s\tStart: %d\tEnd: %d\tLenth: %d' %(maxchrom,
maxstart, maxend, maxend - maxstart)

#----------------------------------
#Problem 2:
#Which of the first 10 records has largest number of 'C' residues?

from collections import Counter
filename ='/Users/loganlangholz/Documents/class_files/Graduate/Spring_2016/Genome_analysis/Projects/data-sets/fastq/SP1.fq'


line_num = 0
num_records = 0
seq_num = 0
max_C = 0
max_qual = 0

C_counts = []
records = []
quals = [] #list of quality scores
seq_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} #dict for rev complement
rev_comps = []

for line in open(filename):

    line_type = line_num % 4

    if line_type == 0:
        if seq_num > 9: continue
        name = line.strip()


    elif line_type == 1:
        seq = line.strip()
        if seq_num < 10:
            counts = Counter(seq)
            num_C = counts['C']
            records.append(num_C) #puts value first in tuple so can sort
            records.append(name) #name is second part of tuple

            #------
            #Generate Reverse Complement:
            rev_seq = ''
            index = len(seq) - 1

            for i in range(0, len(seq)):
                rev_seq += seq[index]
                index  -= 1

            rev_comp = ''.join([seq_dict[base] for base in rev_seq])
            rev_comps.append(rev_comp)
            #------

        C_counts.append(records) #append the record tuple to C_counts list
        records = [] #emptys the record list
        seq_num += 1

    elif line_type == 3:
        qual = line.strip()
        qual_sum = 0

        for char in qual:
            qual_sum += ord(char)

        if qual_sum >= max_qual:
            max_qual = qual_sum

    line_num += 1

#Sort the counts, name tuple:
C_counts.sort()
C_counts.reverse()
Max_record = C_counts[0] #count and name of record with most 'C' counts

print 'answer3- ' + Max_record[1]
print 'answer4- ' + str(max_qual)
print 'answer5- ' 
for i in range(1, len(rev_comps) + 1):
    print 'reverse_complement_%s:' %i, '\t', rev_comps[i-1] 
#-------------

#Use "python 'filename' > answers.yml" to print answers to file

