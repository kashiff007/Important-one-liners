# Important-one-liners

**Download all reads from a NCBI project.**
```
module load sratoolkit edirect
esearch -db sra -query PRJNA40075  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files
```
Or use prefetch for individual SRR* file
```
prefetch -v SRR5956435
fastq-dump --split-files SRR5956435
```
NCBI documents for downloading data: https://www.ncbi.nlm.nih.gov/books/NBK179288/

**Perl one liner to extract sequences by their ID from a FASTA file**
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(id1 id2)}print if $c' fasta.file
```
**SAM to sorted BAM**
```
samtools view -@ 40 -u Wheat_E1_t1S_S8.sam | samtools sort -@ 40 -o Wheat_E1_t1S_S8.bam 
```
**CpG location on the genome**
```
grep -v "^>" genome.fa  | tr -d "\n" | tr -c "[GCgc]" "\n" | grep -v '^$' | grep -o -i 'CG' | wc -l
```
**Extract all lines from a file which is present al list in gene.list**
```
grep -Fw -f gene.list Tair10_GFF3_genes.gff
```
**Get Multi Fasta Length**
```
samtools faidx sample.fa

cut -f1-2 sample.fa.fai
```
**Sum/average a column awk**
```
awk -F',' '{sum+=$57;} END{print sum;}' file.txt
awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' file.txt
```
## R dataframe column/row selection 

|    Selections                                  |    Commands  |
|------------------------------------------------|-------------:|
|All Rows and All Columns                        |df[,]         |
|First row and all columns                       |df[1,]        |
|First two rows and all columns                  |df[1:2,]      |
|First and third row and all columns             |df[ c(1,3), ] |
|First Row and 2nd and third column              |df[1, 2:3]    |
|First, Second Row and second and third column   |df[1:2, 2:3]  |
|Just First Column with All rows                 |df[, 1]       |
|First and Third Column with All rows-           |df[,c(1,3)]   |
