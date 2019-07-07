# Important-one-liners

**Perl one liner to extract sequences by their ID from a FASTA file**
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(id1 id2)}print if $c' fasta.file
```
**SAM to sorted BAM**

samtools view -u N22.sam | samtools sort -o N22.bam

**CpG location on the genome**

grep -v "^>" ../Aiptasia_Genome/New_genome/aiptasia_genome_final.fa  | tr -d "\n" | tr -c "[GCgc]" "\n" | grep -v '^$' | grep -o -i 'CG' | wc -l


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
