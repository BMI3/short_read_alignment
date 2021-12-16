# introduction

With the fast developing DNA sequencing technologies, the amount DNA sequencing reads generated has been increasing and price of sequencing has been affordable for labs around the around world. It calls for a fast and also memory-efficient sequence alignment tool to handle with bursting number of sequencing data. To perform quick query, the most commonly used method is to build an index for the reference first. However, unlike virus or other simple-structure organisms, the reference genome of human is quite huge (over 3 billion base pairs, file size around 2.3 GB). If we want to build index for it with traditional methods like suffix tree or suffix array, it will take over 45 GB memory for suffix tree and more than 12 GB memory for suffix array. To solve this problem, many alignment tools developed  their tools using Burrows-Wheeler Transform (BWT) (Burrows  and  Wheeler,  1994), a famous compression algorithm used in string matching. For example, a popular sequencing aligning tools, Bowtie 2, realizes ultrafast and memory-efficient (around 1.5 GB required memory) aligning using BWT algorithm. Inspired by it, we also tried to utilize BWT algorithm and developed a mini-software that can map DNA sequences against a reference genome, such as human chromosome 1 genome with high speed and low memory occupied. To test the performance of our alignment tool, we use a biological sequence simulation tool called Mason to generate random genomic sequences and NGS reads from random generated genomic sequences as input to our algorithm and other naive approaches. The results showed that our alignment tool outperformed naive approaches with a faster speed and lower memory occupation. 

# alogrithm pseudocode & implemention

1. generate bwt string
2.  create FM index
3.  define the function of query

- 思路：先介绍LF Mapping，（same occurence）好处-略过其他字母， 如何实现？（ranks，tots）--- 然而这样

After constructing the BWT, we can build a index called FM index for querying. There are 4 main component of the FM index: First column (F), Last column (L), the position, and the rank of characters.  The first three component has been generated in the first step mentioned above. The rank of characters can be calculated by counting through the last column and first column. 

After building the FM index, we can match reads to the reference genome. Fig X gives an example of match a pattern (P) 'aba' to the reference string 'abaaba' (T). Starting from the first character 'a', we can locate rows in the F having value 'a', according to how we perform BWT, we know the characters before those 'a's are stored in the L. Comparing with the second last character in the P, we can further narrow the scope of the possible position of P in the T (Figxa). A method called LF mapping, which is based the property that the ith occurence of a character c in L will be the same one of the ith occurence of c in the F (i.e. they are from the same position in the original reference string),  will be used to find the corresponding row in the F from L (Fig xb). The remaining character in p can be matched the same way to get finally match the whole pattern in the T(Fig xc), and whenever a character cannot be matched, the querying will quit and giving a not aligned return (-1 in our software).

However, using the naive FM index mentioned above will scan through the matching row in the L to find the preceding character, whose complexity is O(m) (figure x). To increase the speed of query, we improve the way we store the ranks. We generate a dictionary called tally which store the number of all characters has been occurred up to every row in L. Fig x gives a example of tally for L 'abba$aa'. With tally, we can find the corresponding row for the preceding character with only 2 values in the tally (fig x). 

- seeding and extension

A problem for query mentioned above is that it can only do the perfect matching (i.e. no mismatches allowed), which is impossible for real world short reads alignment. Therefore, seeding has been on the reads. According to the pigeonhole rule, if there are k mismatches in the reads, there must be [n/k+1] consecutive characters that can be perfectly matched in that reads. After finding the perfect match, matched patterned were further extended until the length is equal to the original read length.

# results

We simulated the genome and its corresponding short read sequencing data using the Mason program to test our software. Since its simulated data, we set the overall quality of reads high (all bases above 38). The benefit of using simulated data is that we know the exact coordinate of every read generated, so we can calculate the accuracy of the alignment.

- 需要一张图，反映accuracy（？）

To show the optimization of memory and speed, we compare the space and time usage during performing BWT of space saving way used in our software with the traditional way (by building suffix array). Fig x shows that the performance of using "arg_sort" way we introduced is not better than the traditional way when the length of sequencing is lower than ~1500000 bps, but it greatly outperforms when the length increases greater. We also test the speed of query either with or without tally and Fig x shows that the speed deed increases and the saving-time increase while increasing the reference genome length.

- 需要两张图，一张比较构建bwt的string的速率图，一张比较query的速率图

# discussion

Given the results above, we can conclude that our software is able to do a job when aligning short reads sequences to the long reference genome. One of the advantage of our software is saving space by its implementation of BWT. Compared to suffix trees, BWT only store the last column of a sorted matrix of all possible rotations of an original string and its corresponding index **(both up to the length of the original string)**. Moreover, we also optimize the way to construct BWT.  In the traditional, BWT was constructed by rotation and sorting, and this process will generate a matrix with O(n2) space. Our novel algorithm avoid this problem by using only necessary reads required to sort instead of generating a large matrix (details described in the algorithm part).  However, though we improve the way to perform BWT, the length of reference genome is still limited. This is understandable as the length of genome increases, the working memory our innovative method needs to generate is also increases, though using much lower memory than that traditional way uses. A potential idea to solve lies in reducing working memory by an algorithm proposed by Hon et al. (2007). It was a space and time efficient algorithm used to generate compressed suffix array, which can be then used to generate FM index. It has been reported that less than 1 GB memory at most was used to construct BWT of human genome (Hon et al., 2007). 



  Illumina Options:
    --illumina-read-length INTEGER
          Read length for Illumina simulation. In range [1..inf]. Default: 100.
    --illumina-error-profile-file INPUT_FILE
          Path to file with Illumina error profile. The file must be a text file with floating point numbers separated
          by space, each giving a positional error rate. Valid filetype is: .txt.
    --illumina-prob-insert DOUBLE
          Insert per-base probability for insertion in Illumina sequencing. In range [0..1]. Default: 0.00005.
    --illumina-prob-deletion DOUBLE
          Insert per-base probability for deletion in Illumina sequencing. In range [0..1]. Default: 0.00005.
    --illumina-prob-mismatch-scale DOUBLE
          Scaling factor for Illumina mismatch probability. In range [0..inf]. Default: 1.0.
    --illumina-prob-mismatch DOUBLE
          Average per-base mismatch probability in Illumina sequencing. In range [0.0..1.0]. Default: 0.004.
    --illumina-prob-mismatch-begin DOUBLE
          Per-base mismatch probability of first base in Illumina sequencing. In range [0.0..1.0]. Default: 0.002.
    --illumina-prob-mismatch-end DOUBLE
          Per-base mismatch probability of last base in Illumina sequencing. In range [0.0..1.0]. Default: 0.012.
    --illumina-position-raise DOUBLE
          Point where the error curve raises in relation to read length. In range [0.0..1.0]. Default: 0.66.
    --illumina-quality-mean-begin DOUBLE
          Mean PHRED quality for non-mismatch bases of first base in Illumina sequencing. Default: 40.0.
    --illumina-quality-mean-end DOUBLE
          Mean PHRED quality for non-mismatch bases of last base in Illumina sequencing. Default: 39.5.
    --illumina-quality-stddev-begin DOUBLE
          Standard deviation of PHRED quality for non-mismatch bases of first base in Illumina sequencing. Default:
          0.05.
    --illumina-quality-stddev-end DOUBLE
          Standard deviation of PHRED quality for non-mismatch bases of last base in Illumina sequencing. Default:
          10.0.
    --illumina-mismatch-quality-mean-begin DOUBLE
          Mean PHRED quality for mismatch bases of first base in Illumina sequencing. Default: 40.0.
    --illumina-mismatch-quality-mean-end DOUBLE
          Mean PHRED quality for mismatch bases of last base in Illumina sequencing. Default: 30.0.
    --illumina-mismatch-quality-stddev-begin DOUBLE
          Standard deviation of PHRED quality for mismatch bases of first base in Illumina sequencing. Default: 3.0.
    --illumina-mismatch-quality-stddev-end DOUBLE
          Standard deviation of PHRED quality for mismatch bases of last base in Illumina sequencing. Default: 15.0.
    --illumina-left-template-fastq INPUT_FILE
          FASTQ file to use for a template for left-end reads. Valid filetypes are: .sam[.*], .raw[.*], .gbk[.*],
          .frn[.*], .fq[.*], .fna[.*], .ffn[.*], .fastq[.*], .fasta[.*], .faa[.*], .fa[.*], .embl[.*], and .bam, where
          * is any of the following extensions: gz and bgzf for transparent (de)compression.
    --illumina-right-template-fastq INPUT_FILE
          FASTQ file to use for a template for right-end reads. Valid filetypes are: .sam[.*], .raw[.*], .gbk[.*],
          .frn[.*], .fq[.*], .fna[.*], .ffn[.*], .fastq[.*], .fasta[.*], .faa[.*], .fa[.*], .embl[.*], and .bam, where
          * is any of the following extensions: gz and bgzf for transparent (de)compression.