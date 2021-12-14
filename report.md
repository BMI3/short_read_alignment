# introduction

With the fast developing DNA sequencing technologies, the amount DNA sequencing reads generated has been increasing and price of sequencing has been affordable for labs around the around world. It calls for a fast and also memory-efficient sequence alignment tool to handle with bursting number of sequencing data. To perform quick query, the most commonly used method is to build an index for the reference first. However, unlike virus or other simple-structure organisms, the reference genome of human is quite huge (over 3 billion base pairs, file size around 2.3 GB). If we want to build index for it with traditional methods like suffix tree or suffix array, it will take over 45 GB memory for suffix tree and more than 12 GB memory for suffix array. To solve this problem, many alignment tools developed  their tools using Burrows-Wheeler Transform (BWT) (Burrows  and  Wheeler,  1994), a famous compression algorithm used in string matching. For example, a popular sequencing aligning tools, Bowtie 2, realizes ultrafast and memory-efficient (around 1.5 GB required memory) aligning using BWT algorithm. Inspired by it, we also tried to utilize BWT algorithm and developed a mini-software that can map DNA sequences against a reference genome, such as human chromosome 1 genome with high speed and low memory occupied. To test the performance of our alignment tool, we use a biological sequence simulation tool called Mason to generate random genomic sequences and NGS reads from random generated genomic sequences as input to our algorithm and other naive approaches. The results showed that our alignment tool outperformed naive approaches with a faster speed and lower memory occupation. 

# alogrithm pseudocode & implemention

1. generate bwt string
2.  create FM index
3.  define the function of query

# results

We simulated the genome and its corresponding short read sequencing data using the Mason program to test our software. The overall quality of reads is xxx. The benefit of using simulated data is that we know the exact coordinate of every read generated, so we can calculate the accuracy of the alignment.

- 需要一张图，反应accuracy

To show the optimization of memory and speed, we compare the space and time usage during performing BWT of space saving way used in our software with the traditional way (by building suffix array). Fig x shows that the performance of our software is not better than the traditional way when the length of sequencing is lower than xxx, but it greatly outperforms when the length increases greater. 