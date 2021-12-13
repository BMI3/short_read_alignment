# introduction

Introduction 

Burrows–Wheeler transform (BWT) is a well-developed compression algorithm which gathers the same characters to compress them. Besides, it can be used in the genomics filed to improve sequence alignment. This is because the reference genome of human is quite huge (over 3 billion base pairs, file size around 2.3 GB). If we want to build index for it with traditional methods like suffix tree or suffix array, it will take over 45 GB memory for suffix tree and more than 12 GB memory for suffix array. Moreover, with the help of next generation sequencing technologies, the number of short reads needed to be aligned is increasing quickly. This calls for a high-speed alignment tool. Such big memory and high-speed requirement can be solved with the help of BWT. For example, a popular sequencing aligning tools, Bowtie 2 and BWA, realize ultrafast and memory-efficient (around 1.5 GB required memory) aligning using BWT algorithm. Inspired by it, we also tried to utilize BWT algorithm and developed a mini-software that can map DNA sequences against a reference genome, such as human chromosome 1 genome with high speed and low memory occupied. To test the performance of our alignment tool, we use a biological sequence simulation tool called Mason to generate random genomic sequences and NGS reads from random generated genomic sequences as input to our algorithm and other naive approaches. The results showed that our alignment tool outperformed naive approaches with a faster speed and lower memory occupation. 

Since BWT groups characters together, it facilitates search the beginning of a string. Compared to suffix trees, BWT only store the last column of a sorted matrix of all possible rotations of an original string and an index (both up to the length of the original string). However, generate the matrix cost O(n2) time and space, which we want to avoid, especially reference gene is always over million base pair. The last column consists of the same set of characters as the original string (plus a “$” for convenience) does and the meaning of generating the N*N matrix is to decide the order of characters in the last column. However, not the whole matrix contributes to the sorting process. Therefore, we came up with an idea to sort the last column using only necessary information. 

 

Generally, we have 2 types of input, real sequencing data from exist database, and simple stimulate strings. 

 

The system detects the biological type of input data as well. If 