# short_read_alignment



## code structure

```shell
short read alignment

-- data

-- MyBWT.py
	#define the MyBWT class

```

* interface.py
  	1. build genome index

```python
	python interface.py -b -r refernce.fasta
		#example
		python interface.py -b -r .\data\simulated_1m.fa
```

​			2. alignment

```python
	python interface.py -q -r refernce.fasta -sr shortRead.fastq
	python interface.py -q -r refernce.fasta.bwt -sr shortRead.fastq	#recommend
		#use the pre-generated MyBWT object for alignment

        python interface.py -q -r .\data\simulated_1m.fa.bwt -s .\data\pairend_read1_with_virants.fq
```