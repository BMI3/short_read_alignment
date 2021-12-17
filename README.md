# short_read_alignment

## How to use interface.py

1. build genome index

```python
python interface.py -b -r refernce.fasta

#example
python interface.py -b -r .\data\simulated_1m.fa
```

​	2. alignment

```python
python interface.py -q -r path\to\reference\file -s path\to\sequencing\file
#or
python interface.py -q -r path\to\reference_index\file -s path\to\sequencing\file
	#recommend,use the pre-generated MyBWT object for alignment

#example
python interface.py -q -r .\data\simulated_1m.fa.bwt -s .\data\singleend_read_with_virants.fq
```

## Structure

```shell
-- data

-- MyBWT.py
	#define the MyBWT class
	
	#attributes
	id				#genome id
	characters		#characters appeared ("$ATCG")
	lc				#the BWT last column
	po				#lc corrsponding position in genome
	tally			#the number of each characters appeared in the last column
	char_range		#the distribution of characters in the first column
	origin			#the sequence of the genome
	
	#methods
	__init__()
	generateLC()
	-- my_sort()
	generateTally()
	charRange()
	query()
	-- findNextWithTally()
	seeding()
	extend()
    
-- interface.py
	#define the actions users can do with the software
	
-- test_cost.py
	evalue()		#decorator
	
	cost_generate_LC()
	graph_generate_LC()
	test_generate_LC()
	
	cost_query()
	graph_query()
	test_query()

-- test_generator.py
	generate_DNA(length)
	extract(seq,lm,n)
```

