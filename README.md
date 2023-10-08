# **Sequence Processing**

The module **"Sequence Processing"** is designed to process various biological sequences
It contains:
* **run_aminoacid_seq**
         |this function obtains various characteristics of your amino acid sequence
* **run_dna_rna_tools**
         |this function performs some procedures with nucleotide sequences:
* **run_fastsq**
         |this function filters reads based on the length of their nucleotide sequence, GC-content and phred33-score

To launch the module you need to download folder in your working directory and perform command ``` import sequence_processing.py ```

## **Next functions will be available to you**:
## run_aminoacid_seq

### Overview:
The idea is to make life easier for experimenters and bioinformaticians working with amino acid sequences (both long and short). **run_aminoacid_seq** will allow you to quickly obtain various characteristics of your amino acid sequence. 

### Instruction: 
The function ‘run_aminoacid_seq’ takes 1 amino acid sequence as input (both in three-letter and one-letter form, str), next you need to specify named arguments: function, record_type.
- The argument **sequence** entirely must be written in capital or lowercase letters without spaces or commas. Output saves register
- The argument **“function “** must be passed a string with the name of the action (see possible actions below) that needs to be performed. 
- The argument **“record_type”** indicates the form in which you present your sequence. If you are using a three-letter sequence, specify **“record_type= 3”**. If you are using a one -letter sequence, specify **“record_type= 1”** (by default). 

#### Example:
```run_aminoacid_seq (‘ALAGLNGLU’, function = ‘count_protein_length’, record_type = 1)```        
```run_aminoacid_seq (‘glyvalala’, function = ‘count_protein_length’, record_type = 3)```

Also, if necessary specify a named argument **percent=True** (default False) for actions: determine_charge, determine_polarity (Look in the description of actions).

#### Example:
```run_aminoacid_seq (‘LLYdD’, function = ‘determine_charge’, record_type = 1, percent=True)```

### Troubleshooting:
The program works with 20 proteinogenic amino acids {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W', 'S', 'T', 'N', 'Q', 'Y', 'C', 'K', 'R', 'H', 'D', 'E'}.
If you use a symbol that does not represent an amino acid, or wrong record type, the program will generate an **error** where you can see the first wrong symbol. 
- correct function launch:
```run_aminoacid_seq('ALALEUILE', function = 'count', record_type = 3)```
- incorrect function launch:
```run_aminoacid_seq('ALAmIle', function = 'count', record_type = 3)```

### Possible actions:
1. **translate** - Translation of a one-letter amino acid sequence into a three-letter one (for better visual perception), and the reverse operation. Output: str
2. **count_protein_length** - obtaining the length of the amino acid sequence. Output: int
3. **count_possible_number_of_disulfide_bonds** - counting the number of possible combinations of two different cysteines to form a disulfide bond. Output: int
4. **count_molecular_weight** - calculating the molecular weight of a protein. Output: int
5. **determine_charge** - counting the number of positive, negative and neutral amino acids in a protein. To get the output in percent, specify percent=True. Output: dict
6. **determine_polarity** - counting hybrophobic and hydrophilic amino acids in a protein. To get the output in percent, specify percent=True. Output: dict
7. **convert_amino_acid_seq_to_dna** - convert an amino acid sequence to the most likely DNA sequence. Output: str
8. **summary** - a summary of all information about the sequence (the result of executing all functions). Output: dict

#### Example:

```run_aminoacid_seq (‘LLYdD’, function = ‘translate’, record_type = 1)```         
```run_aminoacid_seq (‘ALAGLYALA’, function = ‘translate’, record_type = 3)```         
```run_aminoacid_seq (‘LLYdD’, function = ‘count_protein_length’, record_type = 1)```          
```run_aminoacid_seq (‘ALAGLYALA’, function = ‘count_protein_length’, record_type = 3)```         
```run_aminoacid_seq (‘LLYdD’, function = ‘count_molecular_weight’, record_type = 1)```          

```run_aminoacid_seq (‘alaglyala’, function = ‘determine_charge', record_type = 3, percent=True)```         
```run_aminoacid_seq (‘LLYdD’, function = ‘determine_charge', record_type = 1, percent=True)```         

```run_aminoacid_seq (‘alaglyala’, function = ‘determine_charge’, record_type = 3)```          
```run_aminoacid_seq (‘LLYdD’, function = ‘determine_charge', record_type = 1)```          

## Development team:

**Iurii Slepov** - team leader, author of main, count, translate and summary functions             
**Veronika Vadekhina** - author of count_possible_number_of_disulfide_bonds, count_molecular_weight and convert_amino_acid_seq_to_dna functions           
**Yulia Nechaeva** - author of determine_charge and determine_polarity functions               

## run_dna_rna_tools

**run_dna_rna_tools.py** is a module that implements the following procedures with nucleotide sequences:
- *transcribe* — print the transcribed sequence
- *reverse* — print a reversed sequence
- *complement* — print the complementary sequence
- *reverse_complement* — print the reverse complementary sequence
- *complement_transcribe* - print the complementary transcribed sequence

### Function launch
To launch a function, you need to pass any number of arguments with DNA or RNA sequences containing 3 nuclуotides (str) in any register, and the name of the procedure to be executed (this is always the last argument, str). All specified elements must be separated by commas and spaces:

``` run_dna_rna_tools('ATG', 'transcribe') # correct function launch ```              
``` run_dna_rna_tools('ATG', 'CTA', 'complement') # correct function launch ```            
``` run_dna_rna_tools('ATgcTA', 'complement') # uncorrect function launch ```                
``` run_dna_rna_tools('ATGCTA', complement) # uncorrect function launch ```              

If your sequence contains the nucleotides "U" and "T", the module cannot distinguish whether it is RNA or DNA. It will result to Error
``` run_dna_rna_tools('ATU', 'transcribe') # uncorrect function launch ```

### Developer:
- Iurii Slepov

## run_fastsq:
Function filters reads based on the length of their nucleotide sequence, GC-content and phred33-score

#### Arguments
* *seqs*: dictionary consisting of fastq sequences
* *gc_bounds*: GC-content interval (percentage) for filtering. Default = 0, 100. If you pass one number as an argument, then it is considered that this is the upper bound.
* *length_bounds*: length interval for filtering. Default = 0, 2**32. If you pass one number as an argument, then it is considered that this is the upper bound.
* *quality_threshold*: threshold value of average read quality for filtering. Default = 0

#### Return
Dictionary consisting of filtered fastq reads matching all filters

#### Exmple
``` run_fastq(dict, gc_bounds = 44.4) # correct launch```             
``` run_fastq(str, gc_bounds = 44.4) # uncorrect launch```            
``` run_fastq(dict, gc_bounds = 44.4, length_bounds = (0, 132)) # correct launch```               
``` run_fastq(dict, gc_bounds = 44.4, length_bounds = 38.53) # uncorrect launch```                 
``` run_fastq(dict, gc_bounds = 44.4, length_bounds = (0, 1023), quality_threshold = (2,20)) # uncorrect launch```                 

#### Troubleshooting:
The function works with **5 nucleotide bases** in any register. If sequence in your read contains any other symbol function will stopped and show you **invalid character**. Please, fixed your sequence and try again.

#### Developer:
- Iurii Slepov








