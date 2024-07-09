# TCRMatch
#### Create a conda env
```
# if you are running in chinese cluster, you can use 

source ~/anaconda3/etc/profile.d/conda.sh
conda activate /share/home/ravin/anaconda3/envs/tcrmatch

### Using yml file
conda env create -f tcrmatch.yml


```
### clone
```
mkdir test_tcrmatch
cd test_tcrmatch
git clone https://gitlab.xbiome.com/RavinP/tcrmatch.git .

```
### Call in tcrmatch tool

```
python tcrmatch.py -h

```

```
usage: tcrmatch.py [-h] [-tcrbetaseq_file TCRBETASEQ_FILE] [-input_genome INPUT_GENOME] [-mismatch MISMATCH]

Run tcrmatch with Docker

options:
  -h, --help            show this help message and exit
  -tcrbetaseq_file TCRBETASEQ_FILE
                        Path to the tcrbeta sequence file
  -input_genome INPUT_GENOME
                        Input genome file
  -mismatch MISMATCH    Number of mismatch allow in kmer match
```


### Example run

```
python tcrmatch.py -tcrbetaseq_file inputdata/GroupA1.txt \
  -input_genome inputdata/GCA_000006785.2_ASM678v2.fa \
  -mismatch 2
  
```
### Output 
  - **GCA_000006785.2_ASM678v2_cross_match_only_crags.txt**: 
    - This has the list of epitopes that has match between TCR database and bacteria. These epitopes represent potential cross-reactive antigens (CRAgs).
  - **output_GroupA1.txt.tsv** 
  This is the result from TCRmatch tool, which has more information such as 
    - score: TCRmatch score -- how good the input TCR sequences match with TCR sequence in the database.
    - receptor id -- TCR repector id. This can be used to search and find more information in IEDB. 
    - epitopes
    - antigen
    - organism
