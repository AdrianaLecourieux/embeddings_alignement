# Short project : Embeddings alignment

The aim of this project is to create an embedding alignment program by dynamic programming. The scalar products of each vector between embeddings are calculated and used as a score matrix. Then, the transformed matrix is filled according to the chosen alignment and gap penalties. The alignments are generated as output and are saved in a .txt file.

## 	:zero: Prerequisites

To use the program you must have python. 
To download python: https://www.python.org/downloads/. The version used for this project is 3.9.12.

Clone the repository:

```SHELL
git clone XX
```

Move to the new directory:

```SHELL
cd embeddings_alignement/
```

Install Miniconda :  https://docs.conda.io/en/latest/miniconda.html#windows-installers.
Once Miniconda is installed, install Mamba :

```SHELL
conda install mamba -n base -c conda-forge
```

Create the environment and load it :

```SHELL
mamba env create -f embeddings.yml
conda activate embeddings
```
If you want to deactivate the environment, use the command :

```SHELL
conda deactivate
```

-----------------------

## :one: Preparing files

Once all prerequisites have been installed, there are a few files that are necessary before starting. At a bare minimum, you need t5emb files for your embeddings of interest as well as fasta files.
1033 embeddings and fasta files were prepared in advance, you can see if your proteins of interest are in it. Otherwise you can get the embeddings by the T5 ProtTrans method (https://github.com/agemagician/ProtTrans).

Path to access the embeddings files:
```SHELL
cd data/embeddings
```
Path to access the fasta files:
```SHELL
cd data/fasta_sequences
```

-----------------------

## :two: Running embeddings alignment

### :point_right: Global Alignment (Needleman and Wunsch)
If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
```
```PYTHON
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m global
```
If you want to run a global alignment with an affine gap penalty (with the penlaties: -1 for a gap opening and 0 for a gap extension):

```SHELL
cd src/
```
```PYTHON
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m global -g yes
```

### :point_right: Local Alignment (Smith and Waterman)
If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
```
```PYTHON
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m local
```


### :point_right: Semi-global Alignment
If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
```
```PYTHON
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m semi_global
```

## :four: Collect output
The outputs are in the following path :

```SHELL
cd ../results
```

## :five: Example

***
