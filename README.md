# Short project : Embeddings alignment

The aim of this project is to create an embedding alignment program by dynamic programming. The scalar products of each vector between embeddings are calculated and used as a score matrix. Then, the transformed matrix is filled according to the chosen alignment and gap penalties. The alignments are generated as output and are saved in a .txt file.  

Three algorithms are available :
* global (Needleman and Wunsch)
* local (Smith and Waterman) 
* semi-global

By default gap penalties are set to zero but you can choose an affine gap penalty (-1 for a gap opening and 0 for a gap extension)

## 	:zero: Prerequisites

To use the program you must have python. 
To download python: https://www.python.org/downloads/. The version used for this project is 3.9.12.

Clone the repository:

```SHELL
git clone git@github.com:AdrianaLecourieux/embeddings_alignment.git
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

Once all prerequisites have been installed, there are a few files that are necessary before starting. At a bare minimum, you need .t5emb files for your embeddings of interest as well as fasta files.
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

If you need help about inputs, you can use the --help command:

```SHELL
cd src/
python main.py --help
```
```
   -h, --help            show this help message and exit
  -emb1 EMBEDDING1, --embedding1 EMBEDDING1
                        Enter Embedding 1 in .t5emb extension
  -emb2 EMBEDDING2, --embedding2 EMBEDDING2
                        Enter Embedding 2 in .t5emb extension
  -f1 FASTA1, --fasta1 FASTA1
                        Enter Fasta 1 in .FASTA extension
  -f2 FASTA2, --fasta2 FASTA2
                        Enter Fasta 2 in .FASTA extension
  -m METHOD, --method METHOD
                        Choose a "global" (Needleman and Wunsch), "local" (Smith and Waterman) or "semi_global" alignment algorithm. -m global default"
  -g GAP_PENALTY, --gap_penalty GAP_PENALTY
                        Use this option to add affine gap penalty (Enter "affine" to used -1 for gap opening and 0 for gap extension). Else, gap penalty is fixed to 0
```

### :point_right: Global Alignment (Needleman and Wunsch)
* If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m global
```
* If you want to run a global alignment with an affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):

```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m global -g affine
```

### :point_right: Local Alignment (Smith and Waterman)
* If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m local
```
* If you want to run a global alignment with an affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):
```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m local -g affine
```


### :point_right: Semi-global Alignment
* If you want to run a global alignment with a gap penalty fixed to 0:

```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m semi_global
```
* If you want to run a global alignment with an affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):
```SHELL
cd src/
python main.py -emb1 embedding1.t5emb -emb2 embedding2.t5emb -f1 fasta1.fasta -f2 fasta2.fasta -m semi_global -g affine
```

## :four: Collect output

When the alignment is complete, you will see the message :

```SHELL
Alignment completed successfully !
```  

The outputs are in the following path in .txt format :

```SHELL
cd ../results
```

## :five: Example
### :point_right: Global Alignment (Needleman and Wunsch)

* gap penalty fixed to 0:
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5em -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta -m global
```
* affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5emb -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta -m global -g affine
```

### :point_right: Local Alignment (Smith and Waterman)
* gap penalty fixed to 0:
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5emb -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/6PF2K_1BIF.fasta -m local
```
* affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5emb -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta -m local -g affine
```

### :point_right: Semi-global alignment

* gap penalty fixed to 0:
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5emb -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta -m semi_global
```

* affine gap penalty (with the penalties: -1 for a gap opening and 0 for a gap extension):
```SHELL
cd src/
python main.py -emb1 ../data/embeddings/6PF2K_1bif.t5emb -emb2 ../data/embeddings/5_3_exonuclease_1bgxt.t5emb -f1 ../data/fasta_sequences/6PF2K_1BIF.fasta -f2 ../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta -m semi_global -g affine
```
***
