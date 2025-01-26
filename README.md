# piPFP
`piPFP` is a tool to gain information about the repetitiveness of a sequence (through the repetitiveness measure $\pi$) and about the openness of a species. 

## Table of contents
<!--ts-->
- [piPFP](#pipfp)
  - [Table of contents](#table-of-contents)
  - [Install](#install)
  - [Usage](#usage)
    - [Calculating pi](#calculating-pi)
    - [Pangenome growth](#pangenome-growth)
<!--te-->

## Install
```bash
git clone https://github.com/simolucaa/piPFP.git
cd piPFP
make
```

## Usage
### Calculating pi
`piPFP` accepts a FASTA file, or a multi-FASTA file as input. You can use the following command:

```bash
./piPFP.py -w [WINDOW] -p [MODULUS] -t [THREADS] -i [INPUT] -o [OUTPUT]
```

If `-t` is not specified, the tool will run in single-thread mode. If `-w` and `-p` are not specified, piPFP will default to $10,100$.

The output is a `.tsv` file containing the following information: Filename, Window Size, Modulus, Dictionary Size, Parse Size, $\pi$.

### Pangenome growth
`piPFP openness` accepts a directory containing only FASTA files as input. You can use the following command:  

```bash
./piPFP -m openness -w [WINDOW] -p [MODULUS] -t [THREADS] -i [INPUT] -o [OUTPUT]
```

If `-t` is not specified, the tool will run in single-thread mode. If `-w` and `-p` are not specified, piPFP will default to $4,10$.

This command will output 3 files:
- A `.hist` file containing the histogram;
- A `.growth` file containing the growth;
- A `.tsv` file containing the following information: Filename, Window Size, Modulus, $\alpha$.