cmaple
=======


Command reference
-----------------

All the options available in CMAPLE are shown below:

| Option | Usage and meaning |
|--------------|------------------------------------------------------------------------------|
|`-h` or `-?`| Print help usage. |
|`-aln <ALIGNMENT>`| Specify an input alignment file in PHYLIP, FASTA, or [MAPLE](MAPLE_FORMAT) format. |
| `-m <MODEL>`   | Specify a model name. See [DNA Models](DNA_MODELS) and [Protein Models](AA_MODELS) for the list of supported models. *DEFAULT: GTR for DNA and LG for Protein data* |
| `-st <SEQ_TYPE>`  | Specify a sequence type as either of `DNA` or `AA` for DNA or amino-acid sequences. *DEFAULT: auto-detect from the alignment or model* |
| `-format <FORMAT>`  | Set the alignment format as either of PHYLIP, FASTA, [MAPLE](MAPLE_FORMAT), or AUTO. *DEFAULT: auto-detect from the alignment* |
| `-t <TREE_FILE>`   | Specify a file containing a starting tree for tree search. Note: the starting tree is not mandatory to consist all taxa in the input alignment. |
| `-blfix`   | Keep the branch lengths of the input tree unchanged (only applicable if the input tree consists all the taxa in the alignment). |
| `-tree-search <TYPE>`   | Specify a [tree search type](TREE_SEARCH) as either of `FAST`, `NORMAL`, or `MORE_ACCURATE`. *DEFAULT: `NORMAL`* |
| `-shallow-search`   | Enable a shallow tree search before a deeper tree search. *DEFAULT: No shallow search* |
| `-branch-support`   | Compute branch supports (aLRT-SH) of the tree. |
| `--replicates <NUM>`   | Set the number of replicates for computing branch supports (aLRT-SH). *DEFAULT: 1000* |
| `-eps <NUM>`   | Set the epsilon value (see [Guindon et al., 2010]) for computing branch supports (aLRT-SH). *DEFAULT: 0.1*|
| `-nt <NUM_THREADS>` | Set the number of CPU cores used for computing branch supports. One can use `-nt AUTO` to use all CPU cores available on the current machine. *DEFAULT: 1* |
| `-pre <PREFIX>` | Specify a prefix for all output files. *DEFAULT: the alignment file name (`-aln`)* |
| `-rep-tree` | Allow CMAPLE to replace the input tree if a better likelihood tree is found when computing branch supports. |
| `-out-mul-tree` | Output the tree in multifurcating format. *DEFAULT: bifurcating tree* |
| `-overwrite` | Overwrite output files if existing. |
| `-ref <FILENAME>,<SEQNAME>` | Specify the reference genome by a sequence named `<SEQNAME>` from an alignment file `<FILENAME>`. |
| `-out-aln <NAME>,<FORMAT>` | Write the input alignment to a file named `<NAME>` in a specific format `<FORMAT>` which could be PHYLIP, FASTA, or [MAPLE](MAPLE_FORMAT). |
| `--min-blength <NUM>` | Set the minimum branch length. *DEFAULT: 0.2 x \<one mutation per site\>* |
| `--threshold-prob <NUM>` | Specify a relative probability threshold, which is used to ignore possible states with very low probabilities. *DEFAULT: 1e-8* |
| `-mut-update <NUM>` | Specify the period (in term of the number of sample placements) to update the substitution rate matrix. *DEFAULT: 25* |
| `-seed <NUM>` | Set a random number seed to reproduce a previous run. *DEFAULT: the CPU clock* |
| `-v <MODE>`   | Set the verbose mode (`QUIET`, `MIN`, `MED`,  `MAX`, or `DEBUG`) to control the amount of messages to screen, which is uselful for debugging purposes. *DEFAULT: `MED`* |

[Gaston et al. 2011]: https://doi.org/10.1093/bioinformatics/btr470
[MAPLE_FORMAT]: #
[DNA_MODELS]: #
[AA_MODELS]: #
[TREE_SEARCH]: #
[Guindon et al., 2010]: https://academic.oup.com/sysbio/article/59/3/307/1702850