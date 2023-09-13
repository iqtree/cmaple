
# What's CMAPLE?


Introduction to CMAPLE...


# How to cite CMAPLE?

[C]MAPLE paper...


# How to use CMAPLE?


## Installation

CMAPLE executables for different platforms such as Linux, OSX, and Windows are available at [CMAPLE-RELEASES]([CMAPLE_RELEASES]). 

## Usage examples

In the release package, we provide two excutables `cmaple` and `cmaple-aa` for analysing DNA and amino acid data, respectively. For simplicity, we use `cmaple` in the following examples. 

As a command-line program, CMAPLE can be run by executing `cmaple ...` from a terminal/console (or a command prompt under Windows). 

### 1. Infer a phylogenetic tree from an alignment

Together with the executables, there is an example alignment (`example.maple` in the `example` directory) in the release package. Once can reconstruct a phylogenetic tree from that alignment (assuming that you are now in the same folder with `example.maple`).

    cmaple -aln example.maple
    
In the above command,

* `-aln` is to specify an input alignment, which could be in FASTA, PHYLIP, or [MAPLE]([MAPLE_FORMAT]) format.

CMAPLE uses the default model, which is [General Time Reversible]([Tavare, 1986]) model for the DNA data in this example. At the end of the run, CMAPLE outputs the following files.

* `example.maple.treefile`: the inferred phylogenetic tree, which can be visualized by many tree viewer programs (e.g., FigTree).

* `example.maple.log`: log file captures all messages printed on the screen during the run. To report bugs, please send this log file and the input alignment (if possible) to the authors.


### 2. Specify a substitution model

CMAPLE supports various common DNA and empirical amino-acid models (replicated from the [IQ-TREE]([IQ_TREE]) software) as shown in [Supported substitution models](#supported-substitution-models). One can specify the [Jukes Cantor]([Jukes and Cantor, 1969]) model for the inference with the `-m` option.

    cmaple -aln example.maple -m JC

### 3. Specify an input tree
One can specify an input tree (e.g., `tree.nwk`) in a NEWICK format for the tree search using the `-t` option. 

	    cmaple -aln example.maple -t tree.nwk

**If the input tree is incomplete** (which doesn't contain all the taxa in the alignment), CMAPLE will:

* Performs placement (i.e., adding missing taxa from the alignment to the tree);
* Applies a NORMAL tree search (which does SPR moves only on newly-added nodes);
* Optimizes all branch lengths.

**If the input tree is complete** (which contains all the taxa in the alignment), CMAPLE will, by default, does neither placment nor tree search, but it optimizes all branch lengths.To keep the branch lengths fixed, one can add `-blfix` to the command.

        cmaple -aln example.maple -t tree.nwk -blfix

To use the input (complete/incomplete) tree as a starting tree to perform placement (for an incomplete tree), then consider SPR moves on all nodes, and optimize branch lengths, one can add `-tree-search MORE_ACCURATE` to the command.

        cmaple -aln example.maple -t tree.nwk -tree-search MORE_ACCURATE
        
### 4. Set the tree search type

We implemented three types of tree search:

 Tree search type        | Explanation |
|------------------------|---------------------------------------------------------------|
| FAST\_TREE\_SEARCH       | No tree search (placement only) |
| NORMAL\_TREE\_SEARCH     | Only consider pruning branches at newly-added nodes when seeking SPR moves |
| MORE\_ACCURATE\_TREE\_SEARCH| Consider all nodes when seeking SPR moves |

If users don't specify an input tree, `NORMAL_TREE_SEARCH` and `MORE_ACCURATE_TREE_SEARCH` perform the same behaviors since all taxa from the alignment are first added to the tree then SPR moves are considered at all those (newly-added) nodes during the tree search.

If users specify an input complete tree, both `FAST_TREE_SEARCH` and `NORMAL_TREE_SEARCH` do nothing since no new taxa is added to the tree.

If users specify an input incomplete tree, those tree search types have completely different behaviors as described in the above table. The runtime and the accuracy increase, in general, when changing the tree search type from the top to the bottom ones. 

By default, CMAPLE applies the `NORMAL_TREE_SEARCH`. One can change it by using the `-tree-search` option.

    cmaple -aln example.maple -t tree.nwk -tree-search FAST_TREE_SEARCH

    
### 5. Assess branch supports with aLRT-SH

CMAPLE implemented the SH-like approximate likelihood ratio test ([Guindon et al., 2010]). To perform this test, one can run:

    cmaple -aln example.maple -branch-support
    
To speed up the assessment, one can employ multithreading by adding `-nt` option

    cmaple -aln example.maple -branch-support -nt 8
    
In the above example, CMAPLE uses 8 threads for computing branch supports. One can use `-nt AUTO` to employ all CPU cores available on the current machine.

One can specify the number of replicates (default, 1000) and the epsilon value (default, 0.1) (see [Guindon et al., 2010]) by using `--replicates` and `-eps` options

    cmaple -aln example.maple -branch-support --replicates 5000 -eps 0.05

### 6. Convert an alignment to a different format

One can use the `-out-aln` option to write an input alignment to a file in a specific format such as PHYLIP, FASTA, or [MAPLE]([MAPLE_FORMAT]) format. For example, the following command

    cmaple -aln example.maple -out-aln aln.phy,PHYLIP
   
writes the input alignment to `aln.phy` file in a PHYLIP format.

    
## Supported substitution models

All the supported substitution models in CMAPLE are listed in the following.

**DNA models**

| Model        | Explanation |
|--------------|---------------------------------------------------------------|
| JC or JC69   | Equal substitution rates and equal base frequencies ([Jukes and Cantor, 1969]). |
| GTR          | General time reversible model with unequal rates and unequal base freq. ([Tavare, 1986]). |
| UNREST       | Unrestricted model with non-reversible, unequal rates and unequal base freq. |

**Amino-acid models**

| Model | Region | Explanation |
|-------|--------|---------------------------------------------------------------|
| Blosum62 | nuclear | BLOcks SUbstitution Matrix ([Henikoff and Henikoff, 1992]). Note that `BLOSUM62` is not recommended for phylogenetic analysis as it was designed mainly for sequence alignments. |
| cpREV    | chloroplast |chloroplast matrix ([Adachi et al., 2000]). |
| Dayhoff  | nuclear | General matrix ([Dayhoff et al., 1978]). |
| DCMut    | nuclear | Revised `Dayhoff` matrix ([Kosiol and Goldman, 2005]). |
| FLAVI    | viral | Flavivirus ([Le and Vinh, 2020]). | 
| FLU      | viral | Influenza virus ([Dang et al., 2010]). |
| GTR20    | general | General time reversible models with 190 rate parameters. *WARNING: Be careful when using this parameter-rich model as parameter estimates might not be stable, especially when not having enough phylogenetic information (e.g. not long enough alignments).* |
| HIVb     | viral | HIV between-patient matrix HIV-B<sub>m</sub> ([Nickle et al., 2007]). |
| HIVw     | viral | HIV within-patient matrix HIV-W<sub>m</sub> ([Nickle et al., 2007]). |
| JTT      | nuclear | General matrix ([Jones et al., 1992]). |
| JTTDCMut | nuclear | Revised `JTT` matrix ([Kosiol and Goldman, 2005]). |
| LG       | nuclear | General matrix ([Le and Gascuel, 2008]). |
| mtART    | mitochondrial | Mitochondrial Arthropoda ([Abascal et al., 2007]). |
| mtMAM    | mitochondrial | Mitochondrial Mammalia ([Yang et al., 1998]). |
| mtREV    | mitochondrial | Mitochondrial Vertebrate ([Adachi and Hasegawa, 1996]). |
| mtZOA    | mitochondrial | Mitochondrial Metazoa (Animals) ([Rota-Stabelli et al., 2009]). |
| mtMet    | mitochondrial | Mitochondrial Metazoa ([Vinh et al., 2017]). |
| mtVer    | mitochondrial | Mitochondrial Vertebrate ([Vinh et al., 2017]). |
| mtInv    | mitochondrial | Mitochondrial Inverterbrate ([Vinh et al., 2017]). |
| NQ.bird   | nuclear | Non-reversible Q matrix ([Dang et al., 2022]) estimated for birds ([Jarvis et al., 2015]). | 
| NQ.insect | nuclear | Non-reversible Q matrix ([Dang et al., 2022]) estimated for insects ([Misof et al., 2014]). | 
| NQ.mammal | nuclear | Non-reversible Q matrix ([Dang et al., 2022]) estimated for mammals ([Wu et al., 2018]). | 
| NQ.pfam   | nuclear | General non-reversible Q matrix ([Dang et al., 2022]) estimated from Pfam version 31 database ([El-Gebali et al., 2018]). | 
| NQ.plant  | nuclear | Non-reversible Q matrix ([Dang et al., 2022]) estimated for plants ([Ran et al., 2018]). | 
| NQ.yeast  | nuclear | Non-reversible Q matrix ([Dang et al., 2022]) estimated for yeasts ([Shen et al., 2018]). | 
| Poisson  | none | Equal amino-acid exchange rates and frequencies. |
| PMB      | nuclear | Probability Matrix from Blocks, revised `BLOSUM` matrix ([Veerassamy et al., 2004]). |
| Q.bird   | nuclear | Q matrix ([Minh et al., 2021]) estimated for birds ([Jarvis et al., 2015]). | 
| Q.insect | nuclear | Q matrix ([Minh et al., 2021]) estimated for insects ([Misof et al., 2014]). | 
| Q.mammal | nuclear | Q matrix ([Minh et al., 2021]) estimated for mammals ([Wu et al., 2018]). | 
| Q.pfam   | nuclear | General Q matrix ([Minh et al., 2021]) estimated from Pfam version 31 database ([El-Gebali et al., 2018]). | 
| Q.plant  | nuclear | Q matrix ([Minh et al., 2021]) estimated for plants ([Ran et al., 2018]). | 
| Q.yeast  | nuclear | Q matrix ([Minh et al., 2021]) estimated for yeasts ([Shen et al., 2018]). | 
| rtREV    | viral | Retrovirus ([Dimmic et al., 2002]). |
| VT       | nuclear | General 'Variable Time' matrix ([Mueller and Vingron, 2000]). |
| WAG      | nuclear | General matrix ([Whelan and Goldman, 2001]). |








# Command reference

All the options available in CMAPLE are shown below:

| Option | Usage and meaning |
|--------------|------------------------------------------------------------------------------|
|`-h` or `-?`| Print help usage. |
|`-aln <ALIGNMENT>`| Specify an input alignment file in PHYLIP, FASTA, or [MAPLE]([MAPLE_FORMAT]) format. |
| `-m <MODEL>`   | Specify a model name. See [DNA Models]([DNA_MODELS]) and [Protein Models]([AA_MODELS]) for the list of supported models. *DEFAULT: GTR for DNA and LG for Protein data* |
| `-st <SEQ_TYPE>`  | Specify a sequence type as either of `DNA` or `AA` for DNA or amino-acid sequences. *DEFAULT: auto-detect from the alignment or model* |
| `-format <FORMAT>`  | Set the alignment format as either of PHYLIP, FASTA, [MAPLE]([MAPLE_FORMAT]), or AUTO. *DEFAULT: auto-detect from the alignment* |
| `-t <TREE_FILE>`   | Specify a file containing a starting tree for tree search. Note: the starting tree is not mandatory to consist all taxa in the input alignment. |
| `-blfix`   | Keep the branch lengths of the input tree unchanged (only applicable if the input tree consists all the taxa in the alignment). |
| `-tree-search <TYPE>`   | Specify a [tree search type]([TREE_SEARCH]) as either of `FAST`, `NORMAL`, or `MORE_ACCURATE`. *DEFAULT: `NORMAL`* |
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
| `-out-aln <NAME>,<FORMAT>` | Write the input alignment to a file named `<NAME>` in a specific format `<FORMAT>` which could be PHYLIP, FASTA, or [MAPLE]([MAPLE_FORMAT]). |
| `--min-blength <NUM>` | Set the minimum branch length. *DEFAULT: 0.2 x \<one mutation per site\>* |
| `--threshold-prob <NUM>` | Specify a relative probability threshold, which is used to ignore possible states with very low probabilities. *DEFAULT: 1e-8* |
| `-mut-update <NUM>` | Specify the period (in term of the number of sample placements) to update the substitution rate matrix. *DEFAULT: 25* |
| `-seed <NUM>` | Set a random number seed to reproduce a previous run. *DEFAULT: the CPU clock* |
| `-v <MODE>`   | Set the verbose mode (`QUIET`, `MIN`, `MED`,  `MAX`, or `DEBUG`) to control the amount of messages to screen, which is uselful for debugging purposes. *DEFAULT: `MED`* |

[CMAPLE_RELEASES]: https://github.com/trongnhanuit/cmaple/releases
[Gaston et al. 2011]: https://doi.org/10.1093/bioinformatics/btr470
[MAPLE_FORMAT]: #
[DNA_MODELS]: #
[AA_MODELS]: #
[TREE_SEARCH]: #
[Guindon et al., 2010]: https://academic.oup.com/sysbio/article/59/3/307/1702850
[IQ_TREE]: http://www.iqtree.org/

[Abascal et al., 2007]: https://doi.org/10.1093/molbev/msl136
[Adachi and Hasegawa, 1996]: https://doi.org/10.1007/BF02498640
[Adachi et al., 2000]: https://doi.org/10.1007/s002399910038
[Bielawski and Gold, 2002]: https://doi.org/10.1093/genetics/161.4.1589
[Dang et al., 2010]: https://doi.org/10.1186/1471-2148-10-99
[Dang et al., 2022]: https://doi.org/10.1093/sysbio/syac007
[Dayhoff et al., 1978]: http://compbio.berkeley.edu/class/c246/Reading/dayhoff-1978-apss.pdf
[Dimmic et al., 2002]: https://doi.org/10.1007/s00239-001-2304-y
[El-Gebali et al., 2018]: https://doi.org/10.1093/nar/gky995
[Felsenstein, 1981]: https://doi.org/10.1007%2FBF01734359
[Goldman and Yang, 1994]: http://mbe.oxfordjournals.org/content/11/5/725.abstract
[Gu et al., 1995]: http://mbe.oxfordjournals.org/content/12/4/546.abstract
[Hasegawa, Kishino and Yano, 1985]: https://dx.doi.org/10.1007%2FBF02101694
[Henikoff and Henikoff, 1992]: https://dx.doi.org/10.1073%2Fpnas.89.22.10915
[Jarvis et al., 2015]: https://doi.org/10.1186/s13742-014-0038-1
[Jones et al., 1992]: https://dx.doi.org/10.1093%2Fbioinformatics%2F8.3.275
[Jukes and Cantor, 1969]: http://doi.org/10.1016/B978-1-4832-3211-9.50009-7
[Kimura, 1980]: https://doi.org/10.1007%2FBF01731581
[Kimura, 1981]: https://doi.org/10.1073/pnas.78.1.454
[Kosiol and Goldman, 2005]: https://doi.org/10.1093/molbev/msi005
[Kosiol et al., 2007]: https://doi.org/10.1093/molbev/msm064
[Lartillot and Philippe, 2004]: https://doi.org/10.1093/molbev/msh112
[Le and Gascuel, 2008]: https://doi.org/10.1093/molbev/msn067
[Le and Vinh, 2020]: https://doi.org/10.1007/s00239-020-09943-3
[Le et al., 2008a]: https://doi.org/10.1093/bioinformatics/btn445
[Le et al., 2008b]: https://doi.org/10.1098/rstb.2008.0180
[Le and Gascuel, 2010]: https://doi.org/10.1093/sysbio/syq002
[Le et al., 2012]: https://doi.org/10.1093/molbev/mss112
[Lewis, 2001]: https://doi.org/10.1080/106351501753462876
[Minh et al., 2021]: https://doi.org/10.1093/sysbio/syab010
[Misof et al., 2014]: https://doi.org/10.1126/science.1257570
[Mueller and Vingron, 2000]: https://doi.org/10.1089/10665270050514918
[Muse and Gaut, 1994]: http://mbe.oxfordjournals.org/content/11/5/715.abstract
[Nickle et al., 2007]: https://dx.doi.org/10.1371/journal.pone.0000503
[Ran et al., 2018]: https://doi.org/10.1098/rspb.2018.1012
[Rota-Stabelli et al., 2009]: https://doi.org/10.1016/j.ympev.2009.01.011
[Schneider et al., 2005]: https://doi.org/10.1186/1471-2105-6-134
[Shen et al., 2018]: https://doi.org/10.1016/j.cell.2018.10.023
[Soubrier et al., 2012]: https://doi.org/10.1093/molbev/mss140
[Tamura and Nei, 1993]: http://mbe.oxfordjournals.org/cgi/content/abstract/10/3/512
[Tavare, 1986]: http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf
[Veerassamy et al., 2004]: https://doi.org/10.1089/106652703322756195
[Vinh et al., 2017]: https://doi.org/10.1186/s12862-017-0987-y
[Wang et al., 2008]: https://doi.org/10.1186/1471-2148-8-331
[Whelan and Goldman, 2001]: https://doi.org/10.1093/oxfordjournals.molbev.a003851
[Woodhams et al., 2015]: https://doi.org/10.1093/sysbio/syv021
[Wu et al., 2018]: https://doi.org/10.1016%2Fj.dib.2018.04.094
[Yang, 1994]: https://doi.org/10.1007/BF00160154
[Yang, 1995]: http://www.genetics.org/content/139/2/993.abstract
[Yang et al., 1998]: http://mbe.oxfordjournals.org/content/15/12/1600.abstract
[Zharkikh, 1994]: https://doi.org/10.1007/BF00160155
[Guindon et al., 2010]: https://doi.org/10.1093/sysbio/syq010