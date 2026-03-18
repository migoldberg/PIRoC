# PIRoC: Phylogeny-Informed Removal of Contaminants
By Mark I. Goldberg, Nickellaus G. Roberts, Kevin M. Kocot

Kocot Lab  
Department of Biological Sciences  
The University of Alabama  
Tuscaloosa, Alabama, USA  

***
### Introduction
A tool that uses phylogenomic trees to classify and remove individual transcript-derives sequences as potentail contaminants. PIRoC is intended to be used after OrthoFinder and standard tree-building workflows. The tool accepts these orthogroups (in fasta format), the constructed trees, and a taxa metadata file. In this metadata file, you should assign each taxa a group. This group can be as specific or broad as desired, based on your dataset. A broader group assignment (ie. phylum) would be more conservative, while a more specific group assignment (ie. family) would result in stricter requirements for a sequence to be kept. The final required input is a "focal group", which is the group you seek to decontaminate. For example, if you include a large number of molluscs to decontaminate, and assign them all the group "Mollusca" in the metadata, you would pass "Mollusca" as your focal group. PIRoC will index all sequences in the phylogenomic trees which belong to this assigned group, and perform the classifications. The tool is designed to recognize common contaminants, which should be included in the dataset, as certain groups. By default, this group assignment is "Contaminant" but this can be changed by the user to another name or multiple names.  
***

### Standard Workflow
For a PIRoC run, the expected workflow would be the following:
1. Preperation of a transcriptome dataset.
2. Identification and inclusion of common contaminants for your focal group.
3. OrthoFinder
4. Standard quality control, MSA, and tree building workflows.
5. Assignment of groups in the metadata file, including your contaminants as "Contaminant".
6. Run PIRoC


***


### Instillation (Development Build)

    git clone https://github.com/migoldberg/PIRoC
    cd PIRoC
    pip install -e .

Required Packages include Numpy and BioPython. The tool was written in Python 3.11.

***

### Example
Here is an example command for running PIRoC using a dataset of 22 molluscs, 2 related taxa, and 17 common contaminants in molluscs. 

`Command:`

    PIRoC  \
	    --tree_dir  $TREE_DIR  \
	    --focal_group  Mollusca  \
	    --metadata  metadata.tsv  \
	    --remove_contaminants

  
  `metadata.tsv`

  ```tsv
species_id	group
adux	Mollusca
npom	Mollusca
obim	Mollusca
esco	Mollusca
aimm	Mollusca
echl	Mollusca
pcan	Mollusca
lgig	Mollusca
csqu	Mollusca
hruf	Mollusca
bgla	Mollusca
scon	Mollusca
cgig	Mollusca
pfuc	Mollusca
cfar	Mollusca
sbro	Mollusca
mmer	Mollusca
mphi	Mollusca
pstr	Mollusca
agra	Mollusca
pver	Mollusca
sdal	Mollusca

lana	Outgroup
eand	Outgroup

Gracilaria_domingensis_GCA_022539475.1	Contaminant
Chlamydomonas_reinhardtii_GCF_000002595.2	Contaminant
Cryptomonas_paramaecium_GCF_000194455.1	Contaminant
Andalucia_godoyi_GCA_009859145.1	Contaminant
Paulinella_micropora_GCA_009731375.1	Contaminant
Cafeteria_roenbergensis_GCA_008330645.1	Contaminant
Candida_tropicalis_GCF_000006335.3	Contaminant
Capsaspora_owczarzaki_GCF_001186125.2	Contaminant
Entamoeba_marina_GCA_051620835.1	Contaminant
Naegleria_fowleri_GCF_008403515.1	Contaminant
Phaeodactylum_tricornutum_GCF_000150955.2	Contaminant
Polarella_glacialis_GCA_905237085.1	Contaminant
Emiliana_huxleyi_GCF_000372725.1	Contaminant
Pelagomonas_calceolata_GCA_918797485.1	Contaminant
Sphaeroforma_arctica_GCF_001186125.1	Contaminant
Thecamonas_trahens_GCF_000142905.1	Contaminant
Paramecium_tetraurelia_GCF_000165425.1	Contaminant
```

***
### Arguments

| Argument                 | Flag(s)                        | Description                                                        | Default        |
| ------------------------ | ------------------------------ | ------------------------------------------------------------------ | -------------- |
| `tree_dir`               | `-t`, `--tree_dir`             | Directory containing Newick tree files and orthogroups in fasta format                              | (required)     |
| `suffix`                 | `-s`, `--suffix`               | Tree filename suffix                                               | `.tre`         |
| `output_dir`             | `-o`, `--output_dir`           | Output directory                                                   | `PIRoC_output` |
| `focal_group`            | `-f`, `--focal_group`          | Name of the focal group (e.g., Mollusca)                           | (required)     |
| `metadata`               | `-m`, `--metadata`             | TSV file with species_id and group columns                         | (required)     |
| `contaminants`           | `--contaminants`               | Comma-separated list of known contaminant group names              | `Contaminant`  |
| `remove_contaminants`    | `-rm`, `--remove_contaminants` | Remove contaminant sequences from orthogroups after classification | `False`        |
| `min_support`            | `--min_support`                | Minimum bootstrap support for confident classification             | `70.0`         |
| `min_target_purity`      | `--min_target_purity`          | The minimum fraction of sequences that must belong to the focal group in a clade for a TARGET classification (upper bound of `max_contaminant_purity`)            | `0.8`          |
| `max_contaminant_purity` | `--max_contaminant_purity`     | The highest fraction of sequences that can belong to the focal group in a clade for a CONTAMINANT classification (lower bound of `min_target_purity`)        | `0.5`          |
| `review_flags`           | `--review-flags`               | Interactively review flagged sequences via CLI (in development)                     | `False`        |
| `quiet`                  | `--quiet`                      | Enable quiet mode (minimal console output)                         | `False`        |

***
### Notes
For help, suggestions, or other feedback please contact migoldberg@crimson.ua.edu. 

Thank you to Dr. Kevin Kocot, Dr. Nickellaus Roberts, and the Kocot Lab for your contributions to and support of this project.
