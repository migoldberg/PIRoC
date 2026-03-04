# PIRoC: Phylogeny-Informed Removal of Contaminants
By Mark I. Goldberg, Nickellaus G. Roberts, Kevin M. Kocot

Kocot Lab
Department of Biological Sciences 
The University of Alabama
Tuscaloosa, Alabama, USA
***
### Abstract
Phylogenomics uses transcriptome-scale datasets to infer evolutionary relationships, particularly for non-model invertebrates. Despite recent advances, contamination remains a significant challenge for transcriptomes of small-bodied taxa, as whole organisms are often used for RNA extraction, allowing biological contaminants such as gut contents, symbionts, or epibionts to be co-extracted. Additional technical contamination may occur during library preparation or sequencing. If left unaddressed, such contamination can compromise the accuracy of downstream analyses and result in strongly supported but erroneous conclusions. To address this issue, we present PIRoC (​​Phylogeny-Informed Removal of Contaminants), a tool that uses phylogenomic trees to classify individual transcript-derived sequences as potential contaminants. PIRoC was evaluated using 42 transcriptomes, including 22 molluscs, 2 related invertebrates, and 18 phylogenomically diverse eukaryotic “lightning rod” taxa. A bivalve transcriptome was artificially spiked with 500 sequences from a diatom. Of these, 42 were missed by routine quality control measures. PIRoC correctly identified all 42 sequences as contaminants, while over 90% of bivalve sequences were retained.
***
### Instillation (Development Build)

    git clone https://github.com/migoldberg/PIRoC
    cd PIRoC
    pip install -e .
***
### Example
This is the example command/files for the run discussed in the abstract above. 

`Command:`

    PIRoC  \
	    --tree_dir  $TREE_DIR  \
	    --focal_species  scon  \
	    --metadata  species_metadata.tsv  \
	    --output_dir  $./PIRoC_output  \
	    --remove_contaminants

  
  `species_metadata.tsv`

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
Gracilaria_domingensis_GCA_022539475.1	Outgroup
Chlamydomonas_reinhardtii_GCF_000002595.2	Outgroup
Cryptomonas_paramaecium_GCF_000194455.1	Outgroup
Andalucia_godoyi_GCA_009859145.1	Outgroup
Paulinella_micropora_GCA_009731375.1	Outgroup
Loxomitra_sp_SRR11344291	Outgroup
Cafeteria_roenbergensis_GCA_008330645.1	Outgroup
Candida_tropicalis_GCF_000006335.3	Outgroup
Capsaspora_owczarzaki_GCF_001186125.2	Outgroup
Entamoeba_marina_GCA_051620835.1	Outgroup
Naegleria_fowleri_GCF_008403515.1	Outgroup
Phaeodactylum_tricornutum_GCF_000150955.2	Outgroup
Polarella_glacialis_GCA_905237085.1	Outgroup
Emiliana_huxleyi_GCF_000372725.1	Outgroup
Pelagomonas_calceolata_GCA_918797485.1	Outgroup
Sphaeroforma_arctica_GCF_001186125.1	Outgroup
Thecamonas_trahens_GCF_000142905.1	Outgroup
Paramecium_tetraurelia_GCF_000165425.1	Outgroup
```
The group column can be as specific or as broad as needed to fit your dataset. Larger and more diverse datasets should be more specific, while smaller and less diverse should be broader (ie. phylum). 
***
### Arguments

| Argument | Flag(s) | Description | Default |
|----------|---------|-------------|---------|
| `tree_dir` | `-t`, `--tree_dir` | Directory containing Newick tree files | (required) |
| `suffix` | `-s`, `--suffix` | Tree filename suffix | `.tre` |
| `output_dir` | `-o`, `--output_dir` | Output directory | `PIRoC_output` |
| `focal_species` | `-f`, `--focal_species` | Species ID of the focal species (e.g., obim) | (required) |
| `metadata` | `-m`, `--metadata` | TSV file with species_id and group columns | (required) |
| `min_support` | `--min_support` | Minimum bootstrap support for confident classification | `70.0` |
| `min_target_purity` | `--min_target_purity` | Minimum focal-group purity for TARGET | `0.8` |
| `max_contaminant_purity` | `--max_contaminant_purity` | Maximum focal-group purity for CONTAMINANT | `0.5` |
| `collapse_threshold` | `--collapse_threshold` | Collapse nodes with bootstrap below this value | `50.0` |
| `outgroups` | `--outgroups` | Comma-separated list of outgroup names | `Outgroup` |
| `remove_contaminants` | `-rm`, `--remove_contaminants` | Produce new tree files with contaminants removed | `False` |
| `debug` | `--debug` | Enable debug mode | `False` |

***
### Notes
For help, suggestions, or other feedback please contact migoldberg@crimson.ua.edu. 

Thank you to Dr. Kevin Kocot, Dr. Nickellaus Roberts, and the Kocot Lab for your contributions to and support of this project.
