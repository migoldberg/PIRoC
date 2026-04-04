# PIRoC: Phylogeny-Informed Removal of Contaminants
By: Mark I. Goldberg, Nickellaus G. Roberts, Kevin M. Kocot

Kocot Lab\
Department of Biological Sciences\
The University of Alabama\
Tuscaloosa, Alabama, USA

---

### Introduction

PIRoC is a tool that uses single-gene trees to screen for exogenous contamination in transcriptomes. Orthologs are identified among taxa of interest plus "lightning rod" taxa representing likely contaminants (e.g., various algae if studying an algae-eater) and gene trees are inferred. PIRoC then evaluates all sequences originating from a defined focal group, and outputs decontaminated sequences.

\
A PIRoC run requires three inputs to start:
1. A directory containing **single-gene trees** and corresponding **sequence alignments**.
2. A metadata file with **group assignments** (e.g., ingroup taxa [Mollusca], ignored outgroup taxa [Outgroup], and ‘lightning rod’ taxa closely related to likely contaminants [Contaminant]) for each taxon.
3. A **focal group** (the group of taxa that will be decontaminated).

\
For the metadata file, each taxon in the dataset should be assigned a "group". This group can be as specific or broad as desired. A broader set of groups (i.e., phylum) would result in stricter conditions for a sequence to be classified as clean. A more specific group assignment (i.e., family) would result in less strict conditions for a clean sequence classification. The species ID should correspond to the taxon name in the gene-tree and sequence alignment headers.

**Example Metadata Files**\
*Example 1*
```
species_id                group
Sinonovacula_constricta   Mollusca
Crassostrea_gigas         Mollusca
Architeuthis_dux          Mollusca
Haliotis_rufescens        Mollusca

Lingula_anatina           Outgroup
Eisenia_andrei            Outgroup

Gracilaria_domingensis    Contaminant
Cryptomonas_paramaecium   Contaminant
Andalucia_godoyi          Contaminant
Phaeodactylum_tricornutum Contaminant
Polarella_glacialis       Contaminant
```
*Example 2*
```
species_id                group
Sinonovacula_constricta   Bivalvia
Crassostrea_gigas         Bivalvia
Architeuthis_dux          Cephalopada
Haliotis_rufescens        Gastropoda

Lingula_anatina           Outgroup
Eisenia_andrei            Outgroup

Gracilaria_domingensis    Contaminant
Cryptomonas_paramaecium   Contaminant
Andalucia_godoyi          Contaminant
Phaeodactylum_tricornutum Contaminant
Polarella_glacialis       Contaminant
```

---

### Installation
PIRoC is still in development and has therefore not yet been published on PyPI for installation via pip. It can be installed by cloning the source code from this GitHub repository. 

```
git clone https://github.com/migoldberg/PIRoC
cd PIRoC
pip install -e .
```
Required packages include Numpy and BioPython. PIRoC was written in Python 3.11. These required packages should automatically be installed when the command above is run.

---

### Example Run
Here is an example command from running PIRoC. These commands were used for evaluation of these methods with a dataset of 22 molluscs, two related outgroup taxa (one annelid and one brachiopod), and 17 diverse bait eukaryotes representing likely contamination in filter-feeding bivalves. 

\
`Command`
```
PIRoC  \
    --tree_dir  $TREE_DIR  \
    --focal_group  Mollusca  \
    --metadata  metadata.tsv  \
    --remove_contaminants
```

\
`metadata.tsv`
```
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

\
Prior to running PIRoC for this dataset, OrthoFinder (2.5.4; Emms and Kelly 2019) was used for the inference of orthogroups. To simulate a standard workflow, an edited version of the pipeline of Kocot et al. (2017) was used for quality control, multiple sequence alignment and tree-building.

---

### Arguments/CLI Options

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
| `min_clean_purity`      | `--min_clean_purity`          | The minimum fraction of sequences that must belong to the focal group in a clade for a CLEAN classification (upper bound of `max_contaminant_purity`)            | `0.8`          |
| `max_contaminant_purity` | `--max_contaminant_purity`     | The highest fraction of sequences that can belong to the focal group in a clade for a CONTAMINANT classification (lower bound of `min_target_purity`)        | `0.5`          |
| `review_flags`           | `--review-flags`               | Interactively review flagged sequences via CLI (in development)                     | `False`        |
| `loud`                   | `--loud`                       | Enable detailed console output                        | `False`        |



---

### Notes
For help, suggestions, or other feedback, please contact migoldberg@crimson.ua.edu.

This project will be presented at the Undergraduate Research and Creative Activities (URCA) conference in April of 2026 at the University of Alabama.

---

### References
1. Kocot, K. M., Struck, T. H., Merkel, J., Waits, D. S., Todt, C., Brannock, P. M., Weese, D. A., Cannon, J. T., Moroz, L. L., Lieb, B., & Halanych, K. M. (2017). **Phylogenomics of Lophotrochozoa with Consideration of Systematic Error.** *Systematic Biology*, 66(2), 256–282.   https://doi.org/10.1093/sysbio/syw079

2. Emms, D. M., & Kelly, S. (2019). **OrthoFinder: Phylogenetic orthology inference for comparative genomics.** *Genome Biology*, 20, 238. https://doi.org/10.1186/s13059-019-1832-y
