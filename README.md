# FishLifeExonCapture

For more information check the following repo: [FishLifeExonCapture](https://github.com/lilychughes/FishLifeExonCapture). This is just a python version and easy-install initiative of the original repo

### Requirement:

* [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/download/)
* [Git](https://git-scm.com/downloads/)

### Installation:

```shell
git clone https://github.com/Ulises-Rosas/FishLifeExonCapture.git && cd FishLifeExonCapture 
conda env create -f environment.yml
source activate flec  
```

# Commands available

Each command has its own helping page (either with `-h` or `--help`) if more information about available options is needed. All commands assume your are working in a same working directory (i.e. directory where your fastq files are initially located), but this can also be modified (see option '-p' on each command).  Everything between brackets (`[]`) should be fulfilled by the user accordingly


### Step 0: Pack fastq files into directories

Given these following files in your working directory:
```
Family_Genus_species_S1234_R1.fastq.gz
Family_Genus_species_S1234_R2.fastq.gz
```
The directory `Family_Genus_species_S1234` is created storing above files with:
```
fishmanager mkdir
```

Note that your working directory may contain more than two files. Check a progress report with: `fishmanager look`.

### Step 1: Run Trimmomatic to quality trim the sequences
```
trimmomatic-loop-PE -a [adapter file] -n [number of cpus]
```

### Step 2: Map raw reads back to representative bait sequences

```
map-exons -n [number of cpus]
```

For the older set of Otophysi markers (see Arcila et al 2017), use:
```
map-exons-otophysi -n [number of cpus]
```

### Step 3: Build initial assemblies in Velvet

```
initialVelvet -n [number of cpus]
```

### Step 4: Run aTRAM

```
runaTRAM -n [number of cpus]
```

To deal with permission issues at clusters, the option `-r` might be useful. Check `runaTRAM -h`


### Step 5: Find reading frames and filter exons

For percomorph fishes:
```
ExonFiltering Percomorph -n [number of cpus]
```

For elopomorph fishes:
```
ExonFiltering Elopomorph -n [number of cpus]
```

For osteoglossomorph fishes:
```
ExonFiltering Osteoglossomorph -n [number of cpus]
```
Otophysi set of markers available (see Arcila et al. 2017). It just calls a different set of markers and reading frames.

```
ExonFiltering Otophysi -n [number of cpus]
```

### Step 5b: Filter Exons and Flanking Introns (Optional)

```
FlankFiltering -n [number of cpus]
```


### Step 6: Reading-Frame Aware Alignment

```
preAlignment -n [number of cpus]
```

An alternative version exists for the older set of Otophysan markers:

```
preAlignment_Otophysi -n [number of cpus]
```

To run MACSE2:

```
run_macse -n [number of cpus] -p Alignments
```
Where `Alignments` is the directory previously created by the `preAlignment` command. 
As with the last step, there is an alternative script for the otophysan markers:
```
run_macse_Otophysi -n [number of cpus] -p Alignments
```

### Step 6b: Alignment for contigs with flanking regions
This will create a folder called Alignments_Flanks

```
 preAlignmentFlanks -n [number of cpus]
```
alternative for the otophysan markers:
```
preAlignmentFlanks_Otophysi -n [number of cpus]
```

### Fishmanager

This is a submodule that allows to look/modify metadata. If the metadata is split into a given number of branches (with `fishmanager split`), each branch can be run in parallel (i.e. submitting multiples job files) with the option `-b` available in almost all commands. This multibatch approach speeds up analyses
