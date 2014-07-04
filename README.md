TEAM: TE Annotation from Methylation
====================================

Transposable elements (TEs) are DNA sequences that can "jump" and replicate throughout their host genome. Because their repetitive and prolific presence increases the difficulty of genome assembly, sequence alignment, and genome annotation, new algorithms are continuously developed. The detection and classification of transposable elements is crucial since they comprise significant portions of eukaryotic genomes and may induce large-scale genome rearrangement. Currently, transposable elements are identified through homology search or de novo repeat discovery. Homology searching relies on previously characterized transposon families and is not appropriate for new genomes. While de novo repeat discovery methods can detect highly repeated novel transposable elements, they may report non-TE repeats, such as tandem repeats, segmental duplications, and satellites. There are also low copy number transposons that are kept silent through DNA methylation, which are difficult to detect through existing methods. To improve the detection of low copy number transposable elements, I propose TEAM, which detects TEs in a reference genome based on its methylation signature. TEAM scans the frequencies of each methylation motif (CG, CHH, and CHG) in a sliding window across the whole genome and detects the unique methylation profiles of TEs, pseudogenes, and protein-coding genes using a hidden markov model. Not only is TEAM more sensitive than existing algorithms, but it also demands less memory and processing time.

Dependencies
------------

- NumPy
- Matplotlib
- Python C Headers

Installation
------------

```Shell
python setup.py install
```

or, for local installations

```Shell
python setup.py install --user
```

Example
-------------

First, TEAM needs to be trained on a set of data. This is done by calculating the model emissions and transitions from a `methExperiment` of methylation and feature files.

### Making methExperiment

```python
from team import methExperiment
from team import annotate

testExp = methExperiment()
testExp.addSample("30-109",["30-109_1_methratio.txt","30-109_2_methratio.txt"])
testExp.addGFF("NC","region_gffs/TAIR10_NC.gff") #non-coding regions
testExp.addGFF("TG","region_gffs/TAIR10_te_gene.gff") #TE-genes
testExp.addGFF("G","region_gffs/TAIR10_genes.gff") #genes
testExp.addGFF("PG","region_gffs/TAIR10_pseudogene.gff") #pseudogenes
testExp.addGFF("TE","region_gffs/TAIR10_te.gff") #TEs
testExp.setReference("TAIR10.fa")
testExp.printExperiment()
```
### Training Model

```python
testExp.makeMethBins()
testExp.makeEmissions()
#testExp.readEmissions() #read previously generated emissions
testExp.printEmissions()
testExp.makeTransitions()
```

Lastly, TEs can be discovered in new data using the generated model.

### Annotate Data

```python
testExp = methExperiment()
testExp.addSample("30-29",["30-29_1_methratio.txt","30-29_2_methratio.txt","30-29_3_methratio.txt"])
annotate.regionFinder(trainExp.emissions,trainExp.transMatrix,trainExp.fa,trainExp.gff,testExp.samples, 500, 50)
```

A more detailed example with test data will be released soon.
