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

```python
from team import methExperiment
import sys

testExp = methExperiment()
testExp.addSample("30-109",["30-109_1_methratio.txt","30-109_2_methratio.txt"])
testExp.addGFF("NC","region_gffs/TAIR10_NC.gff")
testExp.addGFF("TE genes","region_gffs/TAIR10_te_gene.gff")
testExp.addGFF("Genes","region_gffs/TAIR10_genes.gff")
testExp.addGFF("Pseudogenes","region_gffs/TAIR10_pseudogene.gff")
testExp.addGFF("TEs","region_gffs/TAIR10_te.gff")
testExp.setReference("TAIR10.fa")
testExp.printExperiment()
##############################
# Make emission probabilities
##############################
testExp.makeMethBins()
testExp.printProbs()
testExp.plotData()
testExp.writeProbs()
##############################
# Load emission probabilities
##############################
#testExp.readProbs()
sys.exit()
```
A more detailed example with test data will be released soon.
