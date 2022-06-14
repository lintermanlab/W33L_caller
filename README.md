# W33L_caller
Call W33L from Sanger sequencing of the VH186.2 locus.

To support Silva-Cayetano _et al._ submitted. The isolation of single B cells, and the generation of PCR products for sequencing are all described in the manuscript.

## Outline

An Rscript that automates sequence trace analysis.


1. sequence traces were read using readsangerseq function with default parameters from the sangerseqR package (https://www.bioconductor.org/packages/release/bioc/html/sangerseqR.html), published here: [Hill _et al._]( https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/dvdy.24183).

2. Simple quality control was applied: sequences shorter than 300 nucleotides, or containing more than one N base call were removed (these are unaligned sequences, so N reflects poor sequencing traces rather than SHM).

3. The W33L locus was identified using matchPattern function - allowing 2 nucleotide mismatches and indels) from the Biostrings package (https://bioconductor.org/packages/Biostrings) - searching for `ACCAGCTACTNNATGCACTGG` in the reverse complemented sequence data.
  
4. a `W` or `L` call was assigned as follows:
- if `TNN` (in the appropriate position in the nucleotide string above) was `TGG` the assignment was `W`
- if `TNN` was `TTA` or `TTG`, this was assigned `L`
- any alternative sequences for `TNN` were assigned `other`

5. Per sample calls were exported as a CSV for downstream analysis.

## Use

The user provides a directory of sequencing trace files (as ab1 files), and a results directory.
Particular user attention is required at the end of the script: the `tabulate the results` section uses ab1 filenames to assign genotypes and replicates, and these will vary with differ labs + setup. It is extremely unlikely that the out-of-the-box code would assign genotypes etc correctly.

