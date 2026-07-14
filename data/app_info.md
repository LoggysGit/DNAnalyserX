# GeneAnalyserX BETA - Quick Guide
GeneAnalyserX is a lightweight tool that compares a patient's DNA sequence for a single gene against a healthy reference from NCBI, and reports exactly what differs. It runs mostly offline: reference data is downloaded once per gene and cached locally, so the tool works on any laptop **without** needing cloud infrastructure.
App uses Biopython and CustomTkinter libraries.

**YOU MUST SET A REAL EMAIL ADDRESS IN CONFIG.CFG BEFORE USE**

NCBI requires a valid email to identify requests to their servers. A placeholder or fake address can get your requests blocked.

---

## What You Need to Load
* **Patient File (`.fasta.gz`):** A DNA sequence file containing exactly **one** gene. The file must not contain more than one sequence, and it should not be a whole chromosome — this tool analyzes one gene at a time.
* **Gene Name:** The official gene symbol you want to check (for example: **TP53** or **BRCA1**).

---

## How It Works
* **Step 1: Automatic Download**: You type in the gene name, and the app fetches the matching reference sequence and its annotation from NCBI. Both files are cached locally, so this only happens once per gene - every following analysis reuses the cached copy. **Keep in mind that you must have enough disk memory for files.**.

* **Step 2: Strand direction**: Engine compares StartPos and EndPos data from NCBI. If StartPos is greater than EndPos, gene is on minus-strand. Else - plus-strand. This information is saved once in {Gene}_chrom_strand.txt.

* **Step 3: DNA Alignment**: Before anything is cut out, the patient's full sequence is aligned against the reference at the DNA level. This step accounts for insertions and deletions, not just simple letter swaps - it's what allows the next step to correctly locate coding regions even if the patient's sequence is a slightly different length than the reference.

* **Step 4: Applying the Annotation**: The tool picks the gene's canonical transcript (preferring the official MANE Select transcript, or the longest coding transcript if none is tagged) and uses the DNA alignment to map its coding regions (exons - CDSs) onto the patient's sequence.

* **Step 5: Protein Translation**: Both the reference and patient coding sequences are translated into protein chains. Translation stops at the first stop codon, the same way it naturally would inside a cell - so a mutation that shifts or removes the stop codon is captured accurately, instead of translating past it into meaningless letters.

* **Step 6: Protein Comparison**: The two protein chains are aligned side by side and compared position by position to find every difference.

* **Step 7: Nucleotide Translation**: Protein difference translates into codons and trims until bare nucleotides.

* **Step 8: Result Report**: All detected changes are listed in the main window. Results can be exported to a standard .vcf file.

---

## Reading the Results
Each row represents one detected difference between the patient and the reference protein. Mutation types include simple substitutions, insertions, deletions, and stop-codon changes (a mutation that creates or removes the protein's stop signal).

**If a row is highlighted:** The mutation matches a known entry in ClinVar, a public clinical variant database. Clicking it opens a detail window showing the disease association and a color-coded significance badge (red/orange for pathogenic, yellow for uncertain, green for benign).

**If a row is dim:** No matching entry was found in the public database. This does not mean the mutation is harmless - it simply means there's no public record to compare it against yet. It could also mean the app didn't find a match due to an issue on its side. Position is local, Ref/Alt is designated as proteins (not nucliotides).

---

## Exporting Results
Results can be saved as a `.vcf` file with a header and standard columns, so they can be opened in other genomics tools.

---

## Current Limitations (BETA)
* **One gene per file.** Whole-chromosome input isn't supported.
* **One transcript per gene.** If a gene has multiple isoforms, only the canonical one is analyzed.
* **Frameshift mutations** are not labeled as their own category yet - they typically show up as a cluster of substitutions after the shift point, or as a moved stop codon. Treat long runs of unexpected substitutions as a possible sign of a frameshift.
* **Large structural rearrangements** (e.g. inversions) are not supported - the alignment step assumes the patient and reference sequences are otherwise in the same order.