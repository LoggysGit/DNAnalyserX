# GeneAnalyzerX - Quick Guide

GeneAnalyzerX is a simple tool that helps you find mutations in a patient's DNA. It takes the patient's genetic data, automatically compares it to a healthy standard template from the internet (NCBI), and shows you exactly what went wrong.

---

## What You Need to Load
* **Patient File (`.fasta.gz`):** The actual DNA sequence file of the patient you want to test.
* **Gene Name:** The official name of the gene you are investigating (for example: **TP53** or **BRCA1**).

---

## How It Works

* **Step 1: Automatic Download** - You type in the gene name, and the app automatically fetches the correct healthy map and reference sequence from official NCBI database servers.

* **Step 2: Exon Extraction** - The engine cuts out only the important coding parts (exons) of the gene from both the healthy template and your patient's file.

* **Step 3: Protein Translation** - The program converts these DNA strings into protein chains so we can see the real-world impact of any changes.

* **Step 4: Comparison & Alignment** - The core compares the healthy protein and the patient's protein side-by-side, letter by letter.

* **Step 5: Result Report** - The app highlights all mutations (swapped letters or shifted frames) in the main window grid and lets you export them into a standard clinical `.vcf` file.