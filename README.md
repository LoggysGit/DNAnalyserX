# GeneAnalyserX v1.0

GeneAnalyserX is a software tool for analyzing local segments of human DNA. It compares a patient's sequence for a single gene against a reference from NCBI, detects point mutations, cross-references them against the ClinVar database, and reports the genetic diseases associated with each match.

* * *

### Features

Key capabilities and design choices:

*   **Purpose** - accurate comparison of a patient's gene sequence against a reference using pairwise sequence alignment, covering both simple substitutions and insertions/deletions.
*   **Reference handling** - reference sequences and annotations are fetched from NCBI automatically and cached locally per gene, so repeated analyses of the same gene don't require a new download.
*   **Interface** - a lightweight, straightforward interface built with CustomTkinter.
*   **Disease database** - local database of known variants and their clinical significance, cross-referenced against ClinVar and refreshed on a weekly basis.
*   **Analysis export** - results can be exported as a strict `.vcf` file, ready to open in other genomics tools.

* * *

### Requirements

Make sure you have installed:

*   Git
*   Python (preferably via VS Code)
*   A package manager (pip)

### Setup Instructions

Follow these steps to run the project locally.

1.  **Clone the repository:**

        git clone https://github.com/LoggysGit/GeneAnalyserX
        cd GeneAnalyserX

2.  **Create a virtual environment:**

        python -m venv .venv

3.  **Activate the virtual environment:**

    Windows PowerShell:

        Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
        .venv\Scripts\Activate.ps1

    Git Bash:

        source .venv/Scripts/activate

4.  **Install dependencies:**

        pip install -r requirements.txt

5.  **Set your real email in `config.cfg`**

6.  **Run the application:**

        python app.py

### Usage

1. Launch the application.
2. Load the patient's gene sequence file in `.fasta.gz` format. The file must contain exactly one sequence - a single gene, not a full chromosome.
3. Enter the official gene name (for example, `TP53`).
4. Click **Analyse**. The app will download and cache the reference automatically, align the patient's sequence, translate it, and detect mutations.
5. Review the results in the main window and export them if needed.

### How the Analysis Works

*   **One gene, one file.** Each analysis run covers a single gene; the input file must contain exactly one sequence.
*   **A real email is required.** `config.cfg` must contain a valid email address - NCBI requires this to identify API requests, and a placeholder address can get requests blocked.
*   **A single canonical transcript is used per gene** - the official MANE Select transcript when tagged, otherwise the transcript with the longest coding sequence.
*   **Translation stops at the first stop codon**, matching how translation actually terminates in a cell, so stop-codon mutations (premature stops, stop-loss) are captured correctly instead of translating through them.
*   **Detected mutation types** currently include substitutions, insertions, deletions, premature stop codons (nonsense), and stop-loss variants. Frameshift mutations are not labeled as a distinct category - they typically surface as a run of substitutions after the shift point or a moved stop codon.
*   **Not supported:** whole-chromosome input, multiple transcripts per gene, and large structural rearrangements (e.g. inversions) - the aligner assumes the patient and reference sequences are otherwise in the same order.

_Note: if errors or crashes occur, the full session history is saved to `logs.log`. To adjust working parameters, edit the `config.cfg` file._

* * *

### Sources Used During Development

* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
* [Wikipedia](https://www.wikipedia.org/)
* [StackOverflow](https://stackoverflow.com/questions)


### Built With
 
* [Biopython](https://biopython.org/)
* [CustomTkinter](https://github.com/TomSchimansky/CustomTkinter)

### Authors

**Metsler Albert** - *Developer*, 2026