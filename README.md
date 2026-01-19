# Phylogeny-Analysis-Tool
Python application for automated sequence alignment and phylogenetic tree construction
# 🧬 Phylogenetic Analysis Tool

A Python application to automate the retrieval of biological sequences from NCBI, align them using MUSCLE, and generate phylogenetic trees with genetic divergence analysis.

## 🚀 Features

* **Automated Download:** Fetches protein sequences from NCBI using accession numbers.
* **Data Cleaning:** Automatically removes duplicate species and formats names.
* **Sequence Alignment:** Uses **MUSCLE** algorithm (included).
* **Analysis:** Calculates genetic divergence percentages matrix.
* **Visualization:** * Draws Phylogenetic Trees (Neighbor Joining).
    * Displays Clustal format alignment text.
* **User Interface:** Full graphical interface (no code required).

## 📥 How to Download (For Windows)

1. Go to the **[Releases Section](../../releases)** on the right.
2. Download the **ZIP file** from the latest version (v1.0).
3. Unzip the folder.
4. Run `main.exe`.
   * *Note: Ensure `muscle.exe` is in the same folder as the application.*

## 🛠️ How to Use

1. **Configuration:**
   * Enter your email (required by NCBI).
   * Select your CSV Data file.
   * Select the `muscle.exe` file path.
   * Click **LOAD**.
2. **Execution:**
   * Select a Reference Species for comparison.
   * Click **RUN ANALYSIS**.
3. **Results:**
   * Select species in the table (Ctrl + Click).
   * Click **DRAW TREE** to see the graph.

## 📋 Prerequisites (For running from source code)

If you want to run the python script `main.py` instead of the exe:
* Python 3.x
* Libraries: `biopython`, `pandas`, `matplotlib`, `tk`
