# 🧬 Genetic Analyzer Toolkit (GAT)

> A professional, modular, menu-driven **bioinformatics toolkit** written in pure C.  

## 📋 Features

| # | Feature | Description |
|---|---------|-------------|
| 1 | **DNA → RNA Transcription** | Converts a DNA sequence to its mRNA transcript (A→U, T→A, G→C, C→G) |
| 2 | **Complementary DNA** | Generates the complementary strand (A↔T, G↔C) and its reverse complement |
| 3 | **GC Content Analysis** | Calculates GC% to 2 decimal places with an ASCII progress bar |
| 4 | **Mutation Detection** | Compares two equal-length sequences and reports every SNP with position and bases |
| 5 | **Codon Translation** | Translates a coding sequence from ATG to a stop codon using the standard genetic code |
| 6 | **Pattern Search** | Finds all occurrences of a motif within a DNA sequence (returns positions) |
| 7 | **Reverse Sequence** | Reverses a sequence and computes the reverse complement |
| 8 | **Nucleotide Frequency** | Counts A, T, G, C with percentages |
| 9 | **Load from File** | Reads a single DNA sequence from a text/FASTA file and runs a quick-analysis report |
| 10 | **Batch File Processing** | Processes every sequence in a file with full analysis on each |

---

## 🗂 Project Structure

```
genetic_analyzer_toolkit/
│
├── main.c           ← Menu driver and all menu-option handlers
├── dna_utils.c      ← Core DNA functions (transformations, analysis, translation, search)
├── dna_utils.h      ← Public API declarations for dna_utils.c
├── file_handler.c   ← File I/O: single and batch sequence loading
├── file_handler.h   ← Public API declarations for file_handler.c
│
├── sample_dna.txt   ← Example input file with annotated sequences
├── Makefile         ← Cross-platform build automation
└── README.md        ← This file
```

---

## — DNA → RNA

```
Enter DNA sequence: ATGCTAGCTA
DNA  (5'→3'): ATGCTAGCTA
mRNA (5'→3'): UACGAUCGAU
```

## — GC Content

```
Enter DNA sequence: ATGCGCATGC
Sequence length : 10 bp
GC Content      : 60.00%
AT Content      : 40.00%

GC  [███████████████░░░░░░░░░] 60.00%
AT  [██████████░░░░░░░░░░░░░░░] 40.00%
```

## — Mutation Detection

```
Enter REFERENCE DNA sequence: ATGCGATCGA
Enter MUTANT   DNA sequence:  ATGCTATCCA

Position   Reference  Mutant
--------   ---------  ------
6          G          T
9          G          C

Total mutations found: 2
```

## — Translation

```
Enter DNA coding sequence: ATGTTTTGGAAATAA
Protein chain:
Methionine (Start) → Phenylalanine → Tryptophan → Lysine → STOP
```

## — Pattern Search

```
Enter DNA sequence: ATGATCATGATCATG
Enter pattern:      ATG

Position     Match
--------     -----
1            ATG
6            ATG
11           ATG

Total occurrences of "ATG": 3
```

## — Load from File

```
Enter file path: sample_dna.txt
Loaded sequence: ATGTAA
Length: 6 bp

── Quick Analysis Report ──
GC Content  : 16.67%
A:2  T:2  G:1  C:0
mRNA        : UACAUU
Complement  : TACATT
Translation : Methionine (Start) → STOP
```

---

## 🧬 Biological Background

| Concept | Explanation |
|---------|-------------|
| **Transcription** | DNA is read as a template to produce mRNA. A→U, T→A, G→C, C→G |
| **Complementarity** | The two DNA strands are anti-parallel and complementary (A↔T, G↔C) |
| **GC Content** | Higher GC% = higher melting temperature; used in primer design and taxonomy |
| **SNP (Mutation)** | A single nucleotide polymorphism is a position where sequences differ |
| **Codon** | A triplet of nucleotides that codes for one amino acid |
| **ORF** | Open Reading Frame — begins at ATG and ends at a stop codon (TAA/TAG/TGA) |
| **Motif Search** | Finding conserved sequence patterns is essential in regulatory element analysis |

---
