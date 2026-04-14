# 🧬 Genetic Analyzer Toolkit (GAT)

> A professional, modular, menu-driven **bioinformatics toolkit** written in pure C.  
> Designed as a resume-worthy GitHub portfolio project demonstrating both C programming mastery and domain knowledge in molecular biology.

---

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

## 🔧 Building & Running

### Requirements
- GCC (≥ 9.0 recommended) or any C11-compliant compiler  
- Windows: [MinGW-w64](https://www.mingw-w64.org/) or [MSYS2](https://www.msys2.org/)  
- Linux/macOS: GCC available via package manager

### Compile (manual)
```bash
gcc -Wall -Wextra -std=c11 -o gat main.c dna_utils.c file_handler.c
```

### Compile (with Makefile)
```bash
make          # build
make run      # build + launch immediately
make clean    # remove compiled output
```

### Run
```bash
./gat          # Linux / macOS
gat.exe        # Windows
```

---

## 📸 Sample Session

```
  ╔══════════════════════════════════════════════════╗
  ║                                                  ║
  ║    G E N E T I C   A N A L Y Z E R              ║
  ║               T O O L K I T   v1.0              ║
  ║                                                  ║
  ║  DNA Analysis · Translation · Mutation Scan      ║
  ║                                                  ║
  ╚══════════════════════════════════════════════════╝

  ╔══════ MAIN MENU ════════════════════════╗
  ║  1.  DNA → RNA Transcription            ║
  ║  2.  Complementary DNA Strand           ║
  ...
```

---

## 🔬 Sample Inputs & Outputs

### Option 1 — DNA → RNA

```
Enter DNA sequence: ATGCTAGCTA
DNA  (5'→3'): ATGCTAGCTA
mRNA (5'→3'): UACGAUCGAU
```

### Option 3 — GC Content

```
Enter DNA sequence: ATGCGCATGC
Sequence length : 10 bp
GC Content      : 60.00%
AT Content      : 40.00%

GC  [███████████████░░░░░░░░░] 60.00%
AT  [██████████░░░░░░░░░░░░░░░] 40.00%
```

### Option 4 — Mutation Detection

```
Enter REFERENCE DNA sequence: ATGCGATCGA
Enter MUTANT   DNA sequence:  ATGCTATCCA

Position   Reference  Mutant
--------   ---------  ------
6          G          T
9          G          C

Total mutations found: 2
```

### Option 5 — Translation

```
Enter DNA coding sequence: ATGTTTTGGAAATAA
Protein chain:
Methionine (Start) → Phenylalanine → Tryptophan → Lysine → STOP
```

### Option 6 — Pattern Search

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

### Option 9 — Load from File

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

## 🛠 Implementation Highlights

- **Modular architecture** — Logic is cleanly separated across `main.c`, `dna_utils.c`, and `file_handler.c`
- **Input validation** — Every sequence is validated before processing; up to 3 retry attempts on bad input
- **Codon table** — Full 64-codon standard genetic code lookup implemented as a struct array
- **File support** — Accepts plain text and FASTA-style files; comment lines (`#`) and headers (`>`) are safely skipped
- **Callback pattern** — Batch file processing uses a function pointer callback for extensibility
- **ANSI colours** — Colour-coded terminal output with graceful Windows 10+ VT100 support
- **No external dependencies** — Compiles with any standard C11 compiler, no libraries required

---

## 📄 License

MIT License — free to use, modify, and distribute with attribution.

---

## 👤 Author

Built as part of a **C Programming & Bioinformatics Portfolio Project**.  
Feel free to fork, extend, and submit pull requests!
