/*
 * dna_utils.h
 * -----------
 * Header file for DNA utility functions used in the Genetic Analyzer Toolkit.
 * Declares all function prototypes operating on DNA/RNA sequences.
 *
 * Author : Genetic Analyzer Toolkit Project
 * License: MIT
 */

#ifndef DNA_UTILS_H
#define DNA_UTILS_H

#include <stddef.h>  /* for size_t */

/* ── Constants ─────────────────────────────────────────────────────────── */
#define MAX_SEQ_LEN   4096   /* Maximum supported sequence length          */
#define MAX_PATTERN   256    /* Maximum pattern length for searching        */
#define MAX_AMINO     512    /* Maximum amino-acid output string length     */

/* ── Validation ────────────────────────────────────────────────────────── */

/**
 * validate_dna - Checks every character in seq is one of A, T, G, C
 * @seq: null-terminated string to validate
 * Returns 1 if valid, 0 otherwise.
 */
int validate_dna(const char *seq);

/* ── Core Transformations ───────────────────────────────────────────────── */

/**
 * dna_to_rna - Converts a DNA strand to its mRNA transcript.
 * Rule: A→U, T→A, G→C, C→G  (reads template strand, produces mRNA)
 * @dna: input DNA sequence (uppercase)
 * @rna: output buffer (must be at least strlen(dna)+1 bytes)
 */
void dna_to_rna(const char *dna, char *rna);

/**
 * complement_dna - Generates the complementary DNA strand.
 * Rule: A↔T, G↔C
 * @dna:   input sequence
 * @comp:  output buffer
 */
void complement_dna(const char *dna, char *comp);

/**
 * reverse_sequence - Reverses a DNA sequence in place copy.
 * @dna: input sequence
 * @rev: output buffer
 */
void reverse_sequence(const char *dna, char *rev);

/* ── Analytical Functions ───────────────────────────────────────────────── */

/**
 * gc_content - Calculates GC content as a percentage.
 * @dna: input sequence
 * Returns GC percentage (0.0 – 100.0).
 */
double gc_content(const char *dna);

/**
 * nucleotide_frequency - Counts occurrences of A, T, G, C.
 * @dna: input sequence
 * @a, @t, @g, @c: output pointers for counts
 */
void nucleotide_frequency(const char *dna, int *a, int *t, int *g, int *c);

/**
 * detect_mutations - Compares two equal-length DNA sequences.
 * Prints position, reference base, and mutant base for each mismatch.
 * @ref: reference sequence
 * @mut: mutant sequence
 * Returns the total number of mutations found.
 */
int detect_mutations(const char *ref, const char *mut);

/* ── Translation ────────────────────────────────────────────────────────── */

/**
 * codon_to_amino - Translates a single 3-letter codon to an amino-acid name.
 * @codon: 3-character string (uppercase)
 * Returns a string literal (amino acid name or "STOP" / "Unknown").
 */
const char *codon_to_amino(const char codon[3]);

/**
 * translate_sequence - Translates a full DNA coding sequence.
 * Begins translation at the first ATG (start codon).
 * Stops at TAA, TAG, or TGA (stop codons).
 * @dna:    input sequence (length must be multiple of 3)
 * @output: buffer to hold comma-separated amino acid names
 */
void translate_sequence(const char *dna, char *output);

/* ── Pattern Search ─────────────────────────────────────────────────────── */

/**
 * pattern_search - Finds all occurrences of a pattern within a DNA sequence.
 * Uses a simple sliding-window (naive) search.
 * @dna:     haystack sequence
 * @pattern: needle pattern
 * Returns the total number of matches found (0 if none).
 */
int pattern_search(const char *dna, const char *pattern);

#endif /* DNA_UTILS_H */
