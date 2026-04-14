/*
 * dna_utils.c
 * -----------
 * Implementation of all DNA/RNA analysis functions for the
 * Genetic Analyzer Toolkit.
 *
 * Each function is self-contained and operates on null-terminated
 * C strings.  All sequences are expected in UPPERCASE (A, T, G, C).
 *
 * Compile together with main.c and file_handler.c:
 *   gcc -Wall -Wextra -o gat main.c dna_utils.c file_handler.c
 *
 * Author : Genetic Analyzer Toolkit Project
 * License: MIT
 */

#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>   /* toupper */

/* ══════════════════════════════════════════════════════════════════════════
 *  INTERNAL HELPERS
 * ══════════════════════════════════════════════════════════════════════════ */

/**
 * to_upper_copy - Copies src into dst converting every char to uppercase.
 * Caller must ensure dst has enough space (strlen(src)+1).
 */
static void to_upper_copy(const char *src, char *dst)
{
    while (*src) {
        *dst++ = (char)toupper((unsigned char)*src++);
    }
    *dst = '\0';
}

/* ══════════════════════════════════════════════════════════════════════════
 *  VALIDATION
 * ══════════════════════════════════════════════════════════════════════════ */

int validate_dna(const char *seq)
{
    if (!seq || *seq == '\0') return 0;          /* empty is invalid */

    for (const char *p = seq; *p; p++) {
        char c = (char)toupper((unsigned char)*p);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
            printf("  [!] Invalid character '%c' found in sequence.\n", *p);
            return 0;
        }
    }
    return 1;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  CORE TRANSFORMATIONS
 * ══════════════════════════════════════════════════════════════════════════ */

/* --- DNA → RNA (mRNA transcript) ------------------------------------------
 * The DNA template strand is read 3'→5' and mRNA is synthesised 5'→3'.
 * For simplicity we directly apply the transcription rules on the coding
 * strand (same direction):
 *   A → U   T → A   G → C   C → G
 * -------------------------------------------------------------------------- */
void dna_to_rna(const char *dna, char *rna)
{
    char upper[MAX_SEQ_LEN];
    to_upper_copy(dna, upper);

    for (size_t i = 0; upper[i]; i++) {
        switch (upper[i]) {
            case 'A': rna[i] = 'U'; break;
            case 'T': rna[i] = 'A'; break;
            case 'G': rna[i] = 'C'; break;
            case 'C': rna[i] = 'G'; break;
            default:  rna[i] = '?'; break;
        }
    }
    rna[strlen(upper)] = '\0';
}

/* --- Complementary DNA ----------------------------------------------------
 * Watson-Crick base pairing:  A↔T  G↔C
 * -------------------------------------------------------------------------*/
void complement_dna(const char *dna, char *comp)
{
    char upper[MAX_SEQ_LEN];
    to_upper_copy(dna, upper);
    size_t len = strlen(upper);

    for (size_t i = 0; i < len; i++) {
        switch (upper[i]) {
            case 'A': comp[i] = 'T'; break;
            case 'T': comp[i] = 'A'; break;
            case 'G': comp[i] = 'C'; break;
            case 'C': comp[i] = 'G'; break;
            default:  comp[i] = '?'; break;
        }
    }
    comp[len] = '\0';
}

/* --- Reverse sequence -------------------------------------------------------
 * Returns a reversed copy of the sequence.
 * Commonly used alongside complement_dna to get the reverse complement.
 * -------------------------------------------------------------------------- */
void reverse_sequence(const char *dna, char *rev)
{
    size_t len = strlen(dna);
    for (size_t i = 0; i < len; i++) {
        rev[i] = dna[len - 1 - i];
    }
    rev[len] = '\0';
}

/* ══════════════════════════════════════════════════════════════════════════
 *  ANALYTICAL FUNCTIONS
 * ══════════════════════════════════════════════════════════════════════════ */

/* --- GC Content ------------------------------------------------------------
 * GC% = (G + C) / total_bases * 100
 * -------------------------------------------------------------------------- */
double gc_content(const char *dna)
{
    if (!dna || *dna == '\0') return 0.0;

    int gc = 0;
    int total = 0;

    for (const char *p = dna; *p; p++) {
        char c = (char)toupper((unsigned char)*p);
        total++;
        if (c == 'G' || c == 'C') gc++;
    }
    return (total == 0) ? 0.0 : (gc * 100.0 / total);
}

/* --- Nucleotide Frequency --------------------------------------------------
 * Counts A, T, G, C occurrences in the sequence.
 * -------------------------------------------------------------------------- */
void nucleotide_frequency(const char *dna, int *a, int *t, int *g, int *c)
{
    *a = *t = *g = *c = 0;

    for (const char *p = dna; *p; p++) {
        switch ((char)toupper((unsigned char)*p)) {
            case 'A': (*a)++; break;
            case 'T': (*t)++; break;
            case 'G': (*g)++; break;
            case 'C': (*c)++; break;
        }
    }
}

/* --- Mutation Detection ----------------------------------------------------
 * Compares ref and mut position by position.
 * Both sequences must have equal length.
 * -------------------------------------------------------------------------- */
int detect_mutations(const char *ref, const char *mut)
{
    size_t ref_len = strlen(ref);
    size_t mut_len = strlen(mut);

    if (ref_len != mut_len) {
        printf("  [!] Sequences must be the same length for mutation "
               "detection.\n"
               "      Reference: %zu bp   Mutant: %zu bp\n",
               ref_len, mut_len);
        return -1;
    }

    int count = 0;

    printf("\n  %-10s %-10s %-10s\n", "Position", "Reference", "Mutant");
    printf("  %-10s %-10s %-10s\n", "--------", "---------", "------");

    for (size_t i = 0; i < ref_len; i++) {
        char r = (char)toupper((unsigned char)ref[i]);
        char m = (char)toupper((unsigned char)mut[i]);
        if (r != m) {
            printf("  %-10zu %-10c %-10c\n", i + 1, r, m);
            count++;
        }
    }

    if (count == 0) {
        printf("  Sequences are identical — no mutations detected.\n");
    }

    return count;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  CODON TRANSLATION
 * ══════════════════════════════════════════════════════════════════════════ */

/* Standard genetic code lookup table (64 codons).
 * Only DNA codons are listed (T instead of U in mRNA).
 * Each entry: codon[3], amino_acid_name.
 */
static const struct {
    char codon[4];        /* 3-letter codon + null terminator */
    const char *amino;    /* Full amino acid name              */
} CODON_TABLE[] = {
    /* Phenylalanine */
    {"TTT", "Phenylalanine"}, {"TTC", "Phenylalanine"},
    /* Leucine */
    {"TTA", "Leucine"},       {"TTG", "Leucine"},
    {"CTT", "Leucine"},       {"CTC", "Leucine"},
    {"CTA", "Leucine"},       {"CTG", "Leucine"},
    /* Isoleucine */
    {"ATT", "Isoleucine"},    {"ATC", "Isoleucine"}, {"ATA", "Isoleucine"},
    /* Methionine / Start */
    {"ATG", "Methionine (Start)"},
    /* Valine */
    {"GTT", "Valine"},        {"GTC", "Valine"},
    {"GTA", "Valine"},        {"GTG", "Valine"},
    /* Serine */
    {"TCT", "Serine"},        {"TCC", "Serine"},
    {"TCA", "Serine"},        {"TCG", "Serine"},
    {"AGT", "Serine"},        {"AGC", "Serine"},
    /* Proline */
    {"CCT", "Proline"},       {"CCC", "Proline"},
    {"CCA", "Proline"},       {"CCG", "Proline"},
    /* Threonine */
    {"ACT", "Threonine"},     {"ACC", "Threonine"},
    {"ACA", "Threonine"},     {"ACG", "Threonine"},
    /* Alanine */
    {"GCT", "Alanine"},       {"GCC", "Alanine"},
    {"GCA", "Alanine"},       {"GCG", "Alanine"},
    /* Tyrosine */
    {"TAT", "Tyrosine"},      {"TAC", "Tyrosine"},
    /* Stop codons */
    {"TAA", "STOP"},          {"TAG", "STOP"},       {"TGA", "STOP"},
    /* Histidine */
    {"CAT", "Histidine"},     {"CAC", "Histidine"},
    /* Glutamine */
    {"CAA", "Glutamine"},     {"CAG", "Glutamine"},
    /* Asparagine */
    {"AAT", "Asparagine"},    {"AAC", "Asparagine"},
    /* Lysine */
    {"AAA", "Lysine"},        {"AAG", "Lysine"},
    /* Aspartate */
    {"GAT", "Aspartate"},     {"GAC", "Aspartate"},
    /* Glutamate */
    {"GAA", "Glutamate"},     {"GAG", "Glutamate"},
    /* Cysteine */
    {"TGT", "Cysteine"},      {"TGC", "Cysteine"},
    /* Tryptophan */
    {"TGG", "Tryptophan"},
    /* Arginine */
    {"CGT", "Arginine"},      {"CGC", "Arginine"},
    {"CGA", "Arginine"},      {"CGG", "Arginine"},
    {"AGA", "Arginine"},      {"AGG", "Arginine"},
    /* Glycine */
    {"GGT", "Glycine"},       {"GGC", "Glycine"},
    {"GGA", "Glycine"},       {"GGG", "Glycine"},
    /* Sentinel */
    {"", NULL}
};

const char *codon_to_amino(const char codon[3])
{
    /* Build a null-terminated string from the 3-char codon */
    char key[4];
    key[0] = (char)toupper((unsigned char)codon[0]);
    key[1] = (char)toupper((unsigned char)codon[1]);
    key[2] = (char)toupper((unsigned char)codon[2]);
    key[3] = '\0';

    for (int i = 0; CODON_TABLE[i].amino != NULL; i++) {
        if (strcmp(CODON_TABLE[i].codon, key) == 0) {
            return CODON_TABLE[i].amino;
        }
    }
    return "Unknown";
}

void translate_sequence(const char *dna, char *output)
{
    size_t len = strlen(dna);

    if (len % 3 != 0) {
        snprintf(output, MAX_AMINO,
                 "[Error] Sequence length (%zu) is not a multiple of 3. "
                 "Cannot translate incomplete codons.", len);
        return;
    }

    /* Scan for the start codon ATG */
    size_t start = len;   /* initialise to 'not found' */
    for (size_t i = 0; i + 2 < len; i += 3) {
        char c0 = (char)toupper((unsigned char)dna[i]);
        char c1 = (char)toupper((unsigned char)dna[i+1]);
        char c2 = (char)toupper((unsigned char)dna[i+2]);
        if (c0 == 'A' && c1 == 'T' && c2 == 'G') {
            start = i;
            break;
        }
    }

    if (start == len) {
        snprintf(output, MAX_AMINO,
                 "[Error] No start codon (ATG) found in sequence.");
        return;
    }

    /* Translate from start codon until stop or end */
    output[0] = '\0';
    int first = 1;

    for (size_t i = start; i + 2 < len; i += 3) {
        const char *amino = codon_to_amino(&dna[i]);

        if (!first) strncat(output, " → ", MAX_AMINO - strlen(output) - 1);
        strncat(output, amino, MAX_AMINO - strlen(output) - 1);
        first = 0;

        if (strcmp(amino, "STOP") == 0) break;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 *  PATTERN SEARCH  (naive sliding-window O(n·m))
 * ══════════════════════════════════════════════════════════════════════════ */

int pattern_search(const char *dna, const char *pattern)
{
    size_t dna_len = strlen(dna);
    size_t pat_len = strlen(pattern);
    int    count   = 0;

    if (pat_len == 0 || pat_len > dna_len) {
        printf("  No occurrences found.\n");
        return 0;
    }

    printf("\n  %-12s %s\n", "Position", "Match");
    printf("  %-12s %s\n",   "--------", "-----");

    for (size_t i = 0; i <= dna_len - pat_len; i++) {
        /* Case-insensitive comparison */
        int match = 1;
        for (size_t j = 0; j < pat_len; j++) {
            if (toupper((unsigned char)dna[i+j]) !=
                toupper((unsigned char)pattern[j])) {
                match = 0;
                break;
            }
        }
        if (match) {
            count++;
            printf("  %-12zu %.*s\n", i + 1, (int)pat_len, dna + i);
        }
    }

    if (count == 0) {
        printf("  Pattern not found in sequence.\n");
    }

    return count;
}
