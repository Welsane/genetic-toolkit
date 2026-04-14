/*
 * main.c
 * ------
 * Entry point and menu driver for the Genetic Analyzer Toolkit (GAT).
 *
 * Menu Options
 * ════════════
 *  1. DNA → RNA Transcription
 *  2. Complementary DNA Strand
 *  3. GC Content Calculation
 *  4. Mutation Detection
 *  5. Codon / Amino Acid Translation
 *  6. Pattern Search
 *  7. Reverse Sequence
 *  8. Nucleotide Frequency
 *  9. Load Sequence from File
 * 10. Batch Process File (all sequences)
 *  0. Exit
 *
 * Compile:
 *   gcc -Wall -Wextra -std=c11 -o gat main.c dna_utils.c file_handler.c
 *
 * Run:
 *   ./gat          (Linux / macOS)
 *   gat.exe        (Windows)
 *
 * Author : Genetic Analyzer Toolkit Project
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "dna_utils.h"
#include "file_handler.h"

/* ── Terminal colour codes (ANSI — fallback gracefully on Windows CMD) ─── */
#define CLR_RESET  "\033[0m"
#define CLR_BOLD   "\033[1m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_RED    "\033[1;31m"
#define CLR_BLUE   "\033[1;34m"
#define CLR_PURPLE "\033[1;35m"

/* ── Utility macros ────────────────────────────────────────────────────── */
#define DIVIDER \
    printf("  " CLR_CYAN \
           "═══════════════════════════════════════════════════\n" \
           CLR_RESET)

#define THIN_DIV \
    printf("  " CLR_BLUE \
           "───────────────────────────────────────────────────\n" \
           CLR_RESET)

/* ══════════════════════════════════════════════════════════════════════════
 *  INPUT HELPERS
 * ══════════════════════════════════════════════════════════════════════════ */

/**
 * read_line - Reads a line from stdin, strips the trailing newline.
 * Returns 1 on success, 0 on EOF / error.
 */
static int read_line(char *buf, size_t size)
{
    if (!fgets(buf, (int)size, stdin)) return 0;
    buf[strcspn(buf, "\r\n")] = '\0';
    return 1;
}

/**
 * get_dna_sequence - Prompts the user, reads a sequence, and validates it.
 * Retries up to 3 times on invalid input.
 * Returns 1 if a valid sequence was captured, 0 otherwise.
 */
static int get_dna_sequence(const char *prompt, char *buf)
{
    for (int attempt = 0; attempt < 3; attempt++) {
        printf("  %s", prompt);
        read_line(buf, MAX_SEQ_LEN);

        /* Convert to uppercase before validation */
        for (char *p = buf; *p; p++) *p = (char)toupper((unsigned char)*p);

        if (validate_dna(buf)) return 1;
        printf("  " CLR_RED "[!] Please enter only A, T, G, C characters.\n"
               CLR_RESET);
    }
    printf("  " CLR_RED "[!] Too many invalid attempts.\n" CLR_RESET);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  ASCII BANNER
 * ══════════════════════════════════════════════════════════════════════════ */

static void print_banner(void)
{
    printf(CLR_CYAN);
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════╗\n");
    printf("  ║                                                  ║\n");
    printf("  ║    " CLR_BOLD "G E N E T I C   A N A L Y Z E R" CLR_CYAN
           "                ║\n");
    printf("  ║               " CLR_BOLD "T O O L K I T   v1.0" CLR_CYAN
           "              ║\n");
    printf("  ║                                                  ║\n");
    printf("  ║  " CLR_YELLOW "DNA Analysis · Translation · Mutation Scan" CLR_CYAN
           "   ║\n");
    printf("  ║                                                  ║\n");
    printf("  ╚══════════════════════════════════════════════════╝\n");
    printf(CLR_RESET "\n");
}

/* ══════════════════════════════════════════════════════════════════════════
 *  MAIN MENU
 * ══════════════════════════════════════════════════════════════════════════ */

static void print_menu(void)
{
    printf(CLR_BOLD CLR_CYAN
           "\n  ╔══════ MAIN MENU ═════════════════════════╗\n" CLR_RESET);
    printf("  ║  " CLR_GREEN "1." CLR_RESET "  DNA → RNA Transcription             ║\n");
    printf("  ║  " CLR_GREEN "2." CLR_RESET "  Complementary DNA Strand             ║\n");
    printf("  ║  " CLR_GREEN "3." CLR_RESET "  GC Content Calculation               ║\n");
    printf("  ║  " CLR_GREEN "4." CLR_RESET "  Mutation Detection                   ║\n");
    printf("  ║  " CLR_GREEN "5." CLR_RESET "  Codon / Amino Acid Translation       ║\n");
    printf("  ║  " CLR_GREEN "6." CLR_RESET "  Pattern Search                       ║\n");
    printf("  ║  " CLR_GREEN "7." CLR_RESET "  Reverse Sequence                     ║\n");
    printf("  ║  " CLR_GREEN "8." CLR_RESET "  Nucleotide Frequency                 ║\n");
    printf("  ║  " CLR_GREEN "9." CLR_RESET "  Load Single Sequence from File       ║\n");
    printf("  ║  " CLR_GREEN "10." CLR_RESET " Batch Process File (all sequences)   ║\n");
    printf("  ║  " CLR_RED "0." CLR_RESET "  Exit                                ║\n");
    printf(CLR_CYAN "  ╚═════════════════════════════════════════╝\n" CLR_RESET);
    printf("  " CLR_BOLD "» " CLR_RESET "Enter choice: ");
}

/* ══════════════════════════════════════════════════════════════════════════
 *  MENU OPTION HANDLERS
 * ══════════════════════════════════════════════════════════════════════════ */

/* Option 1 ── DNA → RNA ─────────────────────────────────────────────────── */
static void menu_dna_to_rna(void)
{
    char dna[MAX_SEQ_LEN];
    char rna[MAX_SEQ_LEN];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "DNA → RNA TRANSCRIPTION\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence: ", dna)) return;

    dna_to_rna(dna, rna);

    printf("\n  " CLR_GREEN "DNA  (5'→3'):" CLR_RESET " %s\n", dna);
    printf("  " CLR_GREEN "mRNA (5'→3'):" CLR_RESET " %s\n", rna);
    printf("\n  " CLR_YELLOW "Rule applied: A→U  T→A  G→C  C→G\n" CLR_RESET);
    DIVIDER;
}

/* Option 2 ── Complementary DNA ─────────────────────────────────────────── */
static void menu_complement(void)
{
    char dna[MAX_SEQ_LEN];
    char comp[MAX_SEQ_LEN];
    char rev_comp[MAX_SEQ_LEN];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "COMPLEMENTARY DNA STRAND\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence (5'→3'): ", dna)) return;

    complement_dna(dna, comp);
    reverse_sequence(comp, rev_comp);

    printf("\n  " CLR_GREEN "Original        (5'→3'):" CLR_RESET " %s\n", dna);
    printf("  " CLR_GREEN "Complement      (3'→5'):" CLR_RESET " %s\n", comp);
    printf("  " CLR_GREEN "Rev. Complement (5'→3'):" CLR_RESET " %s\n", rev_comp);
    printf("\n  " CLR_YELLOW "Rule applied: A↔T  G↔C\n" CLR_RESET);
    DIVIDER;
}

/* Option 3 ── GC Content ────────────────────────────────────────────────── */
static void menu_gc_content(void)
{
    char dna[MAX_SEQ_LEN];
    int  a, t, g, c;

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "GC CONTENT ANALYSIS\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence: ", dna)) return;

    double gc = gc_content(dna);
    nucleotide_frequency(dna, &a, &t, &g, &c);
    int total = a + t + g + c;

    printf("\n  Sequence length : %d bp\n", total);
    printf("  GC Content      : " CLR_BOLD CLR_GREEN "%.2f%%\n" CLR_RESET, gc);
    printf("  AT Content      : " CLR_BOLD CLR_YELLOW "%.2f%%\n" CLR_RESET,
           100.0 - gc);

    /* Simple ASCII bar chart */
    printf("\n  GC  [");
    int filled = (int)(gc / 4);   /* 25 chars = 100% */
    for (int i = 0; i < 25; i++) printf(i < filled ? "█" : "░");
    printf("] %.2f%%\n", gc);

    printf("  AT  [");
    int at_fill = 25 - filled;
    for (int i = 0; i < 25; i++) printf(i < at_fill ? "█" : "░");
    printf("] %.2f%%\n", 100.0 - gc);

    DIVIDER;
}

/* Option 4 ── Mutation Detection ────────────────────────────────────────── */
static void menu_mutations(void)
{
    char ref[MAX_SEQ_LEN];
    char mut[MAX_SEQ_LEN];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "MUTATION DETECTION\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter REFERENCE DNA sequence: ", ref)) return;
    if (!get_dna_sequence("Enter MUTANT   DNA sequence: ", mut)) return;

    int count = detect_mutations(ref, mut);

    if (count > 0) {
        printf("\n  " CLR_RED "Total mutations found: %d\n" CLR_RESET, count);
    } else if (count == 0) {
        printf("\n  " CLR_GREEN "Sequences are identical.\n" CLR_RESET);
    }
    DIVIDER;
}

/* Option 5 ── Translation ───────────────────────────────────────────────── */
static void menu_translate(void)
{
    char dna[MAX_SEQ_LEN];
    char amino[MAX_AMINO];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "CODON / AMINO ACID TRANSLATION\n" CLR_RESET);
    THIN_DIV;
    printf("  " CLR_YELLOW "Note: Length must be a multiple of 3.\n" CLR_RESET);
    printf("  " CLR_YELLOW "Translation starts at the first ATG (start codon).\n\n"
           CLR_RESET);

    if (!get_dna_sequence("Enter DNA coding sequence: ", dna)) return;

    translate_sequence(dna, amino);

    printf("\n  " CLR_GREEN "Protein chain:\n" CLR_RESET);
    printf("  %s\n", amino);
    DIVIDER;
}

/* Option 6 ── Pattern Search ────────────────────────────────────────────── */
static void menu_pattern(void)
{
    char dna[MAX_SEQ_LEN];
    char pattern[MAX_PATTERN];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "PATTERN SEARCH\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence (haystack): ", dna)) return;

    printf("  Enter pattern to search (needle): ");
    read_line(pattern, sizeof(pattern));

    /* Uppercase the pattern */
    for (char *p = pattern; *p; p++) *p = (char)toupper((unsigned char)*p);

    if (!validate_dna(pattern)) {
        printf("  " CLR_RED "[!] Invalid pattern — must contain only A, T, G, C.\n"
               CLR_RESET);
        DIVIDER;
        return;
    }

    int hits = pattern_search(dna, pattern);

    printf("\n  " CLR_GREEN "Total occurrences of \"%s\": %d\n" CLR_RESET,
           pattern, hits);
    DIVIDER;
}

/* Option 7 ── Reverse Sequence ──────────────────────────────────────────── */
static void menu_reverse(void)
{
    char dna[MAX_SEQ_LEN];
    char rev[MAX_SEQ_LEN];
    char rev_comp[MAX_SEQ_LEN];
    char comp[MAX_SEQ_LEN];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "REVERSE SEQUENCE\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence: ", dna)) return;

    reverse_sequence(dna, rev);
    complement_dna(dna, comp);
    reverse_sequence(comp, rev_comp);

    printf("\n  " CLR_GREEN "Original             :" CLR_RESET " %s\n", dna);
    printf("  " CLR_GREEN "Reversed             :" CLR_RESET " %s\n", rev);
    printf("  " CLR_GREEN "Reverse Complement   :" CLR_RESET " %s\n", rev_comp);
    DIVIDER;
}

/* Option 8 ── Nucleotide Frequency ──────────────────────────────────────── */
static void menu_frequency(void)
{
    char dna[MAX_SEQ_LEN];
    int  a, t, g, c;

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "NUCLEOTIDE FREQUENCY\n" CLR_RESET);
    THIN_DIV;

    if (!get_dna_sequence("Enter DNA sequence: ", dna)) return;

    nucleotide_frequency(dna, &a, &t, &g, &c);
    int total = a + t + g + c;

    printf("\n  %-12s %-8s %s\n",  "Nucleotide", "Count", "Percentage");
    printf("  %-12s %-8s %s\n",    "──────────", "─────", "──────────");
    printf("  %-12s %-8d %.2f%%\n", "Adenine (A)", a,
           total ? a * 100.0 / total : 0.0);
    printf("  %-12s %-8d %.2f%%\n", "Thymine (T)", t,
           total ? t * 100.0 / total : 0.0);
    printf("  %-12s %-8d %.2f%%\n", "Guanine (G)", g,
           total ? g * 100.0 / total : 0.0);
    printf("  %-12s %-8d %.2f%%\n", "Cytosine (C)", c,
           total ? c * 100.0 / total : 0.0);
    printf("  %-12s %-8d\n", "Total", total);
    DIVIDER;
}

/* Option 9 ── Load single sequence from file ────────────────────────────── */
static void menu_load_file(void)
{
    char filepath[512];
    char dna[MAX_SEQ_LEN];
    char rna[MAX_SEQ_LEN];
    char comp[MAX_SEQ_LEN];
    char amino[MAX_AMINO];
    int  a, t, g, c;

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "LOAD SEQUENCE FROM FILE\n" CLR_RESET);
    THIN_DIV;
    printf("  Enter file path (e.g. sample_dna.txt): ");
    read_line(filepath, sizeof(filepath));

    if (!load_sequence_from_file(filepath, dna)) {
        DIVIDER;
        return;
    }

    if (!validate_dna(dna)) {
        printf("  " CLR_RED "[!] File contains an invalid DNA sequence.\n"
               CLR_RESET);
        DIVIDER;
        return;
    }

    printf("\n  " CLR_GREEN "Loaded sequence:" CLR_RESET " %s\n", dna);
    printf("  Length: %zu bp\n\n", strlen(dna));

    /* Run quick-analysis suite */
    dna_to_rna(dna, rna);
    complement_dna(dna, comp);
    nucleotide_frequency(dna, &a, &t, &g, &c);
    double gc = gc_content(dna);
    translate_sequence(dna, amino);

    THIN_DIV;
    printf("  " CLR_YELLOW "── Quick Analysis Report ──\n" CLR_RESET);
    printf("  GC Content     : %.2f%%\n", gc);
    printf("  A:%d  T:%d  G:%d  C:%d\n", a, t, g, c);
    printf("  mRNA           : %s\n", rna);
    printf("  Complement     : %s\n", comp);
    printf("  Translation    : %s\n", amino);
    DIVIDER;
}

/* Option 10 ── Batch process all sequences in file ─────────────────────── */
static void batch_callback(const char *seq, int index, void *user_data)
{
    (void)user_data;   /* unused */

    char rna[MAX_SEQ_LEN];
    int  a, t, g, c;

    printf("\n  " CLR_BOLD CLR_CYAN "── Sequence #%d ──────────────────────────\n"
           CLR_RESET, index + 1);
    printf("  DNA  : %s\n", seq);
    printf("  Len  : %zu bp\n", strlen(seq));

    dna_to_rna(seq, rna);
    printf("  mRNA : %s\n", rna);

    nucleotide_frequency(seq, &a, &t, &g, &c);
    double gc = gc_content(seq);
    printf("  GC%%  : %.2f%%\n", gc);
    printf("  Freq : A=%d  T=%d  G=%d  C=%d\n", a, t, g, c);
}

static void menu_batch_file(void)
{
    char filepath[512];

    DIVIDER;
    printf("  " CLR_BOLD CLR_PURPLE "BATCH PROCESS FILE\n" CLR_RESET);
    THIN_DIV;
    printf("  " CLR_YELLOW "Each non-comment line is treated as one DNA sequence.\n"
           CLR_RESET);
    printf("  Enter file path: ");
    read_line(filepath, sizeof(filepath));

    int count = load_all_sequences_from_file(filepath, batch_callback, NULL);

    if (count >= 0) {
        printf("\n  " CLR_GREEN "Processed %d sequence(s) from '%s'.\n"
               CLR_RESET, count, filepath);
    }
    DIVIDER;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  MAIN
 * ══════════════════════════════════════════════════════════════════════════ */

int main(void)
{
    /* Enable VT100 colour codes on Windows 10+ */
#ifdef _WIN32
    system("chcp 65001 > nul");   /* UTF-8 console */
#endif

    print_banner();

    char input[16];
    int  choice;
    int  running = 1;

    while (running) {
        print_menu();

        if (!read_line(input, sizeof(input))) break;

        /* Attempt to parse the choice as an integer */
        if (sscanf(input, "%d", &choice) != 1) {
            printf("  " CLR_RED "[!] Invalid input. Please enter a number.\n"
                   CLR_RESET);
            continue;
        }

        printf("\n");

        switch (choice) {
            case 1:  menu_dna_to_rna();   break;
            case 2:  menu_complement();   break;
            case 3:  menu_gc_content();   break;
            case 4:  menu_mutations();    break;
            case 5:  menu_translate();    break;
            case 6:  menu_pattern();      break;
            case 7:  menu_reverse();      break;
            case 8:  menu_frequency();    break;
            case 9:  menu_load_file();    break;
            case 10: menu_batch_file();   break;
            case 0:
                printf("  " CLR_GREEN
                       "Thank you for using the Genetic Analyzer Toolkit!\n"
                       CLR_RESET);
                running = 0;
                break;
            default:
                printf("  " CLR_RED
                       "[!] Unknown option. Please choose 0–10.\n"
                       CLR_RESET);
        }
    }

    return 0;
}
