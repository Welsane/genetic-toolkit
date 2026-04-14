/*
 * file_handler.c
 * --------------
 * File I/O utilities for the Genetic Analyzer Toolkit.
 *
 * Supports:
 *  - Plain text files (one sequence per line)
 *  - FASTA-like files (header lines starting with '>' are skipped)
 *  - Comment lines starting with '#' are ignored
 *
 * Author : Genetic Analyzer Toolkit Project
 * License: MIT
 */

#include "file_handler.h"
#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ── Internal: strips trailing whitespace/newline from a string ────────── */
static void rtrim(char *s)
{
    if (!s) return;
    size_t len = strlen(s);
    while (len > 0 &&
           (s[len-1] == '\n' || s[len-1] == '\r' ||
            s[len-1] == ' '  || s[len-1] == '\t')) {
        s[--len] = '\0';
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 *  load_sequence_from_file
 *  Reads exactly one DNA sequence from a file (first valid line).
 * ══════════════════════════════════════════════════════════════════════════ */
int load_sequence_from_file(const char *filepath, char *buffer)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        printf("  [!] Error: Cannot open file '%s'.\n", filepath);
        return 0;
    }

    char line[MAX_SEQ_LEN];
    int  found = 0;

    while (fgets(line, sizeof(line), fp)) {
        rtrim(line);

        /* Skip comments and FASTA headers */
        if (line[0] == '#' || line[0] == '>' || line[0] == '\0') continue;

        /* Convert to uppercase in-place */
        for (char *p = line; *p; p++) *p = (char)toupper((unsigned char)*p);

        strncpy(buffer, line, MAX_SEQ_LEN - 1);
        buffer[MAX_SEQ_LEN - 1] = '\0';
        found = 1;
        break;
    }

    fclose(fp);

    if (!found) {
        printf("  [!] No valid DNA sequence found in '%s'.\n", filepath);
        return 0;
    }
    return 1;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  load_all_sequences_from_file
 *  Iterates every valid sequence line and fires a callback.
 * ══════════════════════════════════════════════════════════════════════════ */
int load_all_sequences_from_file(const char     *filepath,
                                 seq_callback_t  callback,
                                 void           *user_data)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        printf("  [!] Error: Cannot open file '%s'.\n", filepath);
        return -1;
    }

    char line[MAX_SEQ_LEN];
    int  index = 0;

    while (fgets(line, sizeof(line), fp)) {
        rtrim(line);

        /* Skip comments, FASTA headers, blank lines */
        if (line[0] == '#' || line[0] == '>' || line[0] == '\0') continue;

        /* Convert to uppercase */
        for (char *p = line; *p; p++) *p = (char)toupper((unsigned char)*p);

        /* Validate before firing callback */
        if (!validate_dna(line)) {
            printf("  [!] Skipping line %d — invalid DNA sequence: %s\n",
                   index + 1, line);
            continue;
        }

        callback(line, index, user_data);
        index++;
    }

    fclose(fp);
    return index;
}
