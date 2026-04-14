/*
 * file_handler.h
 * --------------
 * Header for file I/O utilities used in the Genetic Analyzer Toolkit.
 * Provides functions to load and validate DNA sequences from text files.
 *
 * Author : Genetic Analyzer Toolkit Project
 * License: MIT
 */

#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H

#include "dna_utils.h"   /* shares MAX_SEQ_LEN */

/**
 * load_sequence_from_file
 * -----------------------
 * Opens a plain-text file and reads the first non-empty, non-comment line
 * as a DNA sequence.  Lines beginning with '#' are treated as comments.
 *
 * @filepath: path to the text file
 * @buffer:   output buffer of at least MAX_SEQ_LEN bytes
 *
 * Returns 1 on success, 0 on failure (file not found, empty, read error).
 */
int load_sequence_from_file(const char *filepath, char *buffer);

/**
 * load_all_sequences_from_file
 * ----------------------------
 * Reads every non-comment, non-empty line from a FASTA-like file.
 * Lines beginning with '>' (FASTA headers) and '#' (comments) are skipped.
 * Each valid line is treated as an independent DNA sequence and processed
 * through a caller-supplied callback.
 *
 * @filepath : path to the text file
 * @callback : function pointer called for each sequence found.
 *             Receives: (sequence_string, sequence_index, user_data)
 * @user_data: arbitrary pointer forwarded to every callback invocation
 *
 * Returns the total number of sequences processed, or -1 on file error.
 */
typedef void (*seq_callback_t)(const char *seq, int index, void *user_data);

int load_all_sequences_from_file(const char     *filepath,
                                 seq_callback_t  callback,
                                 void           *user_data);

#endif /* FILE_HANDLER_H */
