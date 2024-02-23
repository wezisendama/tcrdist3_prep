#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 08:37:20 2024

@author: wezi
"""

import pandas as pd
import subprocess

# Read in metadata
sarcoid_metadata = pd.read_csv('/home/wezi/work/sarcoidosis/sarcoidhealthymetadata.txt', sep = '\t')

# Pull filenames and patient IDs from metadata
multi_organ_subset = sarcoid_metadata.loc[sarcoid_metadata['source_name'] == "Multi-organ_PBMC"]
multi_organ_filenames = multi_organ_subset['Filename']

# Copy files in question to tcrdist_test folder
for i in multi_organ_filenames.index.values:
    copyname = "/home/wezi/work/sarcoidosis/" + multi_organ_filenames[i]
    copycommand = "cp " + copyname + " /home/wezi/work/tcrdist_test"
    subprocess.call(copycommand, shell=True)

# Add patient IDs as "subject" column to each file
for i in multi_organ_filenames.index.values:
    edit_target_filename = "/home/wezi/work/tcrdist_test/" + multi_organ_filenames[i]
    target_dataframe = pd.read_csv(edit_target_filename, sep = '\t')
    subject_value = multi_organ_subset["Patient"][i]
    target_dataframe = target_dataframe.assign(subject = subject_value)
    target_dataframe.to_csv(edit_target_filename, sep = '\t', index=False)
    
# Stick files together in one table
full_tcrtable = pd.DataFrame()

for i in multi_organ_filenames.index.values:
    table_filename = "/home/wezi/work/tcrdist_test/" + multi_organ_filenames[i]
    loaded_table = pd.read_csv(table_filename, sep = "\t")
    full_tcrtable = pd.concat([full_tcrtable, loaded_table], ignore_index = True)

# Remove the out-of-frame CDR3aa rows
filtered_tcrtable = full_tcrtable.loc[full_tcrtable['CDR3aa'] != "out_of_frame"]

# Strip out non-CDR3b data
filtered_tcrtable = filtered_tcrtable.loc[filtered_tcrtable['V'].str.contains('TRB')]

# Rename columns to ones recognised by tcrdist3 (cdr3_b_aa, cdr3_b_nt, v_b_gene, j_b_gene)
filtered_tcrtable = filtered_tcrtable.rename(columns = {
    '#count':'count',
    'CDR3nt':'cdr3_b_nt',
    'CDR3aa':'cdr3_b_aa',
    'V':'v_b_gene',
    'D':'d_b_gene',
    'J':'j_b_gene'})

# Write to new .tsv file
filtered_tcrtable.to_csv("/home/wezi/work/tcrdist_test/multi_organ_PBMC_tcrb.tsv", sep = '\t', index=False)