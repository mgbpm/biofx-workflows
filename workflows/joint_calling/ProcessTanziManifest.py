#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Ingest cram manifest tsv file and produce batched output
# for sampletracker ingest and text output for the jointCalling
# terra workflow input fields.

import csv
import logging
from datetime import datetime

# check inputs
try:
    MANIFEST_TSV_FILE = (
        "/Users/AJR87_1/Documents/WDL_DEV/JointGenotyping/Tanzi-CRAM-manisfet.tsv"
    )
    # MANIFEST_TSV_FILE = sys.argv[1]
except:
    raise Exception("Usage: python ProcessTanziManifest.py tanzi_manifest.tsv")


# open input file and ingest tanzi data
print("Processing tanzi manifest file: " + MANIFEST_TSV_FILE)
logging.info("Processing tanzi manifest file: " + MANIFEST_TSV_FILE)
with open(MANIFEST_TSV_FILE, newline="") as csvfile:
    # skip first line
    next(csvfile)
    spamreader = csv.DictReader(csvfile, delimiter="\t", quotechar="|")
    spamreader.fieldnames = (
        "SIZE",
        "DATE",
        "WASABI_CRAM_PATH",
        "FILENAME",
        "ID_1",
        "ID_2",
    )
    line_number = 0
    records_removed = 0
    # define top level dictionary
    tanzi_cram_dict = {}
    print(spamreader)
    # iterate through all tanzi manifest records
    for row in spamreader:
        # incriment line num count
        line_number += 1
        # define each record dictionary
        record_dict = {}
        # populate record
        record_dict["SIZE"] = row["SIZE"]
        # add date as datetime object
        try:
            newdate = datetime.strptime(str(row["DATE"]), "%m/%d/%y")
        except:
            newdate = str(row["DATE"])
        record_dict["DATE"] = newdate
        record_dict["WASABI_CRAM_PATH"] = row["WASABI_CRAM_PATH"]
        record_dict["FILENAME"] = row["FILENAME"]
        record_dict["ID_1"] = row["ID_1"]
        record_dict["ID_2"] = row["ID_2"]

        # attempt to add record to dictionary
        if record_dict["ID_1"] in tanzi_cram_dict.keys():
            print("Attempting to overwtire: " + record_dict["ID_1"])
            print("Removing both records.")
            # decided to remove duplicate cram entries
            del tanzi_cram_dict[record_dict["ID_1"]]
            records_removed += 1
            # should we only use crams with a later date?
            """
            date1 = tanzi_cram_dict["ID_1"]["DATE"]
            date2 = record_dict["DATE"]
            if date1 < date2:
                print(str(date1) + " is earlier than " + str(date2))
            elif date1 > date2:
                print(str(date1) + " is later than  " + str(date2))
            else:
                print(str(date1) + " is the same as  " + str(date2))
            """
        else:
            # add record dict to tanzi dict by ID_1
            tanzi_cram_dict[record_dict["ID_1"]] = record_dict
print("Removed " + str(records_removed) + " duplicated records.")

# sample tracker column names
sample_sheet_header = "Subject ID,Sample ID,Project,Data Location"
# example data location s3://prod-biobank-cram-2023-1/biobank/cram/tanzi/genome/220816_AHLYJNDSX3/10000055_ds.ded30d1c5ac8456cb00258cc0d33a590

# set initial states
records_per_batch = 50
current_batch = 1
line_number = 0
record_number_in_batch = 1
# open initial batch csv
csv_name = "tanzi_cram_sampletracker_load_batch_" + str(current_batch) + ".csv"
OUTPUT_CSV = open(csv_name, "wt")
OUTPUT_CSV.write(sample_sheet_header)
# iterate through all tanzi crams writing load files in batches of 50
for key in tanzi_cram_dict.keys():
    # track total number of records processed
    line_number += 1
    # build csv entry for record
    csv_string = (
        str(tanzi_cram_dict[key]["ID_1"]).strip(' "')
        + ","
        + str(tanzi_cram_dict[key]["ID_2"]).strip(' "')
        + ","
        + "BiobankCramToVcf"
        + ","
        + "s3://prod-biobank-cram-2023-1/biobank/cram/"
        + str(tanzi_cram_dict[key]["WASABI_CRAM_PATH"]).strip(' "')
    )

    if record_number_in_batch > records_per_batch:
        # close current output file
        OUTPUT_CSV.close()
        # increase batch number by 1
        current_batch += 1
        # update output csv to new batch name
        csv_name = "tanzi_cram_sampletracker_load_batch_" + str(current_batch) + ".csv"
        OUTPUT_CSV = open(csv_name, "wt")
        # reset record_number_in_batch = 1
        record_number_in_batch = 1
        # add current record to new output file
        OUTPUT_CSV.write(sample_sheet_header)
        OUTPUT_CSV.write("\n" + csv_string)
    else:
        # write csv entry for row
        OUTPUT_CSV.write("\n" + csv_string)
        record_number_in_batch += 1
