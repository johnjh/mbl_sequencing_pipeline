#!/usr/local/www/vamps/software/python/bin/python

##!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import os
from stat import * # ST_SIZE etc
import sys
import shutil
import types
from time import sleep
from pipeline.utils import *
from pipeline.sample import Sample
from pipeline.runconfig import RunConfig
from pipeline.run import Run
from pipeline.chimera import Chimera
from pipeline.gast import Gast
from pipeline.pipelinelogging import logger
from pipeline.trim_run import TrimRun
import logging
import argparse
    

TRIM_STEP = "trim"
CHIMERA_STEP = "chimera"
GAST_STEP = "gast"
VAMPSUPLOAD = "vampsupload"

existing_steps = [TRIM_STEP, CHIMERA_STEP, GAST_STEP, VAMPSUPLOAD]

# perform trim step
# TrimRun.trimrun() does all the work of looping over each input file and sequence in each file
# all the stats are kept in the trimrun object
#
# when complete...write out the datafiles for the most part on a lane/runkey basis
#
def trim(run):
    mytrim = TrimRun(run)        
    # pass True to write out the straight fasta file of all trimmed non-deleted seqs
    # Remember: this is before chimera checking
    trim_codes = mytrim.trimrun(True)
    if trim_codes[0] == 'SUCCESS':
        new_lane_keys = trim_codes[2]
        logger.debug("Trimming finished successfully")
        run.status_file_h.write("TRIM SUCCESS\n")
        run.status_file_h.write("new_lane_keys="+','.join(new_lane_keys)+"\n")            
        # write_data_files: names, unique and abund files
        mytrim.write_data_files(new_lane_keys)
    else:
        print ("Trimming Failed")
        run.status_file_h.write("TRIM ERROR: "+trim_codes[1]+" "+trim_codes[2]+"\n")
        sys.exit()
        
def chimera(run):
    chimera_cluster_ids = [] 
    logger.debug("Starting Chimera Checker")
    mychimera = Chimera(run, outputdir, args)
    c_den    = mychimera.chimera_denovo(new_lane_keys)
    if c_den[0] == 'SUCCESS':
        chimera_cluster_ids += c_den[2]
        chimera_code='PASS'
    elif c_den[0] == 'NOREGION':
        chimera_code='NOREGION'
    elif c_den[0] == 'FAIL':
        chimera_code = 'FAIL'
    else:
        chimera_code='FAIL'
    
    c_ref    = mychimera.chimera_reference(new_lane_keys)
    
    if c_ref[0] == 'SUCCESS':
        chimera_cluster_ids += c_ref[2]
        chimera_code='PASS'
    elif c_ref[0] == 'NOREGION':
        chimera_code = 'NOREGION'
    elif c_ref[0] == 'FAIL':
        chimera_code='FAIL'
    else:
        chimera_code='FAIL'
    
    #print chimera_cluster_ids
    
    if chimera_code == 'PASS':  
        
        chimera_cluster_code = wait_for_cluster_to_finish(chimera_cluster_ids,args) 
        if chimera_cluster_code[0] == 'SUCCESS':
            logger.info("Chimera checking finished successfully")
            status_file_h.write("CHIMERA SUCCESS\n")
            
            
        else:
            logger.info("Chimera checking Failed")
            status_file_h.write("CHIMERA ERROR: "+str(chimera_cluster_code[1])+" "+str(chimera_cluster_code[2])+"\n")
            sys.exit()
            
    elif chimera_code == 'NOREGION':
        logger.info("No regions found that need chimera checking")
        status_file_h.write("CHIMERA CHECK NOT NEEDED\n")
        
    elif chimera_code == 'FAIL':
        logger.info("Chimera checking Failed")
        status_file_h.write("CHIMERA ERROR: \n")
        sys.exit()
    else:
        logger.info("Chimera checking Failed")
        status_file_h.write("CHIMERA ERROR: \n")
        sys.exit()
    sleep(2)   
    if  go_chimera and chimera_code == 'PASS' and  chimera_cluster_code[0] == 'SUCCESS':
        mychimera.write_chimeras_to_deleted_file(new_lane_keys)
        # should also recreate fasta
        # then read chimera files and place (or replace) any chimeric read_id
        # into the deleted file.
        
        mymblutils = MBLPipelineFastaUtils(outputdir=outputdir,lane_keys=new_lane_keys)
        
        # write new cleaned files that remove chimera if apropriate
        # these are in fasta_mbl_pipeline.py
        # the cleaned file are renamed to the original name:
        # lane_key.unique.fa
        # lane_key.trimmed.fa
        # lane_key.names        -- 
        # lane_key.abund.fa     -- this file is for the uclust chimera script
        # lane_key.deleted.txt  -- no change in this file
        # THE ORDER IS IMPORTANT HERE:
        mymblutils.write_clean_fasta_file()
        mymblutils.write_clean_names_file()
        mymblutils.write_clean_uniques_file()
        mymblutils.write_clean_abundance_file()
        # write keys file for each lane_key - same fields as db table? for easy writing
        # write primers file for each lane_key
 
        
        # Write new clean files to the database
        # rawseq table not used
        # trimseq
        # runkeys
        # primers
        # run primers
        mymblutils.write_clean_files_to_database()

def gast(run):  
    pass
#        mygast = Gast(run, outputdir, args)
#        mygast.clustergast(new_lane_keys)
#        mygast.gast_cleanup(new_lane_keys)
#        mygast.gast2tax(new_lane_keys)

    
def vampsupload(run):
    pass  
#        myvamps = Vamps(run, outputdir, args)
#        myvamps.info(new_lane_keys)
#        myvamps.projects(new_lane_keys)
#        myvamps.taxonomy(new_lane_keys)
#        myvamps.sequences(new_lane_keys)        
#        myvamps.exports(new_lane_keys)
        
def process(run, steps):
    # create output directory:
    requested_steps = steps.split(",")            
    
    if not os.path.exists(run.output_dir):
        logger.debug("Creating directory: "+run.output_dir)
        os.makedirs(run.output_dir)      
        run.status_file_h = open(run.status_file,"a")
    elif("trim" in requested_steps):
        logger.debug("Removing and recreating directory: "+run.output_dir)
        shutil.rmtree(run.output_dir)        
        os.makedirs(run.output_dir)
        run.status_file_h = open(run.status_file,"a")
        run.status_file_h.write("TRIM STARTING\n")
        logger.debug("Output directory exists: overwriting")
    else:
        logger.debug("Keeping directory: "+run.output_dir)
        # reading and writing 'r+'
        run.status_file_h = open(run.status_file,"r+")
        logger.debug("found status file",run.status_file)
        for line in run.status_file_h.readlines():
            line = line.strip()
            if line.split('=')[0] == 'new_lane_keys':
                new_lane_keys = line.split('=')[1].split(',')
                #logger.debug('new_lane_keys',new_lane_keys
                if type(new_lane_keys) is not types.ListType:
                    run.status_file_h.write("READ ERROR: No lane_keys found\n")
                    logger.debug("READ ERROR: No lane_keys found\n")
                    sys.exit()
                    
    # loop through official list...this way we execute the
    # users requested steps in the correct order                
    for step in requested_steps:
        if step not in existing_steps:
            print "Invalid processing step: " + step
            sys.exit()
        else:
            # call the method in here
            step_method = globals()[step]
            step_method(run)
    
    run.status_file_h.close()
