from pipeline.trim_run import TrimRun
from pipeline.pipelinelogging import logger
import logging
import sys

class PipelineSteps:
    
    ####################################################################
    #
    #  1) Start Trimming
    # The RunTrim class should look for keys and primers
    # and trim the sequences if found.
    # Output should be (one for each lane_key):
    #   Fasta file: trimmed sequences
    #          names: 
    #   Fasta file: unique trimmed sequences
    #   Names file: 2 cols: 1st = id of unique, 
    #                       2nd = list of others it represents (include 1st),csl
    #   Text file:  list of deleted read_ids and reasons
    ####################################################################
    @classmethod
    def trim(cls, run):
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
            
        
    @classmethod    
    def chimera(cls, run):
        '''
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
        '''
    @classmethod
    def gast(cls, run):  
        pass
#        mygast = Gast(run, outputdir, args)
#        mygast.clustergast(new_lane_keys)
#        mygast.gast_cleanup(new_lane_keys)
#        mygast.gast2tax(new_lane_keys)

        
        
    @classmethod
    def vampsupload(cls, run):
        pass  
#        myvamps = Vamps(run, outputdir, args)
#        myvamps.info(new_lane_keys)
#        myvamps.projects(new_lane_keys)
#        myvamps.taxonomy(new_lane_keys)
#        myvamps.sequences(new_lane_keys)        
#        myvamps.exports(new_lane_keys)
       
    
    

    
    
    