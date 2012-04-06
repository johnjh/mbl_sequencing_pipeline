import os, sys
from subprocess import check_output
import constants as C
from pipeline.primer_utils import expand
from time import sleep
from pipeline.pipelinelogging import logger
from string import maketrans
import collections

base_complement_translator = maketrans("ACGT", "TGCA")

# the json expected files get loaded and parsed into Unicode strings
# but the asserts won't work comparing unicode to ascii so we need change them
# to plain strings
def convert_unicode_dictionary_to_str(data):
    if isinstance(data, unicode):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert_unicode_dictionary_to_str, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert_unicode_dictionary_to_str, data))
    else:
        return data        


########################################################    
def check_for_Ns(seq):
    """Doc string here.."""
    nCount = seq.count('N') 
    if( nCount > 0 ):
        return nCount
    else:
        return 0

def remove_runkey(seq,runkeys):
    """Doc string here.."""
    found_key = ''
    for key in runkeys:
        # find at first position BUG need to add TCAG to front to find key
        # is this a problem with the sff_extraction script??
        #key_plus = 'TCAG' + key
        #print seq.find(key)
        if (seq.find(key) == 0):
            found_key = key
            seq = seq[len(key):]
            break
#         elif(seq.find(key) > 0 and seq.find(key) < 5): 
#             found_key = key
#             seq = seq[seq.find(key) + len(key):]   
#             break   
        else:
            continue 
    return found_key, seq
    
def find_sequence_direction( direction = '' ):
    """Doc string here.."""
    seqFwd = 0
    seqRev = 0
    if(direction == 'F'):
        seqFwd = 1
    if(direction == "R"):
        seqRev = 1
        
    if ( not direction or (seqFwd + seqRev) == 0 ):
        return 0
    elif ( (seqFwd + seqRev) == 2 ):
        return "B"
    elif (seqFwd):
        return "F"
    else:
        return "R"

def check_for_quality(rawseq, trimseq, quality_scores):
    """Doc string here.."""
    start = rawseq.find(trimseq)
    end   = start + len(trimseq)
    
    scores = quality_scores[start:end]
    
    return sum(scores,0.0) / len(scores)

def revcomp(sequence):
    reversed = str(sequence[::-1])
    return reversed.translate(base_complement_translator)

def set_trim1():
    return (True,False,False,False)
 
def set_trim2():
    return (True,True,False,False)
 
def set_trim3():
    return (True,True,True,False)
 
def set_chim4():
    return (False,True,False,False)
    
def set_chim5():
    return (False,True,True,False)
 
def set_chim6():
    return (False,True,True,True)
 
def set_gast7():
    return (False,False,True,False)
 
def set_gast8():
    return (False,False,True,True)

def set_vamps9():
    return (False,False,False,True)
 
def set_all10():
    return (True,True,True,True)
    
options = {
        1 : set_trim1,
        2 : set_trim2,
        3 : set_trim3,
        4 : set_chim4,
        5 : set_chim5,
        6 : set_chim6,
        7 : set_gast7,
        8 : set_gast8,
        9 : set_vamps9,
        10 : set_all10,
}
 
def wait_for_cluster_to_finish(my_running_id_list):
    #print 'My IDs',running_id_list
    logger.debug('Max run time set to ' + str(C.cluster_max_wait) + ' seconds')
    logger.debug('These are my running qsub IDs ' + my_running_id_list)
    my_working_id_list = my_running_id_list

    counter =0

    sleep(C.cluster_initial_check_interval)
    
    while my_working_id_list:
    
    
        qstat_codes = get_qstat_id_list()
        if not qstat_codes['id']:
            #print 'No qstat ids'
            logger.debug("id list not found: may need to increase initial_interval if you haven't seen running ids.")
            return ('SUCCESS','id list not found','',)
        if 'Eqw' in qstat_codes['code']:
            logger.debug( "Check cluster: may have error code(s), but they may not be mine!")
        
        
        got_one = False
    
        #print 'working ids',my_working_id_list
        if my_working_id_list[0] in qstat_codes['id']:
            
            got_one = True
            name = qstat_codes['name'][qstat_codes['id'].index(my_working_id_list[0])]
            user = qstat_codes['user'][qstat_codes['id'].index(my_working_id_list[0])]
            code = qstat_codes['code'][qstat_codes['id'].index(my_working_id_list[0])]
            
            
            if code == 'Eqw':
                return ('FAIL','Found Eqw code',my_working_id_list[0])
            elif code == 'qw':
                logger.debug("id is still queued: " +  my_working_id_list[0] + " " + code)
            elif code == 'r':
                logger.debug("id is still running: " + my_working_id_list[0] + " " + code)
            else:
                logger.debug('Unknown qstat code ' + code)
        else:
            my_working_id_list = my_working_id_list[1:]
            logger.debug('id finished ' + my_working_id_list[0])
 
        if not my_working_id_list:
            return ('SUCCESS','not my_working_id_list','')
        #if not got_one:
            #print 'IN not got one',
        #    return ('SUCCESS','not got one','')
                
        sleep(C.cluster_check_interval)
        counter = counter + C.cluster_check_interval
        if counter >= C.cluster_max_wait:
            return ('FAIL','Max Time exceeded',C.cluster_max_wait)
    
    return ('FAIL','Unknown','Unknown')
    
def get_qstat_id_list():
    
    # ['139239', '0.55500', 'usearch', 'avoorhis', 'r', '01/22/2012', '09:00:39', 'all.q@grendel-07.bpcservers.pr', '1']
    # 1) id
    # 2) 
    # 3) name
    # 4) username
    # 5) code r=running, Ew=Error
    qstat_cmd = 'qstat'
    qstat_codes={}
    output = check_output(qstat_cmd)
    #print output
    output_list = output.strip().split("\n")[2:]
    qstat_codes['id'] = [n.split()[0] for n in output_list]
    qstat_codes['name'] = [n.split()[2] for n in output_list]
    qstat_codes['user'] = [n.split()[3] for n in output_list]
    qstat_codes['code'] = [n.split()[4] for n in output_list]
    #print 'Found IDs',qstat_ids
 
    
    
    return qstat_codes


def find_key(dic, val):
    """return the first key of dictionary dic given the value"""
    return [k for k, v in dic.iteritems() if v == val][0]
    
def mysort(uniques,names):
    """ Sorts the uniques using the uniques and names hashes:
    
    uniques[lane_tag][trimmed_sequence] = read_id
    names[lane_tag][read_id1] = [read_id1, read_id2, read_id3, read_id4]
    
    returns a list of tuples (read_id, count, sequence) highest to lowest
    """
    sorted_uniques = []
    
    # sorted_names should be list of ids with the highest number at the top
    sorted_names = sorted(names.items(), key=lambda x: len(x[1]), reverse=True)

    for n in sorted_names:
       
        seq = find_key(uniques, n[0])
        sorted_uniques.append( (n[0], len(n[1]),seq) )
    
    return sorted_uniques     






