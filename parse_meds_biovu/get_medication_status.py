'''
Compare the medication file (/data100t1/share/BioVU/phenos/official_release_0619/Below_MEDS_20190531.csv) with a list of drugs.
Find if each individual used any of the drug in the list (1) or not (0)
Output file has 3 columns: ID, RECORD_Of_DRUG_USAGE, MATCHED
- RECORD_Of_DRUG_USAGE: number of records matched in the lookup list

Example runs:
python get_medication_status.py \
--med_record_fn ./example_data/meds_record.csv \
--drug_name_col DRUG_NAME \
--med_list ./meds_list/lipid_lowering.txt \
--output_path ./output \
--output_prefix test_run


python /data100t1/home/wanying/BioVU/202405_PAGE_PRS/code/utils/get_medication_status.py \
--med_record_fn /data100t1/share/BioVU/phenos/official_release_0619/Below_MEDS_20190531.csv \
--drug_name_col DRUG_NAME \
--med_list /data100t1/home/wanying/BioVU/202405_PAGE_PRS/supporting_files/medication_ref/lipid_lowering.txt \
--output_path /data100t1/home/wanying/BioVU/202405_PAGE_PRS/outputs/medication_status \
--output_prefix output

'''

import os
import argparse
import logging
import datetime

# #################### Helper funcitons ####################
def setup_log(fn_log, mode='w'):
    '''
    Print log message to console and write to a log file.
    Will overwrite existing log file by default
    Params:
    - fn_log: name of the log file
    - mode: writing mode. Change mode='a' for appending
    '''
    # f string is not fully compatible with logging, so use %s for string formatting
    logging.root.handlers = [] # Remove potential handler set up by others (especially in google colab)
    logging.basicConfig(level=logging.DEBUG,
                        handlers=[logging.FileHandler(filename=fn_log, mode=mode),
                                  logging.StreamHandler()], format='%(message)s')

def process_args(log_args=True):
    '''
    Process arguments
    - log_args: If true, save arguments into log file
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--med_record_fn', help='File name of the medication record. Must have IDs in the 1st column, and a column of drug names', type=str,
                        default='/data100t1/share/BioVU/phenos/official_release_0619/Below_MEDS_20190531.csv')
    parser.add_argument('--drug_name_col', type=str, help='Column name of drug name to compare with the list of drug names', default='DRUG_NAME')
    parser.add_argument('--med_list', type=str, help='File name of a list of drug names to check. One drug name per line',
                        default='/data100t1/home/wanying/BioVU/202405_PAGE_PRS/supporting_files/medication_ref/lipid_lowering.txt')
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='output')
    args = parser.parse_args()
    
    if not os.path.isdir(args.output_path): # Create output folder if not exists
        os.makedirs(args.output_path)

    # Record arguments used
    if log_args:
        fn_log = os.path.join(args.output_path, args.output_prefix+'.log')
        setup_log(fn_log, mode='w')
        
    # Record script used
    cmd_used = 'python ' + os.path.basename(__file__)
    logging.info('# Arguments used:')
    for arg in vars(args):
        cmd_used += f' --{arg} {getattr(args, arg)}'
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)
    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    return args
        
def process_line(line, sep, meds, indx_drug_name):
    '''
    Process a single line
    Params:
    - line: a single line from the medical record
    - sep: separator of the line
    - meds: an array of medicines to check against
    - indx_drug_name: index of drug name column
    Return:
    - Sample id of current line
    - drug_mark: 1 if drug_name from current line can be found meds, else 0
    - Name of the drug fund in the list
    '''
    # Split the line to get the columns
    lst_line = line.split(sep)
    sample_id = lst_line[0]
    drug_name = lst_line[indx_drug_name]
        
    # Check if the drug_name is in the drug list
    if drug_name.upper() in meds:
        return sample_id, 1, drug_name.upper()
    else:
        return sample_id, 0, ''
    
# #################### Process arguments ####################
args = process_args()

# Load medication info
# Read drug names from the medication list into a set
print('# Load medicine reference file')
meds = set()
with open(args.med_list, 'r') as ref_file:
    for line in ref_file:
        meds.add(line.strip().upper()) # Convert to uppercase for comparisons

output_fh = open(f'{args.output_path}/{args.output_prefix}.txt', 'w')

if args.med_record_fn.endswith('.csv'): # Determine delimiter of the med record file to split
    sep = ','
else:
    sep = '\t'

start = datetime.datetime.now()

# errors='replace': to avoid UnicodeDecodeError error
with open(args.med_record_fn, mode='r', errors='replace') as records_file:
    # Read the header, find index of the "DRUG_NAME"
    line = records_file.readline().strip()
    lst_line = line.split(sep)
    indx_drug_name = lst_line.index(args.drug_name_col)
    sample_id_col = lst_line[0]
    output_fh.write(f'{sample_id_col}\tRECORD_Of_DRUG_USAGE\tMATCHED\n')
    
    # Process rest of the lines
    # prev_sample_id: previous sample id
    # n_records: track the number of records from the sample individual
    prev_sample_id, drug_mark, count, n_records = None, 0, 0, 0
    matched_drugs = set() # store drug names matched in the list
    for line in records_file:
        count += 1
        line = line.strip()
        if not line:
            continue # skip empty lines
            
        # Current sample_id and drug_mark
        cur_sample_id, cur_drug_mark, drug = process_line(line=line, sep=sep, meds=meds,
                                                          indx_drug_name=indx_drug_name)
        n_records += 1
        if prev_sample_id is not None and prev_sample_id!=cur_sample_id: # Start a new sample
            # Write the previous result to file
            output_fh.write(prev_sample_id + '\t' + str(drug_mark) + '\t' + ','.join(matched_drugs) + '\n')
            output_fh.flush()
            # Re-process current line with the correct prev_drug_mark
            matched_drugs = set()
            # Reset drug_mark
            prev_sample_id, drug_mark, drug = process_line(line=line, sep=sep, meds=meds, indx_drug_name=indx_drug_name)
            n_records = 1 # Reset record
        else:
            # Set current sample id and drug mark to "previous"
            prev_sample_id = cur_sample_id
            drug_mark += cur_drug_mark
        if drug != '': matched_drugs.add(drug)
        if count%10000==0:
            print(f'\r# Process lines: {count}        ', end='', flush=True)
    print(f'\r# Process lines: {count}        ', flush=True)
    # Write result of the last sample
    if n_records == 1:
        # If last sample only has one record, need to re-process the line using the correct drug_mark to get correct cur_*
        cur_sample_id, drug_mark, cur_drug_mark = process_line(line=line, sep=sep, meds=meds,
                                                    indx_drug_name=indx_drug_name)
        
    output_fh.write(cur_sample_id + '\t' + str(drug_mark) + '\t' + ','.join(matched_drugs)+ '\n')
    
duration = datetime.datetime.now() - start
if duration.total_seconds() > 60:
    logging.info('# Finished in %.4f min' % (duration.total_seconds()/60))
else:
    logging.info('# Finished in %.4f s' % duration.total_seconds())
logging.info('# DONE')
output_fh.close()
