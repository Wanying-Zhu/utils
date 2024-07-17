def merge_files(lst_input_fn, output_fn, header=True):
    '''
    merge_files(lst_input_fn, output_fn, header=True)
        This function reads in files with names in lst_input_fn, merge lines and outptut to output_fn.
        If input files have header line, then the output file will retain the same header (only once)
    
    Parameters
    ----------
        - lst_input_fn: a list or list-like collection of file names to be merged
        - output_fn: File name of output file to write merged result
        - header: Boolean, True if input files have header line. Only write the header line once in output file
                  There can only be one line as the header line in each file.
    
    Return
    ----------
        - Mereged content of all input files (list_input_fn) is written to output file with output_fn as file name
    '''
    count = 1 # check if this is the first file
    fh_output = open(output_fn, 'w')

    for fn in lst_input_fn:
        with open(fn) as fh:
            line = fh.readline().strip()
            
            if header: # If input files have header line
                if count == 1: # Write the header line from the first file
                    fh_output.write(line+'\n')

                # Skip the first line if this is not the first file, read in an extra line directlry
                line = fh.readline().strip()
            else: # If no header line, then move forward to the next step
                pass

            while line != '':
                # Write to output file
                fh_output.write(line+'\n')

                count += 1
                line = fh.readline().strip()

    fh_output.close()
    
    