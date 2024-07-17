# This function removes white space at the start and end of each line
# Returns a cleaned file
def remove_start_end_white_space(fn):
    with open(fn ,'r') as fh:
        line = fh.readline().strip()
        # Write the first line
        with open(fn+'.whitespace_removed', 'w+') as output_fh:
            output_fh.write(line)

        # Write the rest lines
        line = fh.readline().strip()
        with open(fn+'.whitespace_removed', 'a') as output_fh:
            while line != '':
                output_fh.write('\n' + line)
                line = fh.readline().strip()



