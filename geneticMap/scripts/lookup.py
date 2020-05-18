#!/usr/bin/env python

from argparse import ArgumentParser
import sys, signal
from sys import stderr
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def usage():
    import sys
    print('lookup(<keys file>, <key column>, <search file>, <search column>, <return column (optional)>')
    print('keys_file: file that contains the items to look up\
key_column (optional): the (1-based) column index that contains the keys in the keys file. default=1\
return_column (optional): the (1-based) column from which to return. Default is "True/False" for membership test\
search_file: the file to search keys for and return a value (test of membership, or another column in the aray)\
search_column: the (1-based) column in which to search for keys)')
    sys.exit(1)

def tx2gene(s):
    '''convert transcripts to systematic gene names'''
    try:
        int(s[-1])
        return s
    except ValueError:
        s = s[:-1]
        return tx2gene(s)

def lookup(keys_file, search_file, search_column, search_delim, search_delim_idx, key_column, key_delim, key_delim_idx, return_column, tx2g):
    """generic lookup generator function. yields line from search file with result appended"""

    key_lines = {}
    with open(keys_file, 'U') as f:
        total = set()
        duplicate_keys = set()
        for l in f:
            if l.startswith('#'):
                l = f.next()
            else:
                col = l.strip('\n').split()
                key = col[key_column-1]
                #print key                
                if key_delim:
                    try:
                        key = key.split(key_delim)[key_delim_idx-1]
                    except IndexError:
                        stderr.write('abort. supplied key column index failed for col\n{}\n'.format(key))
                        usage()
                if tx2g: key = tx2gene(key)
                total.add(key)
                if not return_column:
                    key_lines[key] = ''
                else:
                    try:
                        r_col = int(return_column-1)
                    except TypeError:
                        stderr.write('abort. if supplied, return column must be an integer\n')
                        usage()
                    if key in key_lines:
                        key_lines[key] = col[return_column-1]
                        duplicate_keys.add(key)
                    else:                          
                        key_lines[key] = col[return_column-1]
        if duplicate_keys != set():
            err = 'warning, {} of {} keys were not unique, last ocurring value was used\n'.format(len(duplicate_keys), len(total))
            sys.stderr.write(err)
              
    with open(search_file, 'U') as f:
        for l in f:
            if l.startswith('#'):
                l = f.next()
            else:
                col = l.strip('\n').split()
                check = col[search_column-1]
                #print check
                if search_delim:
                    try:
                        check = check.split(search_delim)[search_delim_idx-1]
                    except IndexError:
                        stderr.write('abort. supplied search column index failed for col\n{}\n'.format(check))
                        usage()

                if tx2g: check = tx2gene(check)
                if not return_column:
                    if check in key_lines:
                        col.append('True')
                        yield('\t'.join(col))
                    else:
                        col.append('False')
                        yield('\t'.join(col))
                else:
                    if check in key_lines:
                        col.append(key_lines[check])
                        yield('\t'.join(col))
                    else:
                        col.append('NA')
                        yield('\t'.join(col))

def main():

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    parser = ArgumentParser(
        description = 'Lookup key from file1 column<kc> in file2 column<sc>, print line + True/False for membership OR column<rc> from file2 to stdout. All indexes 1-based. Splits on space.')
    parser.add_argument('-k', dest = 'keys_file', type = str, help = 'target file with value to return (file1)')
    parser.add_argument('-s', dest = 'search_file', type = str, help = 'query file to which lookup value will be appended (file2)')
    parser.add_argument('-sc', dest = 'search_column', type = int, help = 'column index in query (file2), default 1', default = 1)
    parser.add_argument('-sd', dest = 'search_column_delim', type = str, help = 'optional, delimiter to split search column before matching match', default = False)
    parser.add_argument('-sdi', dest = 'search_column_delim_index', type = int, help = 'optional, index for delimiter-split search column', default = False)
    parser.add_argument('-kc', dest = 'key_column', type = int, help = 'column index in target file containing keys (file1), default 1',default = 1)
    parser.add_argument('-kd', dest = 'key_column_delim', type = str, help = 'optional, delimiter to split key column before match', default = False)
    parser.add_argument('-kdi', dest = 'key_column_delim_index', type = int, help = 'optional, index for delimiter-split key column', default = False)
    parser.add_argument('-rc', dest = 'return_column', type = int, help = 'optional, column index in target file to append to search line, default is True/False for membership test', default = False)
    parser.add_argument('-tx2g', dest = 'convert_transcripts', type = bool, help = 'optional, if search/target keys are transcripts, convert to gene names', default = False)
    args = parser.parse_args()

    for l in lookup(args.keys_file, args.search_file, args.search_column, args.search_column_delim, args.search_column_delim_index, args.key_column, args.key_column_delim, args.key_column_delim_index, args.return_column, args.convert_transcripts):
        print(l)

main()
