import sys
import os
import re
import traceback
from optparse import OptionParser
from IndentedHelpFormatterWithNL import *

def error(msg):
    print("ERROR: %s\n" % msg)
    parser.print_help()
    sys.exit(-1)

def warn(msg):
    print("WARNING: %s" % msg)

def fatal(msg):
    print(msg)
    sys.exit(-1)


def join_fields(fields):
    return OUTPUT_FORMAT_DELIMITER.join(fields)


def split_line(line):
    if delimiter:
        return line.split(delimiter)
    else:
        return line.split()

def line_key(line):
    return chrpos_to_n( split_line(line) )


# Computes an integer key for this line. These keys can be used to sort the lines by reference
def chrpos_to_n(lsplit):
    # Get chr, pos from line

    chr_value, start_value = None, None # Init in case of error
    try:
        split1 = lsplit[0].split(":") # Get chr:start-stop out of the 1st column.
        chr_value = split1[0].lower().strip()
        split2 = split1[1].split("-")
        start_value = split2[0].lower().strip()
    except:
        sys.stderr.write("chrom: %s, start: %s. Couldn't parse line: %s \n" % (chr_value, start_value, line))
        raise

    # Covert them to N
    a = 0
    if sequence_build == "UCSC" and chr_value.count("_random"):
        chr_value = chr_value.replace("_random", "")
        a = 30 # Offset so that "random" chromosomes go last

    if sequence_build == "UCSC":
        chr_value = chr_value.replace("chrm", "chr0")
    else:
        chr_value = chr_value.replace("chrm", "chr25")

    chr_n = a + int(chr_value.replace("chrx", "chr23").replace("chry", "chr24").replace("chr",""))
    start_n = int(start_value)

    N = (chr_n * 10**11) + start_n

    #print("N: " + str(N) + " line: " +  line)
    return N

def is_valid_chrpos(line):
    try:
        # Compute the line key
        line_key(line)
        return True
    except Exception, e:
        #print(str(e))
        return False




# Init cmd-line args
description = """
This script takes a text-based tabular INPUT-FILE, validates it, and does whatever is necessary to convert it into the format required by the GenomicAnnotator.
The output format has:
 TODO write this
"""
parser = OptionParser( description=description, usage="usage: %prog [options] INPUT-FILE", formatter=IndentedHelpFormatterWithNL())
parser.add_option("-v", "--verbose", action="store_true", default=False,                help="Verbose.")
parser.add_option("-d", "--delimiter",                                                  help="The delimiter that separates values in a line of INPUT-FILE. Set to 'tab' to make it use tab [Default: spaces].")
parser.add_option("-c", "--location-columns", metavar="COLUMNS",                        help="""The (1-based) column number(s) of the columns in INPUT-FILE that contain coordinates. \n
For example, '-c 2,3' means column #2 and column #3 contain coordinate info. COLUMNS can be set to one, two, or three comma-separated numbers:\n
 1 number means column1 is of the form 'choromosome:position' or 'chromosome:start-stop'\n
 2 numbers means column1 = choromosome, column2 = position.\n
 3 numbers means column1 = choromosome, column2 = start position, column3 = stop position.""")

parser.add_option("-t", "--coordinate-type", dest="coordinates", metavar="COORD-TYPE", help="""Specifies the coordinate system of INPUT-FILE's chromosome/position column(s). COORD-TYPE can be:\n
  * ONE-BASED-HALF-OPEN   1-based half-open.\n
  * POSITIONAL            same as ONE-BASED-HALF-OPEN
  * ZERO-BASED-HALF-OPEN  0-based half-open.\n
  * OFFSET                same as ZERO-BASED-HALF-OPEN.\n
  * DONT-CHANGE           coordinates will be left as-is.\n
Note: This setting is used to convert all coordinates into 1-based half-open for the output.""")
parser.add_option("-s", "--sequence", dest="sequence_build", metavar="BUILD",           help="Sets the output file's reference build type to either UCSC or NCBI. This should be set based on what reference file will be used when running the GenomicAnnotator. UCSC builds can be specified as either 'hgXX' (eg. hg18) or 'UCSC'. NCBI builds can be specified as 'bXX' (eg. b36) or 'NCBI'. The build type determines chromosome order and naming convention (eg. 'chr1' or '1').")
#parser.add_option("-i", "--include-columns", dest="include_fields", metavar="COLUMNS",  help="A comma-separated listing of (1-based) column numbers of all columns to include in the outptut file. Any columns not in this list will be discarded.")
#parser.add_option("-e", "--exclude-columns", dest="exclude_fields", metavar="COLUMNS",  help="A comma-separated listing of (1-based) column numbers of the columns to include in the outptut file. Any columns not in this list will be discarded.")
parser.add_option("-o", "--output-filename",                                           help="Output file path [Default: %default]", default="stdout")


(options, args) = parser.parse_args()

#print(args) # List of positional args.
#print(options.output_filename)
#print(options.coordinates) # required
#print(options.location_columns)
#print(options.delimiter)
#print(options.verbose)


# Validate and process cmd-line args

verbose = options.verbose

if verbose:
    print("%s v0.9" % sys.argv[0])

delimiter = options.delimiter
if delimiter and delimiter.lower() == "tab":
    delimiter = "\t"

if len(args) < 1 or not os.access(args[0], os.R_OK):
    error("Requires a valid INPUT-FILE")
input_filename = args[0]


if options.coordinates == "POSITIONAL" or options.coordinates == "ONE-BASED-HALF-OPEN" or options.coordinates == "DONT-CHANGE":
    coords_type = 1
elif options.coordinates == "OFFSET" or options.coordinates == "ZERO-BASED-HALF-OPEN":
    coords_type = 0
else:
    if not options.coordinates:
        error("-t arg must be specified")
    else:
        error("Invalid -t value: %s" % str(options.coordinates))

if not options.location_columns:
        error("-c arg must be specified")

loc_columns = options.location_columns.split(",")
if len(loc_columns) < 1 or len(loc_columns) > 3:
    error("-c COLUMNS must specify a comma-separated list of between 1 and 3 numbers.")

#if verbose:
#    print("Parsed -c: " +  str(loc_columns))

try:
    chr_column = int(loc_columns[0]) - 1
    start_column = None
    stop_column = None
    pos_columns_sorted = [chr_column]
    if len(loc_columns) > 1:
        start_column = int(loc_columns[1]) - 1
        pos_columns_sorted += [start_column]

        if len(loc_columns) > 2:
            stop_column = int(loc_columns[2]) - 1
            pos_columns_sorted += [stop_column]

    pos_columns_sorted.sort(reverse=True)
except:
    error("-c COLUMNS - all elements in the comma-separated list must be integers.")

if (chr_column and chr_column < 0) or  (start_column and start_column < 0) or (stop_column and stop_column < 0):
    error("-c COLUMNS - all elements in the comma-separated list must be >= 1")


if not options.sequence_build:
    error("-s arg must be specified")

sequence_build = options.sequence_build.lower()
if sequence_build.startswith("b") or sequence_build == "ncbi":
    sequence_build = "NCBI"
elif sequence_build.startswith("hg") or sequence_build == "ucsc":
    sequence_build = "UCSC"
else:
    error("-s arg must be one of these: 'hgXX' (eg. hg18), 'UCSC', 'bXX' (eg. b36), 'NCBI'")
   


output_filename = options.output_filename
if output_filename and output_filename != "stdout" and output_filename != "-" and os.access(output_filename, os.F_OK) and not os.access(output_filename, os.W_OK):
    error("Unable to write to: %s" % str(options.output_filename))

if verbose:
    print("  Input file: " + input_filename)
    print("  Output file: " + output_filename)



# Commence processing
counter = 0
skipped_lines_counter = 0
header_line_found = False
header_fields = []
prepend_lines = []
data_lines = []
previous_n = -1 # Checks whether data is in order
need_to_sort = False

OUTPUT_FORMAT_DELIMITER = "\t"

for line in open(input_filename):
    line = line.strip()
    if counter % 1000000 == 0:
        print("Processed %d records" % counter )
        
    counter+=1
    if not header_line_found:
        if line.startswith("#") or line == "":
            prepend_lines += [line]
        else:
            header_line_found = True
            header_fields = split_line(line)

            # Remove all the existing positional columns, and make a new 1st column: 'chrpos'
            for c in pos_columns_sorted:
                if c >= len(header_fields):
                    error("Found only %d headers. -c arg is out of range." % len(header_fields) )

                header_fields.pop(c)
            header_fields.insert(0, "chrpos")

            if len(header_fields) < 2:
                error("Header appears to have only 1 column in it.")

            if verbose:
                print("Found header containing %d columns: [%s]. Changed it to: [%s]." % (len(split_line(line)), "  ".join(split_line(line)), "  ".join(header_fields)))

    else:
        # This is a data line
        line_fields = split_line(line)
        # Parse out the chrpos and make it the 1st column if its not already
        chrpos_value = ""
        if start_column:
            # There is more than 1 column of position info in this line
            # TODO error check line_fields[chr_column] somehow.
            try: start_int = int(line_fields[start_column])
            except: error("Line #%d, Column %d: start coordinate value '%s' is not an integer." % (counter, start_column, line_fields[start_column]))

            if coords_type == 0:
                start_int += 1 # Convert to 1-based coords

            chrpos_value = "%s:%d" % ( line_fields[chr_column], start_int )

            if stop_column:
                try: stop_int = int(line_fields[stop_column])
                except: error("Line #%d, Column %d: stop coordinate value '%s' is not an integer" % (counter, stop_column, line_fields[stop_column]))

                #if coords_type == 0:
                #    stop_int += 1
                # Converting to 1-based inclusive, so don't need to add 1 here after all

                if stop_int != start_int: # If they are equal, chr1:x is the same as chr1:x-x
                    chrpos_value += "-%d" % stop_int

            # Remove these old position columns and make the valid chrpos column the 1st
            for c in pos_columns_sorted:
                line_fields.pop(c)

            line_fields.insert(0, chrpos_value)

        else:
            # There is only 1 column of position info in this line
            if not re.match(".+\:\d+([-]\d+)?", line_fields[chr_column]):
                error("Line #%d: Invalid chrpos [%s] in column %d" % (counter, line_fields[chr_column], chr_column ))

            # Move it to be line 1
            chrpos_value = line_fields.pop(chr_column)
            line_fields.insert(0, chrpos_value)

        # Validate

        if len(line_fields) < len(header_fields):
            warn("Line #%d: Has %d columns [%s] while header has %d columns [%s]. The missing fields will be treated as empty." % (counter, len(line_fields), "  ".join(line_fields),  len(header_fields), "  ".join(header_fields), ))
            while len(line_fields) < len(header_fields):
                line_fields += [OUTPUT_FORMAT_DELIMITER + ""] # Append '' as filler. TODO - make this behavior a cmd-line switchable

        elif len(line_fields) > len(header_fields):
            warn("Line #%d: Has %d columns [%s] while header has %d columns [%s]. Skippping..." % (counter, len(line_fields), "  ".join(line_fields),  len(header_fields), "  ".join(header_fields), ))
            continue


        try:
            n = chrpos_to_n(line_fields)
            if not need_to_sort and n < previous_n:
                need_to_sort = True
                warn("Line %d is out of order. Will need to sort all lines." % counter)
            previous_n = n
        except Exception, e:
            warn("Couldn't parse line: " + "  ".join(line_fields) + ". " +str(e) + ". Skipping...")
            if verbose: traceback.print_exc()
            skipped_lines_counter += 1
            continue




        data_lines += [ join_fields(line_fields) ]

if verbose and skipped_lines_counter:
    print("Skipped %d / %d lines. (%f%%)" % (skipped_lines_counter, counter, skipped_lines_counter/float(counter)))
if need_to_sort:
    if verbose:
        print("Sorting %d lines..." % len(data_lines))

    data_lines.sort(key=line_key)


if verbose:
    print("Writing data to: " + output_filename)

# Write output file
if output_filename == "stdout" or output_filename == "-":
    output_file = sys.stdout
else:
    output_file = open(output_filename, "w+")

for line in prepend_lines + [ join_fields(header_fields) ]: 
    output_file.write(line + "\n")

for line in data_lines:
    if sequence_build == "NCBI" and line.lower().startswith("chr"):
        if line.lower().startswith("chrm"):
            output_file.write("MT" + line[4:] + "\n")
        else:
            output_file.write(line[3:] + "\n")
    else:
        output_file.write(line + "\n")

output_file.close()
