import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *

# Init cmd-line args
description = """
This script takes a text-based tabular INPUT-FILE, validates it, and converts it into the format expected by the GenomicAnnotator.
More details can found here: http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator#Tabular_Data_Format
"""
parser = OptionParser( description=description, usage="usage: %prog [options] INPUT-FILE", formatter=IndentedHelpFormatterWithNL())
parser.add_option("-l", "--location-columns", metavar="COLUMNS",                        help="""The (1-based) column number(s) of the columns in INPUT-FILE that contain coordinates. \n
For example, '-l 2,3' means column #2 and column #3 contain coordinate info. COLUMNS can be set to one, two, or three comma-separated numbers:\n
 1 number means column1 is of the form 'choromosome:position' or 'chromosome:start-stop'\n
 2 numbers means column1 = choromosome, column2 = position.\n
 3 numbers means column1 = choromosome, column2 = start position, column3 = stop position.""")

parser.add_option("-c", "--input-coords", dest="coordinates", metavar="COORD-TYPE", help="""Specifies the coordinate system of INPUT-FILE's chromosome/position column(s). COORD-TYPE can be:\n
  * ONE-BASED-HALF-OPEN   1-based half-open.\n
  * POSITIONAL            same as ONE-BASED-HALF-OPEN
  * ZERO-BASED-HALF-OPEN  0-based half-open.\n
  * OFFSET                same as ZERO-BASED-HALF-OPEN.\n
  * DONT-CHANGE           coordinates will be left as-is.\n
Note: This setting is used to convert all coordinates into 1-based half-open for the output.""")
parser.add_option("-t", "--output-style", dest="sequence_build", metavar="BUILD",           help="Sets the output file's reference build type to either UCSC or NCBI. This should be set based on what reference file will be used when running the GenomicAnnotator. UCSC builds can be specified as either 'hgXX' (eg. hg18) or 'UCSC'. NCBI builds can be specified as 'bXX' (eg. b36) or 'NCBI'. The build type determines chromosome order and naming convention (eg. 'chr1' or '1').")
#parser.add_option("-i", "--include-columns", dest="include_fields", metavar="COLUMNS",  help="A comma-separated listing of (1-based) column numbers of all columns to include in the outptut file. Any columns not in this list will be discarded.")
#parser.add_option("-e", "--exclude-columns", dest="exclude_fields", metavar="COLUMNS",  help="A comma-separated listing of (1-based) column numbers of the columns to include in the outptut file. Any columns not in this list will be discarded.")

group = OptionGroup(parser, "Optional Args", " ")

group.add_option("-o", "--output-filename",                                                                  help="Output file path [Default: %default]", default="stdout")
group.add_option("-r", "--haplotype-reference-column", metavar="COLUMN", dest="haplotype_reference_column",  help="1-based column number of the column to use as haplotypeReference. Specifying this will rename the column to 'haplotypeReference' in the header.")
group.add_option("-a", "--haplotype-alternate-column", metavar="COLUMN", dest="haplotype_alternate_column",  help="1-based column number of the column to use as haplotypeAlternate. Specifying this will rename the column to 'haplotypeAlternate' in the header.")
group.add_option("-s", "--haplotype-strand-column", metavar="COLUMN", dest="haplotype_strand_column",        help="1-based column number of the haplotypeStrand. Specifying this will rename the column to 'haplotypeStrand' in the header.")
group.add_option("-k", "--keep-original-columns", action="store_true", default=False, dest="keep_copy",      help="This flag makes it so that the columns passed to -l, -r, -a, and -s args are not removed when their contents is used to generate the special columns (eg. 'chrpos', 'haplotypeReference', etc..).")
group.add_option("-m", "--other-start-columns", metavar="COLUMNS", dest="other_start_columns",               help="Comma-separated list of 1 or more column numbers (1-based) representing other columns that contain start coordinates and need to be converted from the coordinate system specified by -c. For example, the refGene table has coordinates for cdsStart which need to be converted along with the chromosome, txStart, and txEnd columns.")
#group.add_option("-n", "--other-end-columns", metavar="COLUMNS", dest="other_stop_columns",                  help="Comma-separated list of 1 or more column numbers (1-based) representing other columns that contain end coordinates and need to be converted from the coordinate system specified by -c. For example, the refGene table has coordinates for cdsEnd which need to be converted along with the chromosome, txStart, and txEnd columns.")
group.add_option("-v", "--verbose", action="store_true", default=False,                                      help="Verbose.")
group.add_option("-d", "--delimiter",                                                                        help="The delimiter that separates values in a line of INPUT-FILE. Set to 'tab' to make it use tab [Default: spaces].")



parser.add_option_group(group)

(options, args) = parser.parse_args()



def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    #parser.print_help()
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
        start_n = long(start_value)
        stop_n = start_n
        if len(split2) > 1:
            stop_value = split2[1].lower().strip()
            stop_n = long(stop_value)
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

    chr_n = a + int(chr_value.replace("chrx", "chr23").replace("chry", "chr24").replace("chr","")) + 1

    N = (chr_n * 10L**23) + (start_n * 10L**11) + stop_n # Combine chr, start, stop into a single numeric key for sorting

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
        error("-c arg must be specified")
    else:
        error("Invalid -c value: %s" % str(options.coordinates))

if not options.location_columns:
        error("-l arg must be specified")

loc_columns = options.location_columns.split(",")
if len(loc_columns) < 1 or len(loc_columns) > 3:
    error("-l COLUMNS must specify a comma-separated list of between 1 and 3 numbers.")

#if verbose:
#    print("Parsed -c: " +  str(loc_columns))

try:
    chr_column = int(loc_columns[0]) - 1
    start_column = None
    stop_column = None
    columns_to_be_moved_in_order = [chr_column]
    if len(loc_columns) > 1:
        start_column = long(loc_columns[1]) - 1
        columns_to_be_moved_in_order += [start_column]

        if len(loc_columns) > 2:
            stop_column = long(loc_columns[2]) - 1
            columns_to_be_moved_in_order += [stop_column]

except:
    error("-l COLUMNS - all elements in the comma-separated list must be integers.")

if (chr_column and chr_column < 0) or  (start_column and start_column < 0) or (stop_column and stop_column < 0):
    error("-l COLUMNS - all elements in the comma-separated list must be >= 1")


if not options.sequence_build:
    error("-t arg must be specified")

sequence_build = options.sequence_build.lower()
if sequence_build.startswith("b") or sequence_build == "ncbi":
    sequence_build = "NCBI"
elif sequence_build.startswith("hg") or sequence_build == "ucsc":
    sequence_build = "UCSC"
else:
    error("-t arg must be one of these: 'hgXX' (eg. hg18), 'UCSC', 'bXX' (eg. b36), 'NCBI'")


haplotype_reference_column = None
if options.haplotype_reference_column:
    try:
        haplotype_reference_column = int(options.haplotype_reference_column) - 1
        columns_to_be_moved_in_order += [haplotype_reference_column]
    except:
        error("-r arg must be an integer")


haplotype_alternate_column = None
if options.haplotype_alternate_column:
    try:
        haplotype_alternate_column = int(options.haplotype_alternate_column) - 1
        columns_to_be_moved_in_order += [haplotype_alternate_column]
    except:
        error("-a arg must be an integer")

haplotype_strand_column = None
if options.haplotype_strand_column:
    try:
        haplotype_strand_column = int(options.haplotype_strand_column) - 1
        columns_to_be_moved_in_order += [haplotype_strand_column]
    except:
        error("-s arg must be an integer")

keep_copy = options.keep_copy

output_filename = options.output_filename
if output_filename and output_filename != "stdout" and output_filename != "-" and os.access(output_filename, os.F_OK) and not os.access(output_filename, os.W_OK):
    error("Unable to write to: %s" % str(options.output_filename))


other_start_columns = []
if options.other_start_columns:
    try:
        for c in options.other_start_columns.split(","):
            other_start_columns += [int(c) - 1]
    except:
        error("-m COLUMNS - all elements in the comma-separated list must be integers.")

if verbose:
    print("  Input file: " + input_filename)
    print("  Output file: " + output_filename)


columns_to_be_moved_in_order.sort(reverse=True)



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

            if len(header_fields) < 2:
                error("Header appears to have only 1 column in it.")

            # Remove all the existing positional columns, and make a new 1st column: 'chrpos'
            if haplotype_reference_column and haplotype_reference_column >= len(header_fields):
                error("Found only %d headers. -r arg is out of range." % len(header_fields) )

            if haplotype_alternate_column and haplotype_alternate_column >= len(header_fields):
                error("Found only %d headers. -a arg is out of range." % len(header_fields) )

            if haplotype_strand_column and haplotype_strand_column >= len(header_fields):
                error("Found only %d headers. -s arg is out of range." % len(header_fields) )

            for c in columns_to_be_moved_in_order:
                if c >= len(header_fields):
                    error("Found only %d headers. Column %d (specified as part of -l or another COLUMN arg) is out of range." % (len(header_fields), c ) )

                if not keep_copy:
                    header_fields.pop(c)



            # Rename columns to haplotypeReference, haplotypeAlternate, haplotypeStrand, and move them so that they are the 2nd,3rd,4th columns:
            if haplotype_strand_column:
                header_fields.insert(0, "haplotypeStrand")

            if haplotype_alternate_column:
                header_fields.insert(0, "haplotypeAlternate")

            if haplotype_reference_column:
                header_fields.insert(0, "haplotypeReference")

            header_fields.insert(0, "chrpos")

            if verbose:
                print("Found header containing %d columns: [%s]. Changed it to: [%s]." % (len(split_line(line)), "  ".join(split_line(line)), "  ".join(header_fields)))

    else:
        # This is a data line
        line_fields = split_line(line)

        # Get the haplotype ref/alt/strand values
        if haplotype_reference_column:
            haplotype_reference_value = line_fields[haplotype_reference_column]

        if haplotype_alternate_column:
            haplotype_alternate_value = line_fields[haplotype_alternate_column]

        if haplotype_strand_column:
            haplotype_strand_value = line_fields[haplotype_strand_column]


        # Compute the chrpos value from the chr,start,stop columns
        chrpos_value = ""
        if start_column:
            # There is more than 1 column of position info in this line
            # TODO error check line_fields[chr_column] somehow.
            try: start_int = long(line_fields[start_column])
            except: error("Line #%d, Column %d: start coordinate value '%s' is not an integer." % (counter, start_column, line_fields[start_column]))

            if coords_type == 0:
                start_int += 1 # Convert to 1-based coords
                line_fields[start_column] = str(start_int) # Change the original column in case keep_copy is True

            chrpos_value = "%s:%d" % ( line_fields[chr_column], start_int )

            if stop_column:
                try: stop_int = long(line_fields[stop_column])
                except: error("Line #%d, Column %d: stop coordinate value '%s' is not an integer" % (counter, stop_column, line_fields[stop_column]))

                #if coords_type == 0:
                #    stop_int += 1
                # Converting to 1-based inclusive, so don't need to add 1 here after all

                if stop_int != start_int: # If they are equal, chr1:x is the same as chr1:x-x
                    chrpos_value += "-%d" % stop_int
        else:
            # There is only 1 column of position info in this line
            if not re.match(".+\:\d+([-]\d+)?", line_fields[chr_column]):
                error("Line #%d: Invalid chrpos [%s] in column %d" % (counter, line_fields[chr_column], chr_column ))

        # Handle the -m arg
        if other_start_columns and coords_type == 0:
            for c in other_start_columns:
                if c >= len(line_fields):
                    error("Line #%d: Found only %d fields. -m arg is out of range." % (counter, len(line_fields)) )

                try:
                    converted_coords_string = ""
                    for coord in line_fields[c].split(","):
                        if coord.strip() == "":
                            continue

                        if len(converted_coords_string) > 0:
                            converted_coords_string += ","
                        converted_coords_string += str(long(coord) + 1)

                    line_fields[c] = converted_coords_string
                except:
                    error( "Line #%d: Processing -m %s arg. Couldn't parse coordinates in column %d: [%s]." % (counter, str(other_start_columns), c, line_fields[c] ) )


        # Move the columns around as needed (eg. so that chrpos is in the 1th column and hap ref/alt/strand are 2nd,3rd,4th):
        if not keep_copy:
            for c in columns_to_be_moved_in_order:
                line_fields.pop(c)


        if haplotype_strand_column:
            line_fields.insert(3, haplotype_strand_value)

        if haplotype_alternate_column:
            line_fields.insert(2, haplotype_alternate_value)

        if haplotype_reference_column:
            line_fields.insert(1, haplotype_reference_value)

        line_fields.insert(0, chrpos_value)



        # Validate
        if len(line_fields) < len(header_fields):
            warn("Line #%d: Has %d columns [%s] while header has %d columns [%s]. The missing fields will be treated as empty." % (counter, len(line_fields), "  ".join(line_fields),  len(header_fields), "  ".join(header_fields), ))
            while len(line_fields) < len(header_fields):
                line_fields += [OUTPUT_FORMAT_DELIMITER + ""] # Append '' as filler. TODO - make this behavior a cmd-line switchable

        elif len(line_fields) > len(header_fields):
            warn("Line #%d: Has %d columns [%s] while header has %d columns [%s].     Skipping..." % (counter, len(line_fields), "  ".join(line_fields),  len(header_fields), "  ".join(header_fields), ))
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
