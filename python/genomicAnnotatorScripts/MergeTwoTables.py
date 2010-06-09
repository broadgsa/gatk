import sys
import os


def print_help():
    sys.stderr.write("\n" + os.path.split(sys.argv[0])[1] + " [file1] [file2] \n" + \
        "     Takes two tab-delimited tables and merges them, so that the output is sorted by genomic position.\n" + \
        "     Both input files must be in AnnotatorInputTable format (http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator#Data_Formats),\n" + \
        "     and must have identical headers.\n")

def read_header(file_obj):       
    for line in file_obj:
        line = line[0:-1] # Remove trailing \n
        if line.strip() != "" and line[0] != "#":
            return  line.split("\t")
    else:
        raise Exception, "Reached the end of the file without finding the header"
            



if len(sys.argv) != 3:
    print_help()
    sys.exit(0)

try:
    file1 = open(sys.argv[1])
    header1 = read_header(file1)
except Exception, e:
    sys.stderr.write("ERROR: While reading header from file \"" + sys.argv[1] + "\": " + str(e) + "\n")
    sys.exit(0)

try:
    file2 = open(sys.argv[2])
    header2 = read_header(file2)
except e:
    sys.stderr.write("ERROR: While reading header from file \"" + sys.argv[1] + "\": " + str(e) + "\n")
    sys.exit(0)


if len(header1) != len(header2):
    sys.stderr.write("ERROR: The two files' headers are of different lengths: \n" + str(header1) + "\n" + str(header2) + "\n")
    sys.exit(0)

if header1 != header2:
    sys.stderr.write("WARNING: The two files' headers are of different lengths: \nHeader1: " + str(header1) + "\nHeader2: " + str(header2) + "\nUsing header1.\n")
print("\t".join(header1))


def get_chrom(line):
    idx = line.find(":")
    if idx == -1:
        raise Exception, "Invalid file format. No ':' found in line, so couldn't parse chromosome name: " + line
    chrom = line[0:idx]
    return chrom
    
# Computes a sort key for chromosome names (UCSC order)
def compute_chrom_key(chr_value):
    a = 0
    chr_value = chr_value.lower()

    if chr_value.count("_random"):
        chr_value = chr_value.replace("_random", "")
        a = 30 # Offset so that "random" chromosomes go last                                                                                                                                                                               

    chr_value = chr_value.replace("chrm", "chr0").replace("chrx", "chr23").replace("chry", "chr24")
    chr_value = chr_value.replace("chr","")
    return a + int(chr_value) + 1

def compute_sort_key(line):
    idx = line.find('\t')
    if idx == -1:
        chrpos = line
    else:        
        chrpos = line[0:idx]

    idx = chrpos.find(":")
    if idx == -1:
        return chrpos
    chrom = chrpos[0:idx]
    pos = chrpos[idx+1:]

    idx = pos.find("-")
    if idx == -1:
        return int(pos)
    else:
        start = pos[0:idx]
        end = pos[idx+1:]
        return int(start)


def read_line(file_obj):
    try:
        line = file_obj.next()[0:-1] # Remove \n
        key = compute_sort_key(line)
        return (line, key)
    except StopIteration:
        return (None, None)
    except Exception, e:
        sys.stderr.write("ERROR: While reading file \"" + sys.argv[1] + "\": " + str(e) + "\n")
        sys.exit(0)


# Read the 1st lines of each file        
line1, key1 = read_line(file1)
line2, key2 = read_line(file2)


# Do a merge sort
while line1 != None or line2 != None: # Iterate over each chromosome
    # Compute the next chromosome
    if line1 != None and line2 != None:
        chrom1 = get_chrom(line1)
        chrom2 = get_chrom(line2)
        if compute_chrom_key(chrom1) < compute_chrom_key(chrom2):
            current_chrom = chrom1
        else:
            current_chrom = chrom2
    elif line1 != None:
        current_chrom = get_chrom(line1)
    elif line2 != None:
        current_chrom = get_chrom(line2)

    # Iterate over lines for that chromosome
    while line1 != None and line2 != None and get_chrom(line1) == current_chrom and get_chrom(line2) == current_chrom: 
        
        if key2 > key1:
            print(line1)
            #print("line1 -  key1: " + str(key1)  + " key2: " + str(key2)) 
            used_line1 = True            
        else:
            #print("line2 -  key1: " + str(key1)  + " key2: " + str(key2)) 
            print(line2)
            used_line1 = False
            
        if used_line1:
            line1, key1 = read_line(file1)
        else:
            line2, key2 = read_line(file2)

        

    # At this point, either line1 or line2 will == None
            
    while line1 != None and get_chrom(line1) == current_chrom:
        print(line1)
        line1, key1 = read_line(file1)

    while line2 != None and get_chrom(line2) == current_chrom:
        print(line2)
        line2, key2 = read_line(file2)
        

            

