#!/usr/bin/env Python

import sys

debug = False

print("Merging:")

# open files for reading
vcfInputFiles = []
vcfOutputFile = ""
for arg in sys.argv:
    if( arg.endswith("Fields.py") ):
        # do nothing
        continue
    else:
        if ( arg.startswith("I=") ):
            input1 = arg.strip().split("I=")[1]
            input2 = input1.split(",")
            for vcfInput in input2:
                print(vcfInput)
                vcfInputFiles.append(vcfInput)
        elif ( arg.startswith("O=") ):
            vcfOutputFile=open(arg.strip().split("O=")[1],'w')
        else:
            print("Unsupported argument: "+arg)
            sys.exit()

lines = 0

# the proper (albeit slower) way to do this is to read the info fields
# into a dict and then iterate through the keys

infodict = dict()

for inFile in vcfInputFiles:
    # for each file
    for line in open(inFile).readlines():
        # for each line
        if ( line.startswith("#") ):
            # do nothing
            continue
        else:
            spline = line.strip().split()
            chrompos = spline[0]+":"+spline[1]
            info = spline[7]
            if ( chrompos in infodict ):
                if ( info == "." ):
                    # do nothing
                    continue
                else:
                    curInfo = infodict[chrompos]
                    if ( curInfo == "." or curInfo == ""):
                        newInfo = info
                    else:
                        # now we need to parse the fields
                        curinfofields = set(curInfo.split(";"))
                        newinfofields = set(info.split(";"))
                        curinfofields.update(newinfofields)
                        newInfo = ";".join(curinfofields)
                    
                    infodict[chrompos] = newInfo
                    #print(newInfo)
            else:
                infodict[chrompos] = info

# dictionary has been constructed; now just iterate through the first vcf
# and update the info field, printing to the output file

if ( debug ):
    for key in infodict:
        print(infodict[key])

for line in open(vcfInputFiles[0]).readlines():
    line = line.strip()
    if ( line.startswith("#") ):
        # this is a header
        vcfOutputFile.write(line+"\n")
    else:
        # this is a real line
        spline = line.split()
        outstring = ""
        fieldno = 0
        for field in spline:
            #print("fieldno="+str(fieldno))
            if ( fieldno == 7 ):
                outstring = outstring+infodict[spline[0]+":"+spline[1]]+"\t"
                fieldno = fieldno + 1
            else:
                outstring = outstring+field+"\t"
                fieldno = fieldno + 1
        # we wrote an extra \t, replace the last one with an \n
        outstring = outstring.strip()+"\n"
        vcfOutputFile.write(outstring)


# and we're done
