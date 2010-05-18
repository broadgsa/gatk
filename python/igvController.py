# Echo server program
import socket, re
from os.path import join
from time import sleep

from farm_commands2 import *
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import faiReader
import math
import shutil
import string
import time
from madPipelineUtils import *

HOST = 'vm0e0-052.broadinstitute.org'       # Symbolic name meaning the local host
PORT = 60151                                # Arbitrary non-privileged port
LOCAL_DIR = "/Users/depristo/Desktop/IGV_screenshots"
SLEEP_TIME = 1

def main():
    global OPTIONS
    usage = """usage: %prog [options] sites
Automatically captures IGV PNG screenshots at each site in the file sites (of the form chrX:start-stop or chrX:pos or a VCF file) by connecting to 
an IGV session.  See http://www.broadinstitute.org/igv/?q=PortCommands.  Make sure you enable ports in the IGV preferences"""
    parser = OptionParser(usage=usage)
    parser.add_option("-w", "--wait", dest="wait",
                        action='store_true', default=False,
                        help="If provided, instead of taking screenshots we will prompt the user on the command line to press return and then jump to the new location")
    parser.add_option("", "--host", dest="host",
                        type='string', default=HOST,
                        help="The host running the port enabled IGV server")
    parser.add_option("-d", "--dir", dest="dir",
                        type='string', default=LOCAL_DIR,
                        help="The local directory on the host machine that screenshots should be written to")
    parser.add_option("-s", "--sort", dest="sort",
                        type='string', default="base",
                        help="The sort order for the bases, currently can be one of base, position, strand, quality, sample, and readGroup.")
    parser.add_option("-p", "--prefix", dest="prefix",
                        type='string', default="",
                        help="A prefix to add before the contig name before sending the command to IGV.  Useful for dealing with b36 -> hg18 issues")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    sites = args[0]
    #It might ignore the first line..
    print "Be sure to turn on ports"
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    def sendCommand(cmd):
        s.send(cmd)
        print cmd.strip(), '=>', s.recv(512).strip()

    s.connect((OPTIONS.host, PORT))
    sendCommand("snapshotDirectory " + OPTIONS.dir + "\n")
    c = 0
    for line in open(sites):
        parts = line.split()
        
        if sites.find(".vcf") != -1 and parts[0][0] != "#":
            site = OPTIONS.prefix + ':'.join(parts[0:2])
        else:
            site = parts[0]
        
        print site
        c+=1
        sendCommand("goto %s\n" % site)
        sendCommand("sort " + OPTIONS.sort + "\n")
        if OPTIONS.wait: 
            raw_input("Enter To Continue") 
        else:
            print 'sleep', SLEEP_TIME, 'secs' 
            time.sleep(SLEEP_TIME)
            sendCommand("snapshot\n") # %s.png\n" % re.sub("-","_",re.sub(':', '_', site)))
    print c

            
if __name__ == "__main__":
    main()


