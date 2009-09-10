import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import itertools

class Status:
    def __init__(self, file, exists, size):
        self.file = file
        self.exists = exists
        self.size = size
        
        if not exists: self.status = "missing"
        if size == 0: self.status = "no-data"
        else: self.status = "exists: bytes=" + str(self.size)
    
    def __str__(self):
        return self.status
        
    def viewSize(self):
        return MergeBAMsUtils.greek(self.size)

class ComparedFiles:
    def __init__(self, file, status, localStat, ftpStat):
        self.file = file
        self.status = status
        self.localStat = localStat
        self.ftpStat = ftpStat
        
    def size(self):
        if self.localStat.size <> 0:
            return self.localStat.size
        if self.ftpStat.size <> 0:
            return self.ftpStat.size
        else:
            return 0

    def modTime(self):
        if self.localStat.exists:
            return os.path.getmtime(self.localStat.file)
        else:
            return 0

def modTimeStr(t):
    return time.strftime("%m/%d/%y", time.localtime(t))

from urlparse import urlparse
from ftplib import FTP

FTPSERVER = None

DEBUG = False

# from directories to lists of lines
CACHED_LIST = dict()
def getSizeForFile(dir, filename):
    global CACHED_LIST
    size = [0]    
    def captureSize(line, cache = True):
        #print line
        if cache: CACHED_LIST[dir].append(line)
        s = line.split()
        if len(s) >= 9 and s[8] == filename:
            size[0] = int(s[4])
            #print 'Found size', s, size

    if dir in CACHED_LIST:
        #print 'cached is', CACHED_LIST[dir]
        map( lambda l: captureSize(l, False), CACHED_LIST[dir] )
    else:
        FTPSERVER.cwd(dir)
        CACHED_LIST[dir] = list()
        result = FTPSERVER.retrlines('LIST', captureSize)

    return size[0]

def ftpStatus( ftpPath ):
    if DEBUG: print 'ftpPath', ftpPath
    dir, filename = os.path.split(ftpPath)
    if DEBUG: print 'listing', dir

    try:
        size = getSizeForFile(dir, filename)   
    except:
        #print 'failing...'
        size = 0
#    finally:
#        pass
        #print 'FTPSERVER', FTPSERVER
        #FTPSERVER.quit()

    if DEBUG: print '  result was', size
    return Status( ftpPath, size <> 0, size )

def localStatus(file):
    exists = os.path.exists(file)
    size = 0
    if exists: size = os.path.getsize(file)
    return Status(file, exists, int(size) )        

def validateFile(relPath, localRoot, ftpRoot):
    localPath = os.path.join(root, relPath)
    ftpPath = os.path.join(ftpRoot, relPath)
    
    # check the local file
    if DEBUG: print 'Checking', relPath
    localStat = localStatus(localPath)
    ftpStat = ftpStatus(ftpPath)
    if DEBUG: print ' local status is', localStat
    if DEBUG: print ' ftp status is  ', ftpStat
    compared = compareFileStatus(localStat, ftpStat)
    
    if not OPTIONS.quiet: 
        print 'STATUS %20s for %s ' % (compared.status, relPath)
    return compared

import MergeBAMsUtils

import time
def compareFileStatus(localStat, ftpStat):
    if localStat.exists:
        if ftpStat.exists:
            if localStat.size == ftpStat.size:
                status = 'in-sync'
            else:
                status = 'size-mismatch'
        else:
            status = 'unknown-local-file'
    else:
        if ftpStat.exists:
            status = 'local-file-missing'
        else:
            status = 'orphaned-file'
    
    return ComparedFiles(localStat.file, status, localStat, ftpStat)


import re
def filesInLocalPath(root, subdir):
    regex = re.compile(".*\.(bam|bai)$")
    localFiles = set()
    
    if subdir <> None:
       for fullroot, dirs, files in os.walk(os.path.join(root, subdir)):
           for file in filter( regex.match, files ):
               fullpath = os.path.join(fullroot, file)
               path = fullpath.split(root)[1]
               #print 'adding relpath=', path, 'fullpath=', fullpath
               localFiles.add(path)
               if OPTIONS.maxLocalFiles <> None and len(localFiles) > OPTIONS.maxLocalFiles: return localFiles
    return localFiles

def readAlignmentIndex(file):
    files = set()
    if file <> None:
       for line in open(file):
           files.add(line.split()[0])
           files.add(line.split()[4])
    return files

if __name__ == "__main__":
    usage = "usage: %prog -l and/or -a root ftpRoot"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--local", dest="scanLocal",
                        type='string', default=None,
                        help="If provided, checks all of the local files against the archive")
    parser.add_option("-a", "--alignmentIndex", dest="alignmentIndex",
                        type='string', default=None,
                        help="If provided, checks all of the files in the alignment.index in the archive")
    parser.add_option("-m", "--maxLocal", dest="maxLocalFiles",
                        type='int', default=None,
                        help="If provided, maximum number of files in the local archive to examine")
    parser.add_option("-q", "--quiet", dest="quiet",
                        action='store_true', default=False,
                        help="If provided, prints out the individual status of all files")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    root, ftpRoot = args

    ftpParsed = urlparse(ftpRoot)
    FTPSERVER = FTP(ftpParsed[1])
    FTPSERVER.login()

    results = dict()
    for file in itertools.chain(readAlignmentIndex(OPTIONS.alignmentIndex), filesInLocalPath(root, OPTIONS.scanLocal )):
        #print line
        #bas = line.split()[6]
        if file not in results:
            compared = validateFile( file, root, ftpParsed[2] )
            results[file] = compared
            #localIndex
        
    print 'SUMMARY: Total files examined', len(results)
    for status in ['in-sync', 'size-mismatch', 'unknown-local-file', 'local-file-missing', 'orphaned-file']:
        print ''.join(['-'] * 80)
        filesOfStatus = filter(lambda x: x.status == status, results.itervalues())
        n = len(filesOfStatus)
        print 'SUMMARY: %s' % ( status )
        print 'SUMMARY: files                    %d (%.2f%% of total)' % ( n, n * 100.0 / len(results))

        if n > 0:
           fileSizes = MergeBAMsUtils.greek(reduce(operator.__add__, map( ComparedFiles.size, filesOfStatus ), 0 ))
           mostRecentMod = apply(max, map( ComparedFiles.modTime, filesOfStatus ))
           if mostRecentMod > 0:
               modTime = modTimeStr(mostRecentMod)
           else:
               modTime = "N/A"
               
           print 'SUMMARY: total size               %s' % ( fileSizes )
           print 'SUMMARY: last modification time   %s' % ( modTime )
