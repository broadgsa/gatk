import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import itertools
from urlparse import urlparse
from ftplib import FTP
import MergeBAMsUtils
import time
import re
import hashlib


FTPSERVER = None
DEBUG = False

CACHED_LIST = dict() # from directories to lists of lines

class Status:
    def __init__(self, file, exists, size):
        self.file = file
        self.exists = exists
        self._size = size
        
        if not exists: self.status = "missing"
        if size == 0: self.status = "no-data"
        else: self.status = "exists: bytes=" + str(self.size())
    
    def __str__(self):
        return self.file + " " + self.status

    def size(self):
        return self._size
        
    def viewSize(self):
        return MergeBAMsUtils.greek(self.size())
        

def md5(file):
    m = hashlib.md5()
    for line in open(file):
        m.update(line)
    return m.hexdigest()

class ComparedFiles:
    def __init__(self, file, status, localStat, ftpStat):
        self.file = file
        self.status = status
        self.localStat = localStat
        self.ftpStat = ftpStat
        
    def size(self):
        if self.localStat.size() <> 0:
            return self.localStat.size()
        if self.ftpStat.size() <> 0:
            return self.ftpStat.size()
        else:
            return 0

    def modTime(self):
        if self.localStat.exists:
            return os.path.getmtime(self.localStat.file)
        else:
            return 0

def modTimeStr(t):
    if t == 0:
        return 'N/A'
    else:
        return time.strftime("%m/%d/%y", time.localtime(int(t)))

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

def fetchFtpFile( file ):
    filename = os.path.split(file)[1]
    destFile = filename + '.fetched.' + date.today().strftime("%m_%d_%y")
    #print 'destFile', destFile
    fd = open(destFile, 'w')
    result = FTPSERVER.retrbinary('RETR ' + file, lambda x: fd.write(x))
    fd.close()
    #print "done"
    return Status(destFile, True, os.path.getsize(destFile))

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

def compareFileStatus(localStat, ftpStat):
    if DEBUG: print 'comparing', localStat, ftpStat
    if localStat.exists:
        if ftpStat.exists:
            if localStat.size() == ftpStat.size():
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


def filesInLocalPath(root, subdir):
    regex = re.compile(".*\.(bam|bai)$")
    localFiles = set()
    
    if subdir <> None:
       for fullroot, dirs, files in os.walk(os.path.join(root, subdir)):
           for file in filter( regex.match, files ):
               #if file <> "NA12761.SLX.WUGSC.Mosaik.SRP000033.2009_08.bam.bai":
               #    continue   
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
           parts = line.split()
           files.add(parts[0])
           if len(parts) > 4: # we have an index, we're not unmapped
               files.add(line.split()[4])
    return files

def compareAlignmentIndices(remoteAlignmentIndex, alignmentIndex):
    if remoteAlignmentIndex <> None and alignmentIndex <> None:
        printHeaderSep()
        print 'Comparing remote and local alignment indices: '
        remotePath = os.path.join(ftpParsed[2], remoteAlignmentIndex)
        remoteAlignmentIndexFile = fetchFtpFile( remotePath )
        print '  Fetched', remotePath, 'to', remoteAlignmentIndexFile.file
        raImd5 = md5(remoteAlignmentIndexFile.file)
        laImd5 = md5(alignmentIndex)
        print '  md5s: local=%s remote=%s' % (raImd5, laImd5)
        if raImd5 <> laImd5:
            print '  [FAIL] -- alignment indices do not have the same hash!'
            sys.exit(1)
        else:
            print '  [PASS] -- alignment indices are the same'

def displayChangeLog( changelog ):
    if changelog <> None:
        printHeaderSep()
        print 'Displaying remote changelog for examination '
        remotePath = os.path.join(ftpParsed[2], changelog)
        remoteChangeLog = fetchFtpFile( remotePath )
        print '  Fetched', remotePath, 'to', remoteChangeLog.file
        
        print 
        for line in itertools.islice(open(remoteChangeLog.file), 20):
            print 'CHANGELOG', line,

def printHeaderSep():
    print
    print ''.join(['-'] * 80)

def sortByName(files):
	return sorted(files, key=lambda x: x.file)

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
    parser.add_option("-M", "--maxFiles", dest="maxFiles",
                        type='int', default=None,
                        help="If provided, maximum number of files in the local archive to examine")
    parser.add_option("-q", "--quiet", dest="quiet",
                        action='store_true', default=False,
                        help="If provided, prints out the individual status of all files")
    parser.add_option("-i", "--remoteAlignmentIndex", dest="remoteAlignmentIndex",
                        type='string', default=None,
                        help="relative path to the FTP's alignment.index file for comparison")
    parser.add_option("-c", "--remoteChangeLog", dest="remoteChangeLog",
                        type='string', default=None,
                        help="relative path to the FTP's CHANGELOG file for display")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    root, ftpRoot = args

    ftpParsed = urlparse(ftpRoot)
    FTPSERVER = FTP(ftpParsed[1])
    FTPSERVER.login()

    displayChangeLog(OPTIONS.remoteChangeLog)
    compareAlignmentIndices(OPTIONS.remoteAlignmentIndex, OPTIONS.alignmentIndex)

    results = dict()
    for file in itertools.chain(readAlignmentIndex(OPTIONS.alignmentIndex), filesInLocalPath(root, OPTIONS.scanLocal )):
        #print line
        #bas = line.split()[6]
        if file not in results:
            compared = validateFile( file, root, ftpParsed[2] )
            results[file] = compared
            #localIndex
        if OPTIONS.maxFiles != None and len(results) > OPTIONS.maxFiles:
            break

    printHeaderSep()
    print 'SUMMARY: Total files examined', len(results)
    for status in ['in-sync', 'size-mismatch', 'unknown-local-file', 'local-file-missing', 'orphaned-file']:
        printHeaderSep()
        filesOfStatus = sortByName(filter(lambda x: x.status == status, results.itervalues()))
        n = len(filesOfStatus)
        print 'SUMMARY: %s' % ( status )
        print 'SUMMARY: Files                    %d (%.2f%% of total)' % ( n, n * 100.0 / len(results))
        
        statusForFileListing = ['size-mismatch', 'local-file-missing', 'orphaned-file']
        maxFilesToList = 20
        if status in statusForFileListing:
            print 'SUMMARY: listing the first', min(maxFilesToList, n), 'of', n
            for file in itertools.islice(filesOfStatus, maxFilesToList):
                if status == 'size-mismatch':
                    print 'SUMMARY: File: ftp=%d bytes local=%d bytes %12s %s' % ( file.ftpStat.size(), file.localStat.size() , modTimeStr(file.modTime()), file.file)
                else:
                    print 'SUMMARY: File: %8s %12s %s' % ( MergeBAMsUtils.greek(file.size()), modTimeStr(file.modTime()), file.file)
        if n > 0:
            fileSizes = MergeBAMsUtils.greek(reduce(operator.__add__, map( ComparedFiles.size, filesOfStatus ), 0 ))
            mostRecentMod = modTimeStr(apply(max, map( ComparedFiles.modTime, filesOfStatus ) + [0]))
                
            print 'SUMMARY: total size               %s' % ( fileSizes )
            print 'SUMMARY: last modification time   %s' % ( mostRecentMod )

    totalSize = MergeBAMsUtils.greek(reduce(operator.__add__, map( ComparedFiles.size, results.itervalues() )))
    print '#### TOTAL PROJECT SIZE:          %s' % ( totalSize )
