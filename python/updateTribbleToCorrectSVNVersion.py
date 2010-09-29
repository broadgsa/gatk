# ====================================================================================================
# update the tribble code to the appropriate version, given the current GATK version.
# this python script looks up the base version of the GATK, checks the date, and cross-references this
# against the Tribble version.
# ====================================================================================================

import os, subprocess, re 
from subprocess import Popen, PIPE, STDOUT
from datetime import datetime, date, time

# test that we can use the SVN tool
print "checking that we can use svn and svnversion commands..."
outputSVN 			= Popen(["svn", 		"--version"], stdout=PIPE).communicate()[0]
outputSVNVersion 	= Popen(["svnversion", 	"--version"], stdout=PIPE).communicate()[0]

# match the regex
m = re.search('version\s+(\d+\.\d+\.\d+)',outputSVN)
m2 = re.search('version\s+(\d+\.\d+\.\d+)',outputSVNVersion)

if ((not m) or (not m2)):
	raise Exception("Unable to find svn and svnversion commands") 
	
if (m.group(1) != m.group(1)):
	raise Exception('Unable to match versions for svn (version ' + m.group(1) + ') and svnversion (version ' + m2.group(1) + ') commands') 

# check that we're in the base GATK directory
outputSVN = Popen(["svn","info"], stdout=PIPE).communicate()[0]
m = re.search('URL:\s+(.+)\s+',outputSVN)

print "checking that we're in the base directory of a SVN check-out of the GATK..."
if (not m):
	raise Exception("Unable to find current working URL for this directory in SVN; make sure your in the base GATK directory when running this script")
if (m.group(1) != "https://svn/Sting/trunk"):
	raise Exception("Not in the correct working directory; we expect to be in a directory pointing to the URL https://svn/Sting/trunk")

# now get the version and date for the GATK checkout
print "getting the GATK check-out version and date..."
dateMatch = re.search('Last Changed Date:\s+(\S+)\s+(\S+)\s+',outputSVN)
versionMatch = re.search('Last Changed Rev:\s+(.+)\s+(.+)\s+',outputSVN)
if ((not dateMatch) or (not versionMatch)):
	raise Exception("Unable to match either the version or the check-out date from the svn -info output: " + outputSVN)	

dt = datetime.strptime(dateMatch.group(1)+" "+dateMatch.group(2), "%Y-%m-%d %H:%M:%S")
print "		The date of the current GATK check-out is = " +  dt.strftime("%A, %d. %B %Y %I:%M%p")

# now look through the tribble logs for the last version before the current GATK date
# the log entries look like: r213 | hanna | 2010-09-22 11:28:16 -0400 (Wed, 22 Sep 2010) | 2 lines
outputSVN = Popen(["svn","log","--quiet","tribble"], stdout=PIPE).communicate()[0]
m = re.findall("r(\d+)\s\|\s\S+\s\|\s(\S+\s\S+)",outputSVN)
diff = None
bestRev = 0
for match in m:
	rev = match
	tribbleDate = datetime.strptime(match[1] , "%Y-%m-%d %H:%M:%S")
	# print tribbleDate.strftime("%A, %d. %B %Y %I:%M%p")
	if (dt > tribbleDate):
		if (diff == None):
			diff = dt - tribbleDate	
			bestRev = match[0]
		elif (dt - tribbleDate < diff):
			diff = dt - tribbleDate
			bestRev = match[0]

if (bestRev == 0):
	raise Exception("Unable to find correct revision that predates the current GATK checkout...failing")
	
# now update the tribble directory to the found revision
print "attempting to update Tribble to the correct version of r" + bestRev
print ""
print "SVN update output:"
print Popen(["svn","update","-r","r"+str(bestRev),"tribble"], stdout=PIPE).communicate()[0]


