#
# Generates BAM lists from Excel and TSV files provided by project managers.  Suitable for input into the pre-QC metrics generation
# script.
#
# To run:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/lib/poi-3.8-beta3.jar:$STING_HOME/lib/poi-ooxml-3.8-beta3.jar:$STING_HOME/lib/poi-ooxml-schemas-3.8-beta3.jar:$STING_HOME/lib/xmlbeans-2.3.0.jar:$STING_HOME/lib/dom4j-1.6.1.jar \
#     parse_pm_input.py <input file.{xls|xlsx|txt|tsv}> > <bam.list>
#
from java.io import FileInputStream
from org.apache.poi.ss.usermodel import Row,Sheet,Workbook,WorkbookFactory

import os,sys

base_path = '/seq/picard_aggregation/%s/%s'

def excel_reader(filename):
    wb = WorkbookFactory.create(FileInputStream(filename));
    for sheet_number in range(wb.getNumberOfSheets()):
        project_column = None
        sample_column = None

        sheet = wb.getSheetAt(sheet_number);

        for cell in sheet.getRow(0):
            column_index = cell.getColumnIndex()
            column_contents = cell.getStringCellValue()
            if column_contents == 'Project':
                project_column = column_index
            if column_contents == 'External ID' or column_contents == 'Individual ID':
                sample_column = column_index

        if project_column != None and sample_column != None:
            for row_number in range(1,sheet.getLastRowNum()+1):
                project = sheet.getRow(row_number).getCell(project_column).getStringCellValue()
                sample = sheet.getRow(row_number).getCell(sample_column).getStringCellValue()
                yield project,sample
            return

def tsv_reader(filename):
    f = open(filename,'rU')
    for line in f:
        tokens =line.split('\t')
        project = tokens[0].strip()
        sample = tokens[1].strip()
        yield project,sample
    f.close()    
        
def create_reader(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.xls' or extension == '.xlsx':
        return excel_reader(filename)
    elif extension == '.tsv' or extension == '.txt':
        return tsv_reader(filename)
    else:
        print 'Unrecognized file extension',extension
        sys.exit(1)

if len(sys.argv) != 2:
    print 'USAGE: %s <input file.{xls|xlsx|tsv|txt}>'
    sys.exit(1)
if not os.path.exists(sys.argv[1]):
    print 'Input file %s not found' % sys.argv[1]
    sys.exit(1)

input_filename = sys.argv[1]

for project,sample in create_reader(input_filename):
    sample_path = base_path % (project,sample)
    versions = []
    for version_path in os.listdir(sample_path):
        if version_path[0] != 'v':
            print 'Hit a path name that cannot be parsed: ',version_path
            sys.exit(1)
        versions.append(int(version_path[1:]))
    if len(versions) == 0:
        continue
    versions = sorted(versions)
    bam_file = '%s/v%d/%s.bam' % (sample_path,versions[-1],sample)
    if not os.path.exists(bam_file):
        print 'Malformed file: tried to find %s, but no such path exists' % bam_file
        sys.exit(1)
    print bam_file
