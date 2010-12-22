#!/bin/tcsh

source /broad/tools/scripts/useuse

reuse Python-2.5
use R-2.11

setenv DIR /humgen/gsa-hpprojects/GATK/reports 
setenv ARCHIVE_DIR $DIR/archive
setenv SUMMARY_DIR $DIR/summaries
setenv DATE `date +"%m_%d_%Y"`
setenv ARCHIVE $ARCHIVE_DIR/$DATE
setenv SUMMARY $SUMMARY_DIR/$DATE
setenv GATK ~/dev/GenomeAnalysisTK/trunk

cd $DIR

echo "Archiving recently submitted jobs"
python $GATK/python/analyzeRunReports.py archive $DIR/submitted -o $ARCHIVE.gz -D

echo "All runs"
python $GATK/python/analyzeRunReports.py summary $ARCHIVE.gz --max_days 1

echo "No-dev"
python $GATK/python/analyzeRunReports.py summary $ARCHIVE.gz --max_days 1 --no-dev
python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE.gz --max_days 1 -E sting --no-dev

echo "Archive directory contents"
ls -ltrh $ARCHIVE_DIR

foreach maxDays ( 7 30 360 ) 
    echo "Creating table"
    setenv table $ARCHIVE.${maxDays}_days.table
    python $GATK/python/analyzeRunReports.py table $ARCHIVE_DIR/*.gz -o $table --max_days $maxDays

    echo "Creating summary"
    Rscript $GATK/R/GATKRunReport.R $table $SUMMARY.${maxDays}_days.pdf "of previous $maxDays days" 

    echo "Creating exception report"
    python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE_DIR/*.gz -o $SUMMARY.${maxDays}_days.sting.exceptions.txt --max_days $maxDays -E sting --no-dev
    python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE_DIR/*.gz -o $SUMMARY.${maxDays}_days.user.exceptions.txt --max_days $maxDays -E user --no-dev

    rm $table
end

echo "GATK daily run report" | mutt -a $SUMMARY.30_days.pdf -a $SUMMARY.360_days.pdf -a $SUMMARY.7_days.pdf -s "GATK Run report PDFs for $DATE" gsamembers
