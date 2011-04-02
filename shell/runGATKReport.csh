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
setenv GATK_RELEASE_VERSION `ls -l /humgen/gsa-hpprojects/GATK/bin/current | sed 's/.*GenomeAnalysisTK-//'`
setenv REPORT_TXT $DIR/report.txt

rm -f $REPORT_TXT 

cd $DIR

echo "\n####################\nArchiving recently submitted jobs" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py archive $DIR/submitted -o $ARCHIVE.gz -D >> $REPORT_TXT

echo "\n####################\nReleased version ($GATK_RELEASE_VERSION), all runs" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py summary $ARCHIVE_DIR/*.gz --rev $GATK_RELEASE_VERSION >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE_DIR/*.gz -E sting --rev $GATK_RELEASE_VERSION >> $REPORT_TXT

echo "\n####################\nLast day, all versions" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py summary $ARCHIVE.gz --max_days 1 --no-dev >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE.gz --max_days 1 -E sting --no-dev >> $REPORT_TXT

#echo "Archive directory contents"
#du -sh $ARCHIVE_DIR

if (1 == 0) then
foreach maxDays ( 30 360 ) 
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
endif

#echo "GATK daily run report" | mutt -a $SUMMARY.30_days.pdf -a $SUMMARY.360_days.pdf -a $SUMMARY.7_days.pdf -s "GATK Run report PDFs for $DATE" gsamembers
#cat $REPORT_TXT | mutt -a $REPORT_TXT -a $SUMMARY.30_days.pdf -a $SUMMARY.360_days.pdf -s "GATK run report for $DATE" gsamembers
cat $REPORT_TXT | mutt -a $REPORT_TXT -s "GATK run report for $DATE" gsamembers

