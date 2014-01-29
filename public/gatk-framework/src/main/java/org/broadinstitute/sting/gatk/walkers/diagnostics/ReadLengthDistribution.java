/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Outputs the read lengths of all the reads in a file.
 *
 *  <p>
 *     Generates a table with the read lengths categorized per sample. If the file has no sample information
 *     (no read groups) it considers all reads to come from the same sample.
 *  </p>
 *
 *
 * <h3>Input</h3>
 *  <p>
 *      A BAM file.
 *  </p>
 *
 * <h3>Output</h3>
 *  <p>
 *      A human/R readable table of tab separated values with one column per sample and one row per read.
 *  </p>
 *
 * <h3>Examples</h3>
 *  <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ReadLengthDistribution
 *      -I example.bam
 *      -R reference.fasta
 *      -o example.tbl
 *  </pre>
 *
 * @author Kiran Garimela
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class ReadLengthDistribution extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    //A map from RG to its column number (its index in an int[] array)
    private Map<SAMReadGroupRecord,Integer> readGroupsLocation;
    //Each line in the table is a read length and each column it the number of reads of a specific RG with that length. Thus a table is a map between read lengths to array of values (one for each RG).
    private Map<Integer,int[]> table;
    private List<SAMReadGroupRecord> readGroups;

    public void initialize() {
        readGroups = getToolkit().getSAMFileHeader().getReadGroups();
        readGroupsLocation = new HashMap<>();
        table = new TreeMap<>();
        int readGroupsNum = 0;

        if (!readGroups.isEmpty()){
            for (SAMReadGroupRecord rg : readGroups){
                readGroupsLocation.put(rg,readGroupsNum);
                readGroupsNum++;
            }
        }
    }

    @Override
    public Integer map(final ReferenceContext referenceContext,final GATKSAMRecord samRecord,final RefMetaDataTracker RefMetaDataTracker) {

        final int length = Math.abs(samRecord.getReadLength());
        final SAMReadGroupRecord rg = samRecord.getReadGroup();

        increment(table,length, rg);

        return null;
    }

    final private void increment(final Map<Integer,int[]> table,final int length,final SAMReadGroupRecord rg){
        if(readGroupsLocation.isEmpty()){
            if(table.containsKey(length))
                table.get(length)[0]++;
            else{
                final int[] newLength = {1};
                table.put(length,newLength);
            }
        }
        else{
            final int rgLocation = readGroupsLocation.get(rg);
            if(table.containsKey(length))
                table.get(length)[rgLocation]++;
            else{
                table.put(length,new int[readGroupsLocation.size()]);
                table.get(length)[rgLocation]++;
            }
        }
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(final Integer integer,final Integer integer1) {
        return null;
    }

    public void onTraversalDone(final Integer sum) {
        final GATKReport report = createGATKReport();
        report.print(out);
    }

    final private GATKReport createGATKReport(){
        final GATKReport report = new GATKReport();
        report.addTable("ReadLengthDistribution", "Table of read length distributions", 1 + (readGroupsLocation.isEmpty() ? 1 : readGroupsLocation.size()));
        final GATKReportTable tableReport = report.getTable("ReadLengthDistribution");

        tableReport.addColumn("readLength");

        if (readGroupsLocation.isEmpty()){
            tableReport.addColumn("SINGLE_SAMPLE");
            int rowIndex = 0;
            for (Integer length : table.keySet()){
                tableReport.set(rowIndex,0,length);
                tableReport.set(rowIndex,1,table.get(length)[0]);
                rowIndex++;
            }
        }
        else{
            for (SAMReadGroupRecord rg : readGroups)
                tableReport.addColumn(rg.getSample());
            int rowIndex = 0;
            for (Integer length : table.keySet()){
                tableReport.set(rowIndex,0,length);
                for (int i=0; i < readGroupsLocation.size(); i++)
                    tableReport.set(rowIndex,i+1,table.get(length)[i]);
                rowIndex++;
            }

        }

        return report;
    }
}
