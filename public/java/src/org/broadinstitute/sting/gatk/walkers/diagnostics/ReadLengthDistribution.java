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
import java.util.List;

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

    private GATKReport report;

    public void initialize() {
        final List<SAMReadGroupRecord> readGroups = getToolkit().getSAMFileHeader().getReadGroups();

        report = new GATKReport();
        report.addTable("ReadLengthDistribution", "Table of read length distributions", 1 + (readGroups.isEmpty() ? 1 : readGroups.size()));
        GATKReportTable table = report.getTable("ReadLengthDistribution");

        table.addColumn("readLength");

        if (readGroups.isEmpty())
            table.addColumn("SINGLE_SAMPLE");
        else
            for (SAMReadGroupRecord rg : readGroups)
                table.addColumn(rg.getSample());
    }

    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return ( !read.getReadPairedFlag() || read.getReadPairedFlag() && read.getFirstOfPairFlag());
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord samRecord, RefMetaDataTracker RefMetaDataTracker) {
        GATKReportTable table = report.getTable("ReadLengthDistribution");

        int length = Math.abs(samRecord.getReadLength());
        String sample = samRecord.getReadGroup().getSample();

        table.increment(length, sample);

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer integer, Integer integer1) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        report.print(out);
    }
}
