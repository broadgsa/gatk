/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.readutils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.io.DirectOutputTracker;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.WalkerName;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Split a BAM file by sample
 *
 * <p>This tool divides the input data set into separate BAM files, one for each sample in the input data set. The split
 * files are named by concatenating the sample name to the end of the provided outputRoot command-line argument.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A single bam file.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A separate bam file for each sample.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SplitSamFile \
 *   -R reference.fasta \
 *   -I input.bam \
 *   --outputRoot myproject_
 * </pre>
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@WalkerName("SplitSamFile")
@Requires({DataSource.READS})
public class SplitSamFile extends ReadWalker<SAMRecord, Map<String, SAMFileWriter>> {
    @Argument(fullName="outputRoot", doc="output BAM file", required=false)
    public String outputRoot = "";

    private static final Logger logger = Logger.getLogger(SplitSamFile.class);
    private static final String VERSION = "0.0.1";

    @Override
    public void initialize() {
        logger.info("SplitSamFile version: " + VERSION);
    }

    @Override
    public SAMRecord map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        return read;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Standard I/O routines
    //
    // --------------------------------------------------------------------------------------------------------------
    @Override
    public void onTraversalDone(Map<String, SAMFileWriter> outputs) {
        for ( SAMFileWriter output : outputs.values() ) {
            output.close();
        }
    }

    @Override
    public Map<String, SAMFileWriter> reduceInit() {
        HashMap<String, SAMFileHeader> headers = new HashMap<>();
        for ( SAMReadGroupRecord readGroup : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            final String sample = readGroup.getSample();
            if ( ! headers.containsKey(sample) ) {
                SAMFileHeader header = duplicateSAMFileHeader(this.getToolkit().getSAMFileHeader());
                logger.debug(String.format("Creating BAM header for sample %s", sample));
                ArrayList<SAMReadGroupRecord> readGroups = new ArrayList<>();
                header.setReadGroups(readGroups);
                headers.put(sample, header);
            }

            SAMFileHeader header = headers.get(sample);
            List<SAMReadGroupRecord> newReadGroups = new ArrayList<>(header.getReadGroups());
            newReadGroups.add(readGroup);
            header.setReadGroups(newReadGroups);
        }

        HashMap<String, SAMFileWriter> outputs = new HashMap<>();
        final OutputTracker outputTracker = new DirectOutputTracker();
        for ( Map.Entry<String, SAMFileHeader> elt : headers.entrySet() ) {
            final String sample = elt.getKey();
            final String filename = outputRoot + sample + ".bam";
            logger.info(String.format("Creating BAM output file %s for sample %s", filename, sample));

            final SAMFileWriter output = SAMFileWriterStub.createSAMFileWriter(filename, getToolkit(), elt.getValue());
            outputs.put(sample, output);
            outputTracker.addOutput( (SAMFileWriterStub) output);
        }

        return outputs;
    }

    /**
     * Write out the read
     */
    @Override
    public Map<String, SAMFileWriter> reduce(SAMRecord read, Map<String, SAMFileWriter> outputs) {
        final String sample = read.getReadGroup().getSample();
        SAMFileWriter output = outputs.get(sample);

        if ( output != null ) {
            output.addAlignment(read);
        } else {
            throw new RuntimeException(String.format("Read group %s not present in header but found in read %s", read.getReadGroup().getReadGroupId(), read.getReadName()));
        }

        return outputs;
    }

    public static SAMFileHeader duplicateSAMFileHeader(SAMFileHeader toCopy) {
        SAMFileHeader copy = new SAMFileHeader();

        copy.setSortOrder(toCopy.getSortOrder());
        copy.setGroupOrder(toCopy.getGroupOrder());
        copy.setProgramRecords(toCopy.getProgramRecords());
        copy.setReadGroups(toCopy.getReadGroups());
        copy.setSequenceDictionary(toCopy.getSequenceDictionary());

        for (Map.Entry<String, String> e : toCopy.getAttributes())
            copy.setAttribute(e.getKey(), e.getValue());

        return copy;
    }

}