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

package org.broadinstitute.sting.gatk.walkers.readutils;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.iterators.ReadTransformersMode;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;
import java.util.*;

/**
 * Renders, in SAM/BAM format, all reads from the input data set in the order in which they appear in the input file.
 *
 * <p>
 * PrintReads can dynamically merge the contents of multiple input BAM files, resulting
 * in merged output sorted in coordinate order.  Can also optionally filter reads based on the
 * --read_filter command line argument.
 * </p>
 *
 * <p>
 * Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow,
 * it takes the --BQSR engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more bam files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A single processed bam file.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T PrintReads \
 *   -o output.bam \
 *   -I input1.bam \
 *   -I input2.bam \
 *   --read_filter MappingQualityZero
 *
 * // Prints the first 2000 reads in the BAM file
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T PrintReads \
 *   -o output.bam \
 *   -I input.bam \
 *   -n 2000
 *
 * // Downsamples BAM file to 25%
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T PrintReads \
 *   -o output.bam \
 *   -I input.bam \
 *   -dfrac 0.25
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@ReadTransformersMode(ApplicationTime = ReadTransformer.ApplicationTime.HANDLED_IN_WALKER)
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = ReadTransformer.ApplicationTime.HANDLED_IN_WALKER)
@Requires({DataSource.READS, DataSource.REFERENCE})
public class PrintReads extends ReadWalker<GATKSAMRecord, SAMFileWriter> implements NanoSchedulable {

    @Output(doc="Write output to this BAM filename instead of STDOUT")
    StingSAMFileWriter out;

    @Argument(fullName = "readGroup", shortName = "readGroup", doc="Exclude all reads with this read group from the output", required = false)
    String readGroup = null;

    /**
     * For example, --platform ILLUMINA or --platform 454.
     */
    @Argument(fullName = "platform", shortName = "platform", doc="Exclude all reads with this platform from the output", required = false)
    String platform = null;

    /**
     * Only prints the first n reads of the file
     */
    @Argument(fullName = "number", shortName = "n", doc="Print the first n reads from the file, discarding the rest", required = false)
    int nReadsToPrint = -1;

    /**
     * Only reads from samples listed in the provided file(s) will be included in the output.
     */
    @Argument(fullName="sample_file", shortName="sf", doc="File containing a list of samples (one per line). Can be specified multiple times", required=false)
    public Set<File> sampleFile = new TreeSet<File>();

    /**
     * Only reads from the sample(s) will be included in the output.
     */
    @Argument(fullName="sample_name", shortName="sn", doc="Sample name to be included in the analysis. Can be specified multiple times.", required=false)
    public Set<String> sampleNames = new TreeSet<String>();

    /**
     * Erase all extra attributes in the read but keep the read group information 
     */
    @Argument(fullName="simplify", shortName="s", doc="Simplify all reads.", required=false)
    public boolean simplifyReads = false;

    @Hidden
    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="", required = false)
    public boolean NO_PG_TAG = false;

    List<ReadTransformer> readTransformers = Collections.emptyList();
    private TreeSet<String> samplesToChoose = new TreeSet<String>();
    private boolean SAMPLES_SPECIFIED = false;

    public static final String PROGRAM_RECORD_NAME = "GATK PrintReads";   // The name that will go in the @PG tag
    
    Random random;


    /**
     * The initialize function.
     */
    public void initialize() {
        final GenomeAnalysisEngine toolkit = getToolkit();

        if  ( platform != null )
            platform = platform.toUpperCase();

        if ( getToolkit() != null )
            readTransformers = getToolkit().getReadTransformers();

        Collection<String> samplesFromFile;
        if (!sampleFile.isEmpty())  {
            samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFile);
            samplesToChoose.addAll(samplesFromFile);
        }

        if (!sampleNames.isEmpty())
            samplesToChoose.addAll(sampleNames);

        if(!samplesToChoose.isEmpty()) {
            SAMPLES_SPECIFIED = true;
        }

        random = GenomeAnalysisEngine.getRandomGenerator();

        final boolean preSorted = true;
        if (getToolkit() != null && getToolkit().getArguments().BQSR_RECAL_FILE != null && !NO_PG_TAG ) {
            Utils.setupWriter(out, toolkit, toolkit.getSAMFileHeader(), preSorted, this, PROGRAM_RECORD_NAME);
        }

    }

    /**
     * The reads filter function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a GATKSAMRecord
     * @return true if the read passes the filter, false if it doesn't
     */
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        // check the read group
        if  ( readGroup != null ) {
            SAMReadGroupRecord myReadGroup = read.getReadGroup();
            if ( myReadGroup == null || !readGroup.equals(myReadGroup.getReadGroupId()) )
                return false;
        }

        // check the platform
        if  ( platform != null ) {
            SAMReadGroupRecord readGroup = read.getReadGroup();
            if ( readGroup == null )
                return false;

            Object readPlatformAttr = readGroup.getAttribute("PL");
            if ( readPlatformAttr == null || !readPlatformAttr.toString().toUpperCase().contains(platform))
                return false;
        }
        if (SAMPLES_SPECIFIED )  {
            // user specified samples to select
            // todo - should be case-agnostic  but for simplicity and speed this is ignored.
            // todo - can check at initialization intersection of requested samples and samples in BAM header to further speedup.
            if (!samplesToChoose.contains(read.getReadGroup().getSample()))
                return false;
        }


        // check if we've reached the output limit
        if ( nReadsToPrint == 0 ) {
            return false;          // n == 0 means we've printed all we needed.
        }
        else if (nReadsToPrint > 0) {
            nReadsToPrint--;       // n > 0 means there are still reads to be printed.
        }

        return true;
    }

    /**
     * The reads map function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param readIn the read itself, as a GATKSAMRecord
     * @return the read itself
     */
    public GATKSAMRecord map( ReferenceContext ref, GATKSAMRecord readIn, RefMetaDataTracker metaDataTracker ) {
        GATKSAMRecord workingRead = readIn;

        for ( final ReadTransformer transformer : readTransformers ) {
            workingRead = transformer.apply(workingRead);
        }

        if ( simplifyReads ) workingRead = workingRead.simplify();

        return workingRead;
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    public SAMFileWriter reduceInit() {
        return out;
    }

    /**
     * given a read and a output location, reduce by emitting the read
     *
     * @param read   the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public SAMFileWriter reduce( GATKSAMRecord read, SAMFileWriter output ) {
        output.addAlignment(read);
        return output;
    }
}
