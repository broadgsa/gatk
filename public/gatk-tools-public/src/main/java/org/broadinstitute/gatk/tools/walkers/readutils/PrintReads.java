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
import org.broadinstitute.gatk.engine.io.NWaySAMFileWriter;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.iterators.ReadTransformersMode;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.File;
import java.util.*;

/**
 * Write out sequence read data (for filtering, merging, subsetting etc)
 *
 * <p>
 * PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically
 * merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can
 * also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf`
 * command line argument (see documentation on read filters for more information).
 * </p>
 *
 * <p>
 * Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow,
 * it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
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
 * <h3>Usage examples</h3>
 * <pre>
 // Prints all reads that have a mapping quality above zero
 * java -jar GenomeAnalysisTK.jar \
 *   -T PrintReads \
 *   -R reference.fasta \
 *   -I input1.bam \
 *   -I input2.bam \
 *   -o output.bam \
 *   --read_filter MappingQualityZero
 *
 * // Prints the first 2000 reads in the BAM file
 * java -jar GenomeAnalysisTK.jar \
 *   -T PrintReads \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -o output.bam \
 *   -n 2000
 *
 * // Downsamples BAM file to 25%
 * java -jar GenomeAnalysisTK.jar \
 *   -T PrintReads \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -o output.bam \
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
    GATKSAMFileWriter out;

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
    public Set<File> sampleFile = new TreeSet<>();

    /**
     * Only reads from the sample(s) will be included in the output.
     */
    @Argument(fullName="sample_name", shortName="sn", doc="Sample name to be included in the analysis. Can be specified multiple times.", required=false)
    public Set<String> sampleNames = new TreeSet<>();

    /**
     * Erase all extra attributes in the read but keep the read group information 
     */
    @Argument(fullName="simplify", shortName="s", doc="Simplify all reads", required=false)
    public boolean simplifyReads = false;

    @Hidden
    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="Don't output a program tag", required = false)
    public boolean NO_PG_TAG = false;

    List<ReadTransformer> readTransformers = Collections.emptyList();
    private Set<String> readGroupsToKeep = Collections.emptySet();

    public static final String PROGRAM_RECORD_NAME = "GATK PrintReads";   // The name that will go in the @PG tag
    
    Random random;


    /**
     * The initialize function.
     */
    public void initialize() {
        final GenomeAnalysisEngine toolkit = getToolkit();

        if ( toolkit != null )
            readTransformers = toolkit.getReadTransformers();

        //Sample names are case-insensitive
        final TreeSet<String> samplesToChoose = new TreeSet<>(new Comparator<String>() {
            @Override
            public int compare(String a, String b) {
                return a.compareToIgnoreCase(b);
            }
        });
        Collection<String> samplesFromFile;
        if (!sampleFile.isEmpty())  {
            samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFile);
            samplesToChoose.addAll(samplesFromFile);
        }

        if (!sampleNames.isEmpty())
            samplesToChoose.addAll(sampleNames);

        random = Utils.getRandomGenerator();

        if (toolkit != null) {
            final SAMFileHeader outputHeader = toolkit.getSAMFileHeader().clone();
            readGroupsToKeep = determineReadGroupsOfInterest(outputHeader, samplesToChoose);

            //If some read groups are to be excluded, remove them from the output header
            pruneReadGroups(outputHeader);

            //Add the program record (if appropriate) and set up the writer
            final boolean preSorted = true;
            if (toolkit.getArguments().BQSR_RECAL_FILE != null && !NO_PG_TAG ) {
                NWaySAMFileWriter.setupWriter(out, toolkit, outputHeader, preSorted, this, PROGRAM_RECORD_NAME);
            } else {
                out.writeHeader(outputHeader);
                out.setPresorted(preSorted);
            }

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
        // check that the read belongs to an RG that we need to keep
        if (!readGroupsToKeep.isEmpty()) {
            final SAMReadGroupRecord readGroup = read.getReadGroup();
            if (!readGroupsToKeep.contains(readGroup.getReadGroupId()))
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

    /**
     * Determines the list of read groups that meet the user's criteria for inclusion (based on id, platform, or sample)
     * @param header         the merged header for all input files
     * @param samplesToKeep  the list of specific samples specified by the user
     * @return               a Set of read group IDs that meet the user's criteria, empty if all RGs should be included
     */
    private Set<String> determineReadGroupsOfInterest(final SAMFileHeader header, final Set<String> samplesToKeep) {
        //If no filter options that use read group information have been supplied, exit early
        if (platform == null && readGroup == null && samplesToKeep.isEmpty())
            return Collections.emptySet();

        if  ( platform != null )
            platform = platform.toUpperCase();

        final Set<String> result = new HashSet<>();
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            // To be eligible for output, a read group must:
            //  NOT have an id that is blacklisted on the command line (note that String.equals(null) is false)
            //  AND NOT have a platform that contains the blacklisted platform from the command line
            //  AND have a sample that is whitelisted on the command line
            if (!rg.getReadGroupId().equals(readGroup) &&
                    (platform == null || !rg.getPlatform().toUpperCase().contains(platform)) &&
                    (samplesToKeep.isEmpty() || samplesToKeep.contains(rg.getSample())))
                result.add(rg.getReadGroupId());
        }

        if (result.isEmpty())
            throw new UserException.BadArgumentValue("-sn/-sf/-platform/-readGroup", "No read groups remain after pruning based on the supplied parameters");

        return result;
    }

    private void pruneReadGroups(final SAMFileHeader header) {
        if (readGroupsToKeep.isEmpty())
            return;

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            if (readGroupsToKeep.contains(rg.getReadGroupId()))
                readGroups.add(rg);
        }
        header.setReadGroups(readGroups);
    }
}
