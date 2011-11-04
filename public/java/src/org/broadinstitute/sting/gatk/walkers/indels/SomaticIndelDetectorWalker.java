/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.*;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.PlatformUnitFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqCodec;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.refseq.Transcript;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.collections.CircularArray;
import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.broadinstitute.sting.utils.interval.OverlappingIntervalIterator;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.*;


/**
 * Tool for calling indels in Tumor-Normal paired sample mode; this tool supports single-sample mode as well,
 * but this latter functionality is now superceded by UnifiedGenotyper.
 *
 * <p>
 * This is a simple, counts-and-cutoffs based tool for calling indels from aligned (preferrably MSA cleaned) sequencing
 * data. Supported output formats are: BED format, extended verbose output (tab separated), and VCF. The latter two outputs
 * include additional statistics such as mismtaches and base qualitites around the calls, read strandness (how many
 * forward/reverse reads support ref and indel alleles) etc. It is highly recommended to use these additional
 * statistics to perform post-filtering of the calls as the tool is tuned for sensitivity (in other words it will
 * attempt to "call" anything remotely reasonable based only on read counts and will generate all the additional
 * metrics for the post-processing tools to make the final decision). The calls are performed by default
 * from a matched tumor-normal pair of samples. In this case, two (sets of) input bam files must be specified using tagged -I
 * command line arguments: normal and tumor bam(s) must be passed with -I:normal and -I:tumor arguments,
 * respectively. Indels are called from the tumor sample and annotated as germline
 * if even a weak evidence for the same indel, not necessarily a confident call, exists in the normal sample, or as somatic
 * if normal sample has coverage at the site but no indication for an indel. Note that strictly speaking the calling
 * is not even attempted in normal sample: if there is an indel in normal that is not detected/does not pass a threshold
 * in tumor sample, it will not be reported.
 *
 * To make indel calls and associated metrics for a single sample, this tool can be run with --unpaired flag (input
 * bam tagging is not required in this case, and tags are completely ignored if still used: all input bams will be merged
 * on the fly and assumed to represent a single sample - this tool does not check for sample id in the read groups).
 *
 * <h2>Input</h2>
 * <p>
 * Tumor and normal bam files (or single sample bam file(s) in --unpaired mode).
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Indel calls with associated metrics.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SomaticIndelDetector \
 *   -o indels.vcf \
 *   -verbose indels.txt
 *   -I:normal normal.bam \
 *   -I:tumor tumor.bam
 * </pre>
 *
 */

@ReadFilters({Platform454Filter.class, MappingQualityZeroFilter.class, PlatformUnitFilter.class})
public class SomaticIndelDetectorWalker extends ReadWalker<Integer,Integer> {
//    @Output
//    PrintStream out;
    @Output(doc="File to write variants (indels) in VCF format",required=true)
    protected VCFWriter vcf_writer = null;

    @Argument(fullName="outputFile", shortName="O", doc="output file name (BED format). DEPRECATED> Use --bed", required=true)
    @Deprecated
    java.io.File output_file;

    @Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print callability metrics output", required = false)
    public PrintStream metricsWriter = null;

//    @Argument(fullName="vcf_format", shortName="vcf", doc="generate output file in VCF format", required=false)
//    boolean FORMAT_VCF = false;

    @Hidden
    @Input(fullName = "genotype_intervals", shortName = "genotype",
        doc = "Calls will be made at each position within the specified interval(s), whether there is an indel or not", required = false)
    public IntervalBinding<Feature> genotypeIntervalsFile = null;

    @Hidden
    @Argument(fullName="unpaired", shortName="unpaired",
                    doc="Perform unpaired calls (no somatic status detection)", required=false)
    boolean call_unpaired = false;
    boolean call_somatic ;

    @Argument(fullName="verboseOutput", shortName="verbose",
                    doc="Verbose output file in text format", required=false)
    java.io.File verboseOutput = null;

    @Argument(fullName="bedOutput", shortName="bed",
        doc="Lightweight bed output file (only positions and events, no stats/annotations)", required=false)
    java.io.File bedOutput = null;

    @Argument(fullName="minCoverage", shortName="minCoverage",
                    doc="indel calls will be made only at sites with tumor coverage of minCoverage or more reads; "+
            "with --unpaired (single sample) option, this value is used for minimum sample coverage", required=false)
    int minCoverage = 6;

    @Argument(fullName="minNormalCoverage", shortName="minNormalCoverage",
                    doc="used only in default (somatic) mode;  normal sample must have at least minNormalCoverage "+
            "or more reads at the site to call germline/somatic indel, otherwise the indel (in tumor) is ignored", required=false)
    int minNormalCoverage = 4;

    @Argument(fullName="minFraction", shortName="minFraction",
                    doc="Minimum fraction of reads with CONSENSUS indel at a site, out of all reads covering the site, required for making a call"+
                    " (fraction of non-consensus indels at the site is not considered here, see minConsensusFraction)", required=false)
    double minFraction = 0.3;

    @Argument(fullName="minConsensusFraction", shortName="minConsensusFraction",
                    doc="Indel call is made only if fraction of CONSENSUS indel observations at a site wrt "+
            "all indel observations at the site exceeds this threshold", required=false)
    double minConsensusFraction = 0.7;

    @Argument(fullName="minIndelCount", shortName="minCnt",
                    doc="Minimum count of reads supporting consensus indel required for making the call. "+
                    " This filter supercedes minFraction, i.e. indels with acceptable minFraction at low coverage "+
                    "(minIndelCount not met) will not pass.", required=false)
    int minIndelCount = 0;

    @Argument(fullName="refseq", shortName="refseq",
                    doc="Name of RefSeq transcript annotation file. If specified, indels will be annotated with "+
            "GENOMIC/UTR/INTRON/CODING and with the gene name", required=false)
    String RefseqFileName = null;

//@Argument(fullName="blacklistedLanes", shortName="BL",
//        doc="Name of lanes (platform units) that should be ignored. Reads coming from these lanes will never be seen "+
//                "by this application, so they will not contribute indels to consider and will not be counted.", required=false)
//PlatformUnitFilterHelper dummy;

    @Hidden
    @Argument(fullName="indel_debug", shortName="idebug", doc="Detailed printout for debugging, do not turn this on",
            required=false) Boolean DEBUG = false;
    @Argument(fullName="window_size", shortName="ws", doc="Size (bp) of the sliding window used for accumulating the coverage. "+
            "May need to be increased to accomodate longer reads or longer deletions. A read can be fit into the "+
            "window if its length on the reference (i.e. read length + length of deletion gap(s) if any) is smaller "+
            "than the window size. Reads that do not fit will be ignored, so long deletions can not be called "+
            "if window is too small",required=false) int WINDOW_SIZE = 200;
    @Argument(fullName="maxNumberOfReads",shortName="mnr",doc="Maximum number of reads to cache in the window; if number of reads exceeds this number,"+
                " the window will be skipped and no calls will be made from it",required=false) int MAX_READ_NUMBER = 10000;



	private WindowContext tumor_context;
	private WindowContext normal_context; 
	private int currentContigIndex = -1;
    private int contigLength = -1; // we see to much messy data with reads hanging out of contig ends...
	private int currentPosition = -1; // position of the last read we've seen on the current contig
	private String refName = null;
	private java.io.Writer output = null;
	private GenomeLoc location = null;
    private long normalCallsMade = 0L, tumorCallsMade = 0L;

    boolean outOfContigUserWarned = false;

    private LocationAwareSeekableRODIterator refseqIterator=null;

//	private Set<String> normalReadGroups; // we are going to remember which read groups are normals and which are tumors in order to be able
//	private Set<String> tumorReadGroups ; // to properly assign the reads coming from a merged stream
    private Set<String> normalSamples; // we are going to remember which samples are normal and which are tumor:
    private Set<String> tumorSamples ; // these are used only to generate genotypes for vcf output

	private int NQS_WIDTH = 5; // 5 bases on each side of the indel for NQS-style statistics

    private Writer bedWriter = null;
    private Writer verboseWriter = null;


	private static String annGenomic = "GENOMIC";
	private static String annIntron = "INTRON";
	private static String annUTR = "UTR";
	private static String annCoding = "CODING";
	private static String annUnknown = "UNKNOWN";

    enum CallType {
        NOCOVERAGE,
        BADCOVERAGE,
        NOEVIDENCE,
        GERMLINE,
        SOMATIC
    };

	private SAMRecord lastRead;
    private byte[] refBases;
    private ReferenceDataSource refData;
    private Iterator<GenomeLoc> genotypeIntervalIterator = null;

    // the current interval in the list of intervals, for which we want to do full genotyping
    private GenomeLoc currentGenotypeInterval = null;
    private long lastGenotypedPosition = -1; // last position on the currentGenotypeInterval, for which a call was already printed;
                                     // can be 1 base before lastGenotyped start


    // "/humgen/gsa-scr1/GATK_Data/refGene.sorted.txt"

    private Set<VCFHeaderLine> getVCFHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // first, the basic info
        headerInfo.add(new VCFHeaderLine("source", "SomaticIndelDetector"));
        headerInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        // FORMAT and INFO fields
//        headerInfo.addAll(VCFUtils.getSupportedHeaderStrings());

        headerInfo.addAll(VCFIndelAttributes.getAttributeHeaderLines());
        if ( call_somatic ) {
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.SOMATIC_KEY, 0, VCFHeaderLineType.Flag, "Somatic event"));
        }  else {
        }

        // all of the arguments from the argument collection
        Set<Object> args = new HashSet<Object>();
        args.add(this);
        args.addAll(getToolkit().getFilters());
        Map<String,String> commandLineArgs = getToolkit().getApproximateCommandLineArguments(args);
        for ( Map.Entry<String, String> commandLineArg : commandLineArgs.entrySet() )
            headerInfo.add(new VCFHeaderLine(String.format("SID_%s", commandLineArg.getKey()), commandLineArg.getValue()));
        // also, the list of input bams
        for ( String fileName : getToolkit().getArguments().samFiles )
            headerInfo.add(new VCFHeaderLine("SID_bam_file_used", fileName));

        return headerInfo;
    }


	@Override
	public void initialize() {

        call_somatic =  (call_unpaired ? false : true);
		normal_context = new WindowContext(0,WINDOW_SIZE);
        normalSamples = new HashSet<String>();

        if ( bedOutput != null && output_file != null ) {
            throw new UserException.DeprecatedArgument("-O", "-O option is deprecated and -bed option replaces it; you can not use both at the same time");
        }

		if ( RefseqFileName != null ) {
            logger.info("Using RefSeq annotations from "+RefseqFileName);

			RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                                                          getToolkit().getGenomeLocParser(),
                                                          getToolkit().getArguments().unsafe);
            RMDTrack refseq = builder.createInstanceOfTrack(RefSeqCodec.class,new File(RefseqFileName));

            refseqIterator = new SeekableRODIterator(refseq.getHeader(),
                                                     refseq.getSequenceDictionary(),
                                                     getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                                                     getToolkit().getGenomeLocParser(),
                                                     refseq.getIterator());
		}

		if ( refseqIterator == null ) logger.info("No gene annotations available");

		int nSams = getToolkit().getArguments().samFiles.size();

        if ( call_somatic ) {
            if ( nSams < 2 ) throw new UserException.BadInput("In default (paired sample) mode at least two bam files (normal and tumor) must be specified");
            tumor_context = new WindowContext(0,WINDOW_SIZE);
            tumorSamples = new HashSet<String>();
        }

        int nNorm = 0;
        int nTum = 0;
        for ( SAMReaderID rid : getToolkit().getReadsDataSource().getReaderIDs() ) {
             Tags tags = rid.getTags() ;
             if ( tags.getPositionalTags().isEmpty() && call_somatic )
                 throw new UserException.BadInput("In default (paired sample) mode all input bam files must be tagged as either 'normal' or 'tumor'. Untagged file: "+
                         getToolkit().getSourceFileForReaderID(rid));
             boolean normal = false;
             boolean tumor = false;
             for ( String s : tags.getPositionalTags() ) { // we allow additional unrelated tags (and we do not use them), but we REQUIRE one of Tumor/Normal to be present if --somatic is on
                 if ( "NORMAL".equals(s.toUpperCase()) ) {
                     normal = true;
                     nNorm++;
                 }
                 if ( "TUMOR".equals(s.toUpperCase()) ) {
                     tumor = true;
                     nTum++ ;
                 }
             }
             if ( call_somatic && normal && tumor ) throw new UserException.BadInput("Input bam file "+
                     getToolkit().getSourceFileForReaderID(rid)+" is tagged both as normal and as tumor. Which one is it??");
             if ( call_somatic && !normal && ! tumor )
                 throw new UserException.BadInput("In somatic mode all input bams must be tagged as either normal or tumor. Encountered untagged file: "+
                    getToolkit().getSourceFileForReaderID(rid));
             if ( ! call_somatic && (normal || tumor) )
                 System.out.println("WARNING: input bam file "+getToolkit().getSourceFileForReaderID(rid)
                         +" is tagged as Normal and/or Tumor, but somatic mode is not on. Tags will ne IGNORED");
            if ( call_somatic && tumor ) {
                for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader(rid).getReadGroups() ) {
                    tumorSamples.add(rg.getSample());
                }
            } else {
                for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader(rid).getReadGroups() ) {
                    normalSamples.add(rg.getSample());
                }
            }
            if ( genotypeIntervalsFile != null ) {

                // read in the whole list of intervals for cleaning
                GenomeLocSortedSet locs = IntervalUtils.sortAndMergeIntervals(getToolkit().getGenomeLocParser(), genotypeIntervalsFile.getIntervals(getToolkit()), IntervalMergingRule.OVERLAPPING_ONLY);
                genotypeIntervalIterator = locs.iterator();

                // wrap intervals requested for genotyping inside overlapping iterator, so that we actually
                // genotype only on the intersections of the requested intervals with the -L intervals
                genotypeIntervalIterator = new OverlappingIntervalIterator(genotypeIntervalIterator, getToolkit().getIntervals().iterator() );

                currentGenotypeInterval = genotypeIntervalIterator.hasNext() ? genotypeIntervalIterator.next() : null;

                if ( DEBUG) System.out.println("DEBUG>> first genotyping interval="+currentGenotypeInterval);

                if ( currentGenotypeInterval != null ) lastGenotypedPosition = currentGenotypeInterval.getStart()-1;
            }

        }

		location = getToolkit().getGenomeLocParser().createGenomeLoc(getToolkit().getSAMFileHeader().getSequence(0).getSequenceName(),1);

        normalSamples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeaders().get(0));

        try {
            // we already checked that bedOutput and output_file are not set simultaneously
            if ( bedOutput != null ) bedWriter = new FileWriter(bedOutput);
            if ( output_file != null ) bedWriter = new FileWriter(output_file);
        } catch (java.io.IOException e) {
            throw new UserException.CouldNotReadInputFile(bedOutput, "Failed to open BED file for writing.", e);
        }
        try {
            if ( verboseOutput != null ) verboseWriter = new FileWriter(verboseOutput);
        } catch (java.io.IOException e) {
            throw new UserException.CouldNotReadInputFile(verboseOutput, "Failed to open BED file for writing.", e);
        }

        vcf_writer.writeHeader(new VCFHeader(getVCFHeaderInfo(), SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader()))) ;
        refData = new ReferenceDataSource(getToolkit().getArguments().referenceFile);
	}


	@Override
	public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {

    //        if ( read.getReadName().equals("428EFAAXX090610:2:36:1384:639#0") ) System.out.println("GOT READ");

            if ( DEBUG ) {
                //            System.out.println("DEBUG>> read at "+ read.getAlignmentStart()+"-"+read.getAlignmentEnd()+
                //                    "("+read.getCigarString()+")");
                if ( read.getDuplicateReadFlag() ) System.out.println("DEBUG>> Duplicated read (IGNORED)");
            }

            if ( AlignmentUtils.isReadUnmapped(read) ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) {
			return 0; // we do not need those reads!
            }

            if ( read.getReferenceIndex() != currentContigIndex ) {
                // we just jumped onto a new contig
                if ( DEBUG ) System.out.println("DEBUG>>> Moved to contig "+read.getReferenceName());
                if ( read.getReferenceIndex() < currentContigIndex ) // paranoidal
                    throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, read, "Read "+read.getReadName()+": contig is out of order; input BAM file is unsorted");

                // print remaining indels from the previous contig (if any);
                if ( call_somatic ) emit_somatic(1000000000, true);
                else emit(1000000000,true);

                currentContigIndex = read.getReferenceIndex();
                currentPosition = read.getAlignmentStart();
                refName = new String(read.getReferenceName());

                location = getToolkit().getGenomeLocParser().createGenomeLoc(refName,location.getStart(),location.getStop());
                contigLength = getToolkit().getGenomeLocParser().getContigInfo(refName).getSequenceLength();
                outOfContigUserWarned = false;

                lastGenotypedPosition = -1;

                normal_context.clear(); // reset coverage window; this will also set reference position to 0
                if ( call_somatic) tumor_context.clear();

                refBases = new String(refData.getReference().getSequence(read.getReferenceName()).getBases()).toUpperCase().getBytes();
            }

            // we have reset the window to the new contig if it was required and emitted everything we collected
            // on a previous contig. At this point we are guaranteed that we are set up properly for working
            // with the contig of the current read.

            // NOTE: all the sanity checks and error messages below use normal_context only. We make sure that normal_context and
            // tumor_context are synchronized exactly (windows are always shifted together by emit_somatic), so it's safe

            if ( read.getAlignmentStart() < currentPosition ) // oops, read out of order?
                throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, read, "Read "+read.getReadName() +" out of order on the contig\n"+
                        "Read starts at "+refName+":"+read.getAlignmentStart()+"; last read seen started at "+refName+":"+currentPosition
                        +"\nLast read was: "+lastRead.getReadName()+" RG="+lastRead.getAttribute("RG")+" at "+lastRead.getAlignmentStart()+"-"
                        +lastRead.getAlignmentEnd()+" cigar="+lastRead.getCigarString());

            currentPosition = read.getAlignmentStart();
            lastRead = read;

            if ( read.getAlignmentEnd() > contigLength  ) {
                if ( ! outOfContigUserWarned ) {
                    System.out.println("WARNING: Reads aligned past contig length on "+ location.getContig()+"; all such reads will be skipped");
                    outOfContigUserWarned = true;
                }
                return 0;
            }

            long alignmentEnd = read.getAlignmentEnd();
            Cigar c = read.getCigar();
            int lastNonClippedElement = 0; // reverse offset to the last unclipped element
            CigarOperator op = null;
            // moving backwards from the end of the cigar, skip trailing S or H cigar elements:
            do {
                lastNonClippedElement++;
                op = c.getCigarElement( c.numCigarElements()-lastNonClippedElement ).getOperator();
            } while ( op == CigarOperator.H || op == CigarOperator.S );

            // now op is the last non-S/H operator in the cigar.

            // a little trick here: we want to make sure that current read completely fits into the current
            // window so that we can accumulate indel observations over the whole length of the read.
            // The ::getAlignmentEnd() method returns the last position on the reference where bases from the
            // read actually match (M cigar elements). After our cleaning procedure, we can have reads that end
            // with I element, which is not gonna be counted into alignment length on the reference. On the other hand,
            // in this program we assign insertions, internally, to the first base *after* the insertion position.
            // Hence, we have to make sure that that extra base is already in the window or we will get IndexOutOfBounds.

            if ( op == CigarOperator.I) alignmentEnd++;

            if ( alignmentEnd > normal_context.getStop()) {

                // we don't emit anything until we reach a read that does not fit into the current window.
                // At that point we try shifting the window to the start of that read (or reasonably close) and emit everything prior to
                // that position. This is legitimate, since the reads are sorted and  we are not gonna see any more coverage at positions
                // below the current read's start.
                // Clearly, we assume here that window is large enough to accomodate any single read, so simply shifting
                // the window to around the read's start will ensure that the read fits...

                if ( DEBUG) System.out.println("DEBUG>> Window at "+normal_context.getStart()+"-"+normal_context.getStop()+", read at "+
                                read.getAlignmentStart()+": trying to emit and shift" );
                if ( call_somatic ) emit_somatic( read.getAlignmentStart(), false );
                else emit( read.getAlignmentStart(), false );

                // let's double check now that the read fits after the shift
                if ( read.getAlignmentEnd() > normal_context.getStop()) {
                    // ooops, looks like the read does not fit into the window even after the latter was shifted!!
                    // we used to die over such reads and require user to run with larger window size. Now we
                    // just print a warning and discard the read (this means that our counts can be slightly off in
                    // th epresence of such reads)
                    //throw new UserException.BadArgumentValue("window_size", "Read "+read.getReadName()+": out of coverage window bounds. Probably window is too small, so increase the value of the window_size argument.\n"+
                    //                         "Read length="+read.getReadLength()+"; cigar="+read.getCigarString()+"; start="+
                    //                         read.getAlignmentStart()+"; end="+read.getAlignmentEnd()+
                    //                         "; window start (after trying to accomodate the read)="+normal_context.getStart()+"; window end="+normal_context.getStop());
                    System.out.println("WARNING: Read "+read.getReadName()+
                             " is out of coverage window bounds. Probably window is too small and the window_size value must be increased.\n"+
                             "  The read is ignored in this run (so all the counts/statistics reported will not include it).\n"+
                                             "  Read length="+read.getReadLength()+"; cigar="+read.getCigarString()+"; start="+
                                             read.getAlignmentStart()+"; end="+read.getAlignmentEnd()+
                                             "; window start (after trying to accomodate the read)="+normal_context.getStart()+"; window end="+normal_context.getStop());
                    return 1;
                }
            }

            if ( call_somatic ) {

                Tags tags =  getToolkit().getReaderIDForRead(read).getTags();
                boolean assigned = false;
                for ( String s : tags.getPositionalTags() ) {
                    if ( "NORMAL".equals(s.toUpperCase()) ) {
                        normal_context.add(read,ref.getBases());
                        assigned = true;
                        break;
                    }
                    if ( "TUMOR".equals(s.toUpperCase()) ) {
                        tumor_context.add(read,ref.getBases());
                        assigned = true;
                        break;
                    }
                }
                if ( ! assigned )
                    throw new StingException("Read "+read.getReadName()+" from "+getToolkit().getSourceFileForReaderID(getToolkit().getReaderIDForRead(read))+
                    "has no Normal/Tumor tag associated with it");

//                String rg = (String)read.getAttribute("RG");
//                if ( rg == null )
//                    throw new UserException.MalformedBam(read, "Read "+read.getReadName()+" has no read group in merged stream. RG is required for somatic calls.");

//                if ( normalReadGroups.contains(rg) ) {
//                    normal_context.add(read,ref.getBases());
//                } else if ( tumorReadGroups.contains(rg) ) {
//                    tumor_context.add(read,ref.getBases());
//                } else {
//                    throw new UserException.MalformedBam(read, "Unrecognized read group in merged stream: "+rg);
//                }

                if ( tumor_context.getReads().size() > MAX_READ_NUMBER ) {
                    System.out.println("WARNING: a count of "+MAX_READ_NUMBER+" reads reached in a window "+
                            refName+':'+tumor_context.getStart()+'-'+tumor_context.getStop()+" in tumor sample. The whole window will be dropped.");
                    tumor_context.shift(WINDOW_SIZE);
                    normal_context.shift(WINDOW_SIZE);
                }
                if ( normal_context.getReads().size() > MAX_READ_NUMBER ) {
                    System.out.println("WARNING: a count of "+MAX_READ_NUMBER+" reads reached in a window "+
                            refName+':'+normal_context.getStart()+'-'+normal_context.getStop()+" in normal sample. The whole window will be dropped");
                    tumor_context.shift(WINDOW_SIZE);
                    normal_context.shift(WINDOW_SIZE);
                }


            } else {
                normal_context.add(read, ref.getBases());
                if ( normal_context.getReads().size() > MAX_READ_NUMBER ) {
                    System.out.println("WARNING: a count of "+MAX_READ_NUMBER+" reads reached in a window "+
                            refName+':'+normal_context.getStart()+'-'+normal_context.getStop()+". The whole window will be dropped");
                    normal_context.shift(WINDOW_SIZE);
                }
            }

            return 1;
	}

    /** An auxiliary shortcut: returns true if position(location.getContig(), p) is past l  */
    private boolean pastInterval(long p, GenomeLoc l) {
        return ( location.getContigIndex() > l.getContigIndex() ||
                 location.getContigIndex() == l.getContigIndex() && p > l.getStop() );
    }

    /** Emit calls of the specified type across genotyping intervals, from position lastGenotypedPosition+1 to
     * pos-1, inclusive.
     * @param contigIndex
     * @param pos
     * @param call
     */
    /*
    private void emitNoCallsUpTo(int contigIndex, long pos, CallType call) {

        if ( contigIndex < currentGenotypeInterval.getContigIndex() ||
             contigIndex == currentGenotypeInterval.getContigIndex() && pos <= currentGenotypeInterval.getStart() ) return;

        if ( contigIndex == currentGenotypeInterval.getContigIndex() && pos >= currentGenotypeInterval.getStart() ) {
            for ( long p = lastGenotypedPosition+1; p < pos; p++ ) {

            }
        }
        while( currentGenotypeInterval != null ) {

            while ( )
        if ( genotypeIntervalIterator.hasNext() ) {
            currentGenotypeInterval = genotypeIntervalIterator.next() ;
            if ( pastInterval(p,currentGenotypeInterval) ) {
                // if we are about to jump over the whole next interval, we need to emit NO_COVERAGE calls there!
                emitNoCoverageCalls(currentGenotypeInterval);
            }
        } else {
            currentGenotypeInterval = null;
        }
        }
    }
*/
    
   /** Output indel calls up to the specified position and shift the window: after this method is executed, the
    * first element of the window maps onto 'position', if possible, or at worst a few bases to the left of 'position' if we may need more
    * reads to get full NQS-style statistics for an indel in the close proximity of 'position'.
    *
    * @param position
    */
   private void emit(long position, boolean force) {

            long adjustedPosition = adjustPosition(position);

            if ( adjustedPosition == -1 ) {
                // failed to find appropriate shift position, the data are probably to messy anyway so we drop them altogether
                normal_context.shift((int)(position-normal_context.getStart()));
                return;
            }
            long move_to = adjustedPosition;

            for ( int pos = normal_context.getStart() ; pos < Math.min(adjustedPosition,normal_context.getStop()+1) ; pos++ ) {

                boolean genotype = false;
                // first let's see if we need to genotype current position:

                final long p = pos - 1; // our internally used positions (pos) are +1 compared to external format spec (e.g. vcf)

                if ( pos <= lastGenotypedPosition ) continue;

                while ( currentGenotypeInterval != null ) {

                    // if we did not even reach next interval yet, no genotyping at current position:
                    if ( location.getContigIndex() < currentGenotypeInterval.getContigIndex() ||
                         location.getContigIndex() == currentGenotypeInterval.getContigIndex() &&
                                 p < currentGenotypeInterval.getStart() ) break;
                    if ( pastInterval(p, currentGenotypeInterval) ) {
                        // we are past current genotyping interval, so we are done with it; let's load next interval:
                        currentGenotypeInterval = genotypeIntervalIterator.hasNext() ? genotypeIntervalIterator.next() : null;
                        continue; // re-enter the loop to check against the interval we just loaded
                    }

                    // we reach this point only if p is inside current genotyping interval; set the flag and bail out:
                    genotype = true;
                    break;
                }

//                if ( DEBUG ) System.out.println("DEBUG>> pos="+pos +"; genotyping interval="+currentGenotypeInterval+"; genotype="+genotype);

                if ( normal_context.indelsAt(pos).size() == 0 && ! genotype ) continue;

                IndelPrecall normalCall = new IndelPrecall(normal_context,pos,NQS_WIDTH);

                if ( normalCall.getCoverage() < minCoverage && ! genotype ) {
                    if ( DEBUG ) {
                        System.out.println("DEBUG>> Indel at "+pos+"; coverare in normal="+normalCall.getCoverage()+" (SKIPPED)");
                    }
                    continue; // low coverage
                }

                if ( DEBUG ) System.out.println("DEBUG>> "+(normalCall.getAllVariantCount() == 0?"No Indel":"Indel")+" at "+pos);

                long left = Math.max( pos-NQS_WIDTH, normal_context.getStart() );
                long right = pos+( normalCall.getVariant() == null ? 0 : normalCall.getVariant().lengthOnRef())+NQS_WIDTH-1;

                if ( right >= adjustedPosition && ! force) {
                    // we are not asked to force-shift, and there is more coverage around the current indel that we still need to collect

                    // we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
                    // instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
                    move_to = adjustPosition(left);
                    if ( move_to == -1 ) {
                        // failed to find appropriate shift position, the data are probably to messy anyway so we drop them altogether
                        normal_context.shift((int)(adjustedPosition-normal_context.getStart()));
                        return;
                    }
                    if ( DEBUG ) System.out.println("DEBUG>> waiting for coverage; actual shift performed to "+ move_to);
                    break;
                }

                // if indel is too close to the end of the window but we need to emit anyway (force-shift), adjust right:
                if ( right > normal_context.getStop() ) right = normal_context.getStop();

    //            location = getToolkit().getGenomeLocParser().setStart(location,pos);
    //            location = getToolkit().getGenomeLocParser().setStop(location,pos); // retrieve annotation data

                location = getToolkit().getGenomeLocParser().createGenomeLoc(location.getContig(), pos);

                boolean haveCall = normalCall.isCall(); // cache the value

                if ( haveCall || genotype) {
                    if ( haveCall ) normalCallsMade++;
                    printVCFLine(vcf_writer,normalCall);
                    if ( bedWriter != null ) normalCall.printBedLine(bedWriter);
                    if ( verboseWriter != null ) printVerboseLine(verboseWriter, normalCall);
                    lastGenotypedPosition = pos;
                }

                normal_context.indelsAt(pos).clear();
                    // we dealt with this indel; don't want to see it again
                    // (we might otherwise in the case when 1) there is another indel that follows
                    // within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)

//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
            }

            if ( DEBUG ) System.out.println("DEBUG>> Actual shift to " + move_to + " ("+adjustedPosition+")");
            normal_context.shift((int)(move_to - normal_context.getStart() ) );
    }

    /** A shortcut. Returns true if we got indels within the specified interval in single and only window context
     * (for single-sample calls) or in either of the two window contexts (for two-sample/somatic calls)
     *
     */
    private boolean indelsPresentInInterval(long start, long stop) {
        if ( tumor_context == null ) return  normal_context.hasIndelsInInterval(start,stop);
        return tumor_context.hasIndelsInInterval(start,stop) ||
              normal_context.hasIndelsInInterval(start,stop);
    }
        /** Takes the position, to which window shift is requested, and tries to adjust it in such a way that no NQS window is broken.
         * Namely, this method checks, iteratively, if there is an indel within NQS_WIDTH bases ahead of initially requested or adjusted 
         * shift position. If there is such an indel,
         * then shifting to that position would lose some or all NQS-window bases to the left of the indel (since it's not going to be emitted
         * just yet). Instead, this method tries to readjust the shift position leftwards so that full NQS window to the left of the next indel
         * is preserved. This method tries thie strategy 4 times (so that it would never walk away too far to the left), and if it fails to find
         * an appropriate adjusted shift position (which could happen if there are many indels following each other at short intervals), it will give up, 
         * go back to the original requested shift position and try finding the first shift poisition that has no indel associated with it.
         */

    private long adjustPosition(long request) {
        long initial_request = request;
        int attempts = 0;
        boolean failure = false;
        while ( indelsPresentInInterval(request,request+NQS_WIDTH)  ) {
            request -= NQS_WIDTH;
            if ( DEBUG ) System.out.println("DEBUG>> indel observations present within "+NQS_WIDTH+" bases ahead. Resetting shift to "+request);
            attempts++;
            if ( attempts == 4 ) {
                if ( DEBUG ) System.out.println("DEBUG>> attempts to preserve full NQS window failed; now trying to find any suitable position.") ;
                failure = true;
                break;
            }
        }

        if ( failure ) {
            // we tried 4 times but did not find a good shift position that would preserve full nqs window
            // around all indels. let's fall back and find any shift position as long and there's no indel at the very
            // first position after the shift (this is bad for other reasons); if it breaks a nqs window, so be it
            request = initial_request;
            attempts = 0;
            while ( indelsPresentInInterval(request,request+1) ) {
                request--;
                if ( DEBUG ) System.out.println("DEBUG>> indel observations present within "+NQS_WIDTH+" bases ahead. Resetting shift to "+request);
                attempts++;
                if ( attempts == 50 ) {
                    System.out.println("WARNING: Indel at every position in the interval "+refName+":"+request+"-"+initial_request+
                            ". Can not find a break to shift context window to; no calls will be attempted in the current window.");
                    return -1;
                }
            }
        }
        if ( DEBUG ) System.out.println("DEBUG>> Found acceptable target position "+request);
        return request;
    }

    /** Output somatic indel calls up to the specified position and shift the coverage array(s): after this method is executed
     * first elements of the coverage arrays map onto 'position', or a few bases prior to the specified position
     * if there is an indel in close proximity to 'position' so that we may get more coverage around it later.
     *
     * @param position
     */
    private void emit_somatic(long position, boolean force) {

        long adjustedPosition = adjustPosition(position);
        if ( adjustedPosition == -1 ) {
            // failed to find appropriate shift position, the data are probably to messy anyway so we drop them altogether
            normal_context.shift((int)(position-normal_context.getStart()));
            tumor_context.shift((int)(position-tumor_context.getStart()));
            return;
        }
        long move_to = adjustedPosition;

        if ( DEBUG ) System.out.println("DEBUG>> Emitting in somatic mode up to "+position+" force shift="+force+" current window="+tumor_context.getStart()+"-"+tumor_context.getStop());

        for ( int pos = tumor_context.getStart() ; pos < Math.min(adjustedPosition,tumor_context.getStop()+1) ; pos++ ) {

            boolean genotype = false;
             // first let's see if we need to genotype current position:

             final long p = pos - 1; // our internally used positions (pos) are +1 compared to external format spec (e.g. vcf)

             if ( pos <= lastGenotypedPosition ) continue;

             while ( currentGenotypeInterval != null ) {

                 // if we did not even reach next interval yet, no genotyping at current position:
                 if ( location.getContigIndex() < currentGenotypeInterval.getContigIndex() ||
                      location.getContigIndex() == currentGenotypeInterval.getContigIndex() &&
                              p < currentGenotypeInterval.getStart() ) break;
                 if ( pastInterval(p, currentGenotypeInterval) ) {
                     // we are past current genotyping interval, so we are done with it; let's load next interval:
                     currentGenotypeInterval = genotypeIntervalIterator.hasNext() ? genotypeIntervalIterator.next() : null;
                     continue; // re-enter the loop to check against the interval we just loaded
                 }

                 // we reach tjis point only if p is inside current genotyping interval; set the flag and bail out:
                 genotype = true;
                 break;
             }
//            if ( DEBUG) System.out.println("DEBUG>> pos="+pos +"; genotyping interval="+currentGenotypeInterval+"; genotype="+genotype);

            if ( tumor_context.indelsAt(pos).size() == 0 && ! genotype ) continue; // no indels in tumor

            if ( DEBUG && genotype ) System.out.println("DEBUG>> Genotyping requested at "+pos);

            IndelPrecall tumorCall = new IndelPrecall(tumor_context,pos,NQS_WIDTH);
            IndelPrecall normalCall = new IndelPrecall(normal_context,pos,NQS_WIDTH);

            if ( tumorCall.getCoverage() < minCoverage && ! genotype ) {
                if ( DEBUG ) {
                    System.out.println("DEBUG>> Indel in tumor at "+pos+"; coverare in tumor="+tumorCall.getCoverage()+" (SKIPPED)");
                }
                continue; // low coverage
            }
            if ( normalCall.getCoverage() < minNormalCoverage && ! genotype ) {
                if ( DEBUG ) {
                    System.out.println("DEBUG>> Indel in tumor at "+pos+"; coverare in normal="+normalCall.getCoverage()+" (SKIPPED)");
                }
                continue; // low coverage
            }

            if ( DEBUG ) {
                System.out.print("DEBUG>> "+(tumorCall.getAllVariantCount() == 0?"No Indel":"Indel")+" in tumor, ");
                System.out.print("DEBUG>> "+(normalCall.getAllVariantCount() == 0?"No Indel":"Indel")+" in normal at "+pos);
            }

            long left = Math.max( pos-NQS_WIDTH, tumor_context.getStart() );
            long right = pos+ ( tumorCall.getVariant() == null ? 0 : tumorCall.getVariant().lengthOnRef() )+NQS_WIDTH-1;

            if ( right >= adjustedPosition && ! force) {
                // we are not asked to force-shift, and there is more coverage around the current indel that we still need to collect

                // we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
                // instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
                move_to = adjustPosition(left);
                if ( move_to == -1 ) {
                    // failed to find appropriate shift position, the data are probably to messy anyway so we drop them altogether
                    normal_context.shift((int)(adjustedPosition-normal_context.getStart()));
                    tumor_context.shift((int)(adjustedPosition-tumor_context.getStart()));
                    return;
                }
                if ( DEBUG ) System.out.println("DEBUG>> waiting for coverage; actual shift performed to "+ move_to);
                break;
            }

            if ( right > tumor_context.getStop() ) right = tumor_context.getStop(); // if indel is too close to the end of the window but we need to emit anyway (force-shift), adjust right

//            location = getToolkit().getGenomeLocParser().setStart(location,pos);
//            location = getToolkit().getGenomeLocParser().setStop(location,pos); // retrieve annotation data

            location = getToolkit().getGenomeLocParser().createGenomeLoc(location.getContig(),pos); // retrieve annotation data

            boolean haveCall = tumorCall.isCall(); // cache the value

            if ( haveCall || genotype ) {
                if ( haveCall ) tumorCallsMade++;

                printVCFLine(vcf_writer,normalCall,tumorCall);

                if ( bedWriter != null ) tumorCall.printBedLine(bedWriter);

                if ( verboseWriter != null ) printVerboseLine(verboseWriter, normalCall, tumorCall );
                lastGenotypedPosition = pos;
            }
            tumor_context.indelsAt(pos).clear();
            normal_context.indelsAt(pos).clear();
                // we dealt with this indel; don't want to see it again
                // (we might otherwise in the case when 1) there is another indel that follows
                // within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)

//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
        }

        if ( DEBUG ) System.out.println("DEBUG>> Actual shift to " + move_to + " ("+adjustedPosition+")");
        tumor_context.shift((int)(move_to - tumor_context.getStart() ) );
        normal_context.shift((int)(move_to - normal_context.getStart() ) );
    }

    private String makeFullRecord(IndelPrecall normalCall, IndelPrecall tumorCall) {
        StringBuilder fullRecord = new StringBuilder();
        if ( tumorCall.getVariant() != null || normalCall.getVariant() == null) {
            fullRecord.append(tumorCall.makeEventString());
        } else {
            fullRecord.append(normalCall.makeEventString());            
        }
        fullRecord.append('\t');
        fullRecord.append(normalCall.makeStatsString("N_"));
        fullRecord.append('\t');
        fullRecord.append(tumorCall.makeStatsString("T_"));
        fullRecord.append('\t');
        return fullRecord.toString();
    }

    private String makeFullRecord(IndelPrecall normalCall) {
        StringBuilder fullRecord = new StringBuilder();
        fullRecord.append(normalCall.makeEventString());
        fullRecord.append('\t');
        fullRecord.append(normalCall.makeStatsString(""));
        fullRecord.append('\t');
        return fullRecord.toString();
    }

    private String getAnnotationString(RODRecordList ann) {
        if ( ann == null ) return annGenomic;
        else {
            StringBuilder b = new StringBuilder();

            if ( RefSeqFeature.isExon(ann) ) {
                if ( RefSeqFeature.isCodingExon(ann) ) b.append(annCoding); // both exon and coding = coding exon sequence
                else b.append(annUTR); // exon but not coding = UTR
            } else {
                if ( RefSeqFeature.isCoding(ann) ) b.append(annIntron); // not in exon, but within the coding region = intron
                else b.append(annUnknown); // we have no idea what this is. this may actually happen when we have a fully non-coding exon...
            }
            b.append('\t');
            b.append(((Transcript)ann.get(0).getUnderlyingObject()).getGeneName()); // there is at least one transcript in the list, guaranteed
//			while ( it.hasNext() ) { //
//				t.getGeneName()
//			}
            return b.toString();
        }

    }

    public void printVerboseLine(Writer verboseWriter, IndelPrecall normalCall) {
        RODRecordList annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(location));
        String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotationList));

        StringBuilder fullRecord = new StringBuilder();
        fullRecord.append(makeFullRecord(normalCall));
        fullRecord.append(annotationString);
        if ( ! normalCall.isCall() && normalCall.getVariant() != null ) fullRecord.append("\tFILTERED_NOCALL");
        try {
            verboseWriter.write(fullRecord.toString());
            verboseWriter.write('\n');
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(verboseOutput, "Write failed", e);
        }

    }


    public void printVerboseLine(Writer verboseWriter, IndelPrecall normalCall, IndelPrecall tumorCall) {
        RODRecordList annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(location));
        String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotationList));

        StringBuilder fullRecord = new StringBuilder();
        fullRecord.append(makeFullRecord(normalCall,tumorCall));

        if ( normalCall.getVariant() == null && tumorCall.getVariant() == null ) {
            // did not observe anything
            if ( normalCall.getCoverage() >= minNormalCoverage && tumorCall.getCoverage() >= minCoverage ) fullRecord.append("REFERENCE");
            else {
                if ( tumorCall.getCoverage() >= minCoverage ) fullRecord.append("REFERENCE"); // no coverage in normal but nothing in tumor
                else {
                    // no coverage in tumor; if we have no coverage in normal, it can be anything; if we do have coverage in normal,
                    // this still could be a somatic event. so either way it is 'unknown'
                    fullRecord.append("UNKNOWN");
                }
            }

        }

        if ( normalCall.getVariant() == null && tumorCall.getVariant() != null ) {
            // looks like somatic call
            if ( normalCall.getCoverage() >= minNormalCoverage ) fullRecord.append("SOMATIC"); // we confirm there is nothing in normal
            else {
                // low coverage in normal
                fullRecord.append("EVENT_T"); // no coverage in normal, no idea whether it is germline or somatic
            }
        }

        if ( normalCall.getVariant() != null && tumorCall.getVariant() == null ) {
            // it's likely germline (with missing observation in tumor - maybe loh?
            if ( tumorCall.getCoverage() >= minCoverage ) fullRecord.append("GERMLINE_LOH"); // we confirm there is nothing in tumor
            else {
                // low coverage in tumor, maybe we missed the event
                fullRecord.append("GERMLINE"); // no coverage in tumor but we already saw it in normal...
            }
        }

        if ( normalCall.getVariant() != null && tumorCall.getVariant() != null ) {
            // events in both T/N, got to be germline!
            fullRecord.append("GERMLINE"); 
        }


        fullRecord.append('\t');
        fullRecord.append(annotationString);

        if ( ! tumorCall.isCall() && tumorCall.getVariant() != null ) fullRecord.append("\tFILTERED_NOCALL");

        try {
            verboseWriter.write(fullRecord.toString());
            verboseWriter.write('\n');
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(verboseOutput, "Write failed", e);
        }
    }

    public void printVCFLine(VCFWriter vcf, IndelPrecall call) {

        long start = call.getPosition()-1;
        // If the beginning of the chromosome is deleted (possible, however unlikely), it's unclear how to proceed.
        // The suggestion is instead of putting the base before the indel, to put the base after the indel.
        // For now, just don't print out that site.
        if ( start == 0 )
            return;

        long stop = start;

        List<Allele> alleles = new ArrayList<Allele>(2); // actual observed (distinct!) alleles at the site
        List<Allele> homref_alleles = null; // when needed, will contain two identical copies of ref allele - needed to generate hom-ref genotype


        if ( call.getVariant() == null ) {
            // we will need to cteate genotype with two (hom) ref alleles (below).
            // we can not use 'alleles' list here, since that list is supposed to contain
            // only *distinct* alleles observed at the site or VCFContext will frown upon us...
            alleles.add( Allele.create(refBases[(int)start-1],true) );
            homref_alleles = new ArrayList<Allele>(2);
            homref_alleles.add( alleles.get(0));
            homref_alleles.add( alleles.get(0));
        } else {
            // we always create alt allele when we observe anything but the ref, even if it is not a call!
            // (Genotype will tell us whether it is an actual call or not!)
            int event_length = call.getVariant().lengthOnRef();
            if ( event_length < 0 ) event_length = 0;
            fillAlleleList(alleles,call);
            stop += event_length;
        }

        Map<String,Genotype> genotypes = new HashMap<String,Genotype>();

        for ( String sample : normalSamples ) {

            Map<String,?> attrs = call.makeStatsAttributes(null);

            if ( call.isCall() ) // we made a call - put actual het genotype here:
                genotypes.put(sample,new Genotype(sample,alleles,Genotype.NO_NEG_LOG_10PERROR,null,attrs,false));
            else // no call: genotype is ref/ref (but alleles still contain the alt if we observed anything at all) 
                genotypes.put(sample,new Genotype(sample, homref_alleles,Genotype.NO_NEG_LOG_10PERROR,null,attrs,false));

        }
        Set<String> filters = null;
        if ( call.getVariant() != null && ! call.isCall() ) {
            filters = new HashSet<String>();
            filters.add("NoCall");
        }
        VariantContext vc = new VariantContext("IGv2_Indel_call", refName, start, stop, alleles, genotypes,
            -1.0 /* log error */,  filters, null, refBases[(int)start-1]);
        vcf.add(vc);
    }

    /** Fills l with appropriate alleles depending on whether call is insertion or deletion
     * (l MUST have a variant or this method will crash). It is guaranteed that the *first* allele added
     * to the list is ref, and the next one is alt.
     * @param l
     * @param call
     */
    private void fillAlleleList(List<Allele> l, IndelPrecall call) {
        int event_length = call.getVariant().lengthOnRef();
        if ( event_length == 0 ) { // insertion

            l.add( Allele.create(Allele.NULL_ALLELE_STRING,true) );
            l.add( Allele.create(call.getVariant().getBases(), false ));

        } else { //deletion:
            l.add( Allele.create(call.getVariant().getBases(), true ));
            l.add( Allele.create(Allele.NULL_ALLELE_STRING,false) );
        }
    }

    public void printVCFLine(VCFWriter vcf, IndelPrecall nCall, IndelPrecall tCall) {

        long start = tCall.getPosition()-1;
        long stop = start;

        // If the beginning of the chromosome is deleted (possible, however unlikely), it's unclear how to proceed.
        // The suggestion is instead of putting the base before the indel, to put the base after the indel.
        // For now, just don't print out that site.
        if ( start == 0 )
            return;

        Map<String,Object> attrsNormal = nCall.makeStatsAttributes(null);
        Map<String,Object> attrsTumor = tCall.makeStatsAttributes(null);

        Map<String,Object> attrs = new HashMap();

        boolean isSomatic = false;
        if ( nCall.getCoverage() >= minNormalCoverage && nCall.getVariant() == null && tCall.getVariant() != null ) {
            isSomatic = true;
            attrs.put(VCFConstants.SOMATIC_KEY,true);
        }
        List<Allele> alleles = new ArrayList<Allele>(2); // all alleles at the site
 //       List<Allele> normal_alleles = null; // all alleles at the site
        List<Allele> homRefAlleles = null;

//        if ( nCall.getVariant() == null || tCall.getVariant() == null ) {
        homRefAlleles = new ArrayList<Allele>(2) ; // we need this for somatic calls (since normal is ref-ref), and also for no-calls
//        }
        boolean homRefT = ( tCall.getVariant() == null );
        boolean homRefN = ( nCall.getVariant() == null );
        if ( tCall.getVariant() == null && nCall.getVariant() == null) {
            // no indel at all  ; create base-representation ref/ref alleles for genotype construction
            alleles.add( Allele.create(refBases[(int)start-1],true) );
        } else {
            // we got indel(s)
            int event_length = 0;
            if ( tCall.getVariant() != null ) {
                // indel in tumor
                event_length = tCall.getVariant().lengthOnRef();
                fillAlleleList(alleles, tCall);
            } else {
                event_length = nCall.getVariant().lengthOnRef();
                fillAlleleList(alleles, nCall);
            }
            if ( event_length > 0 ) stop += event_length;
        }
        homRefAlleles.add( alleles.get(0));
        homRefAlleles.add( alleles.get(0));

        Map<String,Genotype> genotypes = new HashMap<String,Genotype>();

        for ( String sample : normalSamples ) {
            genotypes.put(sample,new Genotype(sample, homRefN ? homRefAlleles : alleles,Genotype.NO_NEG_LOG_10PERROR,null,attrsNormal,false));
        }

        for ( String sample : tumorSamples ) {
            genotypes.put(sample,new Genotype(sample, homRefT ? homRefAlleles : alleles,Genotype.NO_NEG_LOG_10PERROR,null,attrsTumor,false) );
        }

        Set<String> filters = null;
        if ( tCall.getVariant() != null && ! tCall.isCall() ) {
            filters = new HashSet<String>();
            filters.add("NoCall");
        }
        if ( nCall.getCoverage() < minNormalCoverage ) {
            if ( filters == null ) filters = new HashSet<String>();
            filters.add("NCov");
        }
        if ( tCall.getCoverage() < minCoverage ) {
            if ( filters == null ) filters = new HashSet<String>();
            filters.add("TCov");
        }

        VariantContext vc = new VariantContext("IGv2_Indel_call", refName, start, stop, alleles, genotypes,
            -1.0 /* log error */, filters, attrs, refBases[(int)start-1]);
        vcf.add(vc);
    }

    @Override
    public void onTraversalDone(Integer result) {
        if ( DEBUG ) {
            System.out.println("DEBUG>> Emitting last window at "+normal_context.getStart()+"-"+normal_context.getStop());
        }
        if ( call_somatic ) emit_somatic(1000000000, true);
        else emit(1000000000,true); // emit everything we might have left

        if ( metricsWriter != null ) {
            metricsWriter.println(String.format("Normal calls made     %d", normalCallsMade));
            metricsWriter.println(String.format("Tumor calls made      %d", tumorCallsMade));
            metricsWriter.close();
        }

        try {
            if ( bedWriter != null ) bedWriter.close();
            if ( verboseWriter != null ) verboseWriter.close();
        } catch (IOException e) {
            System.out.println("Failed to close output BED file gracefully, data may be lost");
            e.printStackTrace();
        }
        super.onTraversalDone(result);
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        if ( value == -1 ) {
            onTraversalDone(sum);
            System.exit(1);
        }
        sum += value;
        return sum;
    }

    @Override
    public Integer reduceInit() {
        return new Integer(0);
    }


        static class IndelVariant {
            public static enum Type { I, D};
            private String bases;
            private Type type;
            private ArrayList<Integer> fromStartOffsets = null;
            private ArrayList<Integer> fromEndOffsets = null;

            private Set<ExpandedSAMRecord> reads = new HashSet<ExpandedSAMRecord>(); // keep track of reads that have this indel
            private Set<String> samples = new HashSet<String>();   // which samples had the indel described by this object

            public IndelVariant(ExpandedSAMRecord read , Type type, String bases) {
                this.type = type;
                this.bases = bases.toUpperCase();
                addObservation(read);
                fromStartOffsets = new ArrayList<Integer>();
                fromEndOffsets = new ArrayList<Integer>();
            }

            /** Adds another observation for the current indel. It is assumed that the read being registered
             * does contain the observation, no checks are performed. Read's sample is added to the list of samples
             * this indel was observed in as well.
             * @param read
             */
            public void addObservation(ExpandedSAMRecord read) {
                if ( reads.contains(read) ) {
                    //TODO fix CleanedReadInjector and reinstate exception here: duplicate records may signal a problem with the bam
                    // seeing the same read again can mean only one thing: the input bam file is corrupted and contains
                    // duplicate records. We KNOW that this may happen for the time being due to bug in CleanedReadInjector
                    // so this is a short-term patch: don't cry, but just ignore the duplicate record

                    //throw new StingException("Attempting to add indel observation that was already registered");
                    return;
                }
                reads.add(read);
                String sample = null;
                if ( read.getSAMRecord().getReadGroup() != null ) sample = read.getSAMRecord().getReadGroup().getSample();
                if ( sample != null ) samples.add(sample);
            }


            /** Returns length of the event on the reference (number of deleted bases
             * for deletions, -1 for insertions.
             * @return
             */
            public int lengthOnRef() {
                if ( type == Type.D ) return bases.length();
                else return 0;
            }


            public void addSample(String sample) {
                if ( sample != null )
                samples.add(sample);
            }

            public void addReadPositions(int fromStart, int fromEnd) {
                fromStartOffsets.add(fromStart);
                fromEndOffsets.add(fromEnd);
            }

            public List<Integer> getOffsetsFromStart() { return fromStartOffsets ; }
            public List<Integer> getOffsetsFromEnd() { return fromEndOffsets; }

            public String getSamples() {
                StringBuffer sb = new StringBuffer();
                Iterator<String> i = samples.iterator();
                while ( i.hasNext() ) {
                    sb.append(i.next());
                    if ( i.hasNext() )
                        sb.append(",");
                }
                return sb.toString();
            }

            public Set<ExpandedSAMRecord> getReadSet() { return reads; }

            public int getCount() { return reads.size(); }

            public String getBases() { return bases; }

            public Type getType() { return type; }

            @Override
            public boolean equals(Object o) {
                if ( ! ( o instanceof IndelVariant ) ) return false;
                IndelVariant that = (IndelVariant)o;
                return ( this.type == that.type && this.bases.equals(that.bases) );
            }

            public boolean equals(Type type, String bases) {
                return ( this.type == type && this.bases.equals(bases.toUpperCase()) );
            }
        }

    /**
     * Utility class that encapsulates the logic related to collecting all the stats and counts required to
     * make (or discard) a call, as well as the calling heuristics that uses those data.
      */
    class IndelPrecall {
//        private boolean DEBUG = false;
        private int NQS_MISMATCH_CUTOFF = 1000000;
        private double AV_MISMATCHES_PER_READ = 1.5;

        private int nqs = 0;
        private IndelVariant consensus_indel = null; // indel we are going to call
        private long pos = -1 ; // position on the ref
        private int total_coverage = 0; // total number of reads overlapping with the event
        private int consensus_indel_count = 0; // number of reads, in which consensus indel was observed
        private int all_indel_count = 0 ; // number of reads, in which any indel was observed at current position

        private int total_mismatches_in_nqs_window = 0; // total number of mismatches in the nqs window around the indel
        private int total_bases_in_nqs_window = 0; // total number of bases in the nqs window (some reads may not fully span the window so it's not coverage*nqs_size)
        private int total_base_qual_in_nqs_window = 0; // sum of qualitites of all the bases in the nqs window
        private int total_mismatching_base_qual_in_nqs_window = 0; // sum of qualitites of all mismatching bases in the nqs window

        private int indel_read_mismatches_in_nqs_window = 0;   // mismatches inside the nqs window in indel-containing reads only
        private int indel_read_bases_in_nqs_window = 0;  // number of bases in the nqs window from indel-containing reads only
        private int indel_read_base_qual_in_nqs_window = 0; // sum of qualitites of bases in nqs window from indel-containing reads only
        private int indel_read_mismatching_base_qual_in_nqs_window = 0; // sum of qualitites of mismatching bases in the nqs window from indel-containing reads only


        private int consensus_indel_read_mismatches_in_nqs_window = 0; // mismatches within the nqs window from consensus indel reads only
        private int consensus_indel_read_bases_in_nqs_window = 0;  // number of bases in the nqs window from consensus indel-containing reads only
        private int consensus_indel_read_base_qual_in_nqs_window = 0; // sum of qualitites of bases in nqs window from consensus indel-containing reads only
        private int consensus_indel_read_mismatching_base_qual_in_nqs_window = 0; // sum of qualitites of mismatching bases in the nqs window from consensus indel-containing reads only


        private double consensus_indel_read_total_mm = 0.0; // sum of all mismatches in reads that contain consensus indel
        private double all_indel_read_total_mm = 0.0; // sum of all mismatches in reads that contain any indel at given position
        private double all_read_total_mm = 0.0; // sum of all mismatches in all reads

        private double consensus_indel_read_total_mapq = 0.0; // sum of mapping qualitites of all reads with consensus indel
        private double all_indel_read_total_mapq = 0.0 ; // sum of mapping qualitites of all reads with (any) indel at current position
        private double all_read_total_mapq = 0.0; // sum of all mapping qualities of all reads

        private PrimitivePair.Int consensus_indel_read_orientation_cnt = new PrimitivePair.Int();
        private PrimitivePair.Int all_indel_read_orientation_cnt = new PrimitivePair.Int();
        private PrimitivePair.Int all_read_orientation_cnt = new PrimitivePair.Int();

        private int from_start_median = 0;
        private int from_start_mad = 0;
        private int from_end_median = 0;
        private int from_end_mad = 0;

        /** Makes an empty call (no-call) with all stats set to 0
         *
         * @param position
         */
        public IndelPrecall(long position) {
            this.pos = position;
        }

        public IndelPrecall(WindowContext context, long position, int nqs_width) {
            this.pos = position;
            this.nqs = nqs_width;
            total_coverage = context.coverageAt(pos,true);
            List<IndelVariant> variants = context.indelsAt(pos);
            findConsensus(variants);

            // pos is the first base after the event: first deleted base or first base after insertion.
            // hence, [pos-nqs, pos+nqs-1] (inclusive) is the window with nqs bases on each side of a no-event or an insertion
            // and [pos-nqs, pos+Ndeleted+nqs-1] is the window with nqs bases on each side of a deletion.
            // we initialize the nqs window for no-event/insertion case
            long left = Math.max( pos-nqs, context.getStart() );
            long right = Math.min(pos+nqs-1, context.getStop());
//if ( pos == 3534096 ) System.out.println("pos="+pos +" total reads: "+context.getReads().size());
            Iterator<ExpandedSAMRecord> read_iter = context.getReads().iterator();


            while ( read_iter.hasNext() ) {
                ExpandedSAMRecord rec = read_iter.next();
                SAMRecord read = rec.getSAMRecord();
                byte[] flags = rec.getExpandedMMFlags();
                byte[] quals = rec.getExpandedQuals();
                int mm = rec.getMMCount();


                if( read.getAlignmentStart() > pos || read.getAlignmentEnd() < pos ) continue;

                long local_right = right; // end of nqs window for this particular read. May need to be advanced further right
                                          // if read has a deletion. The gap in the middle of nqs window will be skipped
                                          // automatically since flags/quals are set to -1 there

                boolean read_has_a_variant = false;
                boolean read_has_consensus = ( consensus_indel!= null && consensus_indel.getReadSet().contains(rec) );
                for ( IndelVariant v : variants ) {
                    if ( v.getReadSet().contains(rec) ) {
                        read_has_a_variant = true;
                        local_right += v.lengthOnRef();
                        break;
                    }
                }

                if ( read_has_consensus ) {
                    consensus_indel_read_total_mm += mm;
                    consensus_indel_read_total_mapq += read.getMappingQuality();
                    if ( read.getReadNegativeStrandFlag() ) consensus_indel_read_orientation_cnt.second++;
                    else consensus_indel_read_orientation_cnt.first++;
                }
                if ( read_has_a_variant ) {
                    all_indel_read_total_mm += mm;
                    all_indel_read_total_mapq += read.getMappingQuality();
                    if ( read.getReadNegativeStrandFlag() ) all_indel_read_orientation_cnt.second++;
                    else all_indel_read_orientation_cnt.first++;
                }

                all_read_total_mm+= mm;
                all_read_total_mapq += read.getMappingQuality();
                if ( read.getReadNegativeStrandFlag() ) all_read_orientation_cnt.second++;
                else all_read_orientation_cnt.first++;

                for ( int pos_in_flags = Math.max((int)(left - read.getAlignmentStart()),0);
                      pos_in_flags <= Math.min((int)local_right-read.getAlignmentStart(),flags.length - 1);
                       pos_in_flags++) {

                        if ( flags[pos_in_flags] == -1 ) continue; // gap (deletion), skip it; we count only bases aligned to the ref
                        total_bases_in_nqs_window++;
                        if ( read_has_consensus ) consensus_indel_read_bases_in_nqs_window++;
                        if ( read_has_a_variant ) indel_read_bases_in_nqs_window++;

                        if ( quals[pos_in_flags] != -1 ) {

                            total_base_qual_in_nqs_window += quals[pos_in_flags];
                            if ( read_has_a_variant ) indel_read_base_qual_in_nqs_window += quals[pos_in_flags];
                            if ( read_has_consensus ) consensus_indel_read_base_qual_in_nqs_window += quals[pos_in_flags];
                        }

                        if ( flags[pos_in_flags] == 1 ) { // it's a mismatch
                            total_mismatches_in_nqs_window++;
                            total_mismatching_base_qual_in_nqs_window += quals[pos_in_flags];

                            if ( read_has_consensus ) {
                                consensus_indel_read_mismatches_in_nqs_window++;
                                consensus_indel_read_mismatching_base_qual_in_nqs_window += quals[pos_in_flags];
                            }
                            
                            if ( read_has_a_variant ) {
                                indel_read_mismatches_in_nqs_window++;
                                indel_read_mismatching_base_qual_in_nqs_window += quals[pos_in_flags];
                            }
                        }
                }
//         if ( pos == 3534096 ) {
//             System.out.println(read.getReadName());
//             System.out.println(" cons nqs bases="+consensus_indel_read_bases_in_nqs_window);
//             System.out.println(" qual sum="+consensus_indel_read_base_qual_in_nqs_window);
//         }

            }

            // compute median/mad for offsets from the read starts/ends
            if ( consensus_indel != null ) {
                from_start_median = median(consensus_indel.getOffsetsFromStart()) ;
                from_start_mad = mad(consensus_indel.getOffsetsFromStart(),from_start_median);
                from_end_median = median(consensus_indel.getOffsetsFromEnd()) ;
                from_end_mad = mad(consensus_indel.getOffsetsFromEnd(),from_end_median);   
            }
        }

        /** As a side effect will sort l!
         *
         * @param l
         * @return
         */
        private int median(List<Integer> l) {
            Collections.sort(l);
            int k = l.size()/2;
            return ( l.size() % 2 == 0 ?
                  (l.get(k-1).intValue()+l.get(k).intValue())/2 :
                   l.get(k).intValue());
        }

        private int median(int[] l) {
            Arrays.sort(l);
            int k = l.length/2;
            return ( l.length % 2 == 0 ?
                  (l[k-1]+l[k])/2 :
                   l[k]);
        }

        private int mad(List<Integer> l, int med) {
            int [] diff = new int[l.size()];
            for ( int i = 0; i < l.size(); i++ ) {
                   diff[i] = Math.abs(l.get(i).intValue() - med);
            }
            return median(diff);
        }

        public long getPosition() { return pos; }

        public boolean hasObservation() { return consensus_indel != null; }

        public int getCoverage() { return total_coverage; }

        public double getTotalMismatches() { return all_read_total_mm; }
        public double getConsensusMismatches() { return consensus_indel_read_total_mm; }
        public double getAllVariantMismatches() { return all_indel_read_total_mm; }

        /** Returns average number of mismatches per consensus indel-containing read */
        public double getAvConsensusMismatches() {
            return ( consensus_indel_count != 0 ? consensus_indel_read_total_mm/consensus_indel_count : 0.0 );
        }

        /** Returns average number of mismatches per read across all reads matching the ref (not containing any indel variants) */
        public double getAvRefMismatches() {
            int coverage_ref = total_coverage-all_indel_count;
            return ( coverage_ref != 0 ? (all_read_total_mm - all_indel_read_total_mm )/coverage_ref : 0.0 );
        }

        public PrimitivePair.Int getConsensusStrandCounts() {
            return consensus_indel_read_orientation_cnt;
        }

        public PrimitivePair.Int getRefStrandCounts() {
            return new PrimitivePair.Int(all_read_orientation_cnt.first-all_indel_read_orientation_cnt.first,
                                         all_read_orientation_cnt.second - all_indel_read_orientation_cnt.second);
        }

        /** Returns a sum of mapping qualities of all reads spanning the event. */
        public double getTotalMapq() { return all_read_total_mapq; }

        /** Returns a sum of mapping qualities of all reads, in which the consensus variant is observed. */
        public double getConsensusMapq() { return consensus_indel_read_total_mapq; }

        /** Returns a sum of mapping qualities of all reads, in which any variant is observed at the current event site. */
        public double getAllVariantMapq() { return all_indel_read_total_mapq; }

        /** Returns average mapping quality per consensus indel-containing read. */
        public double getAvConsensusMapq() {
            return ( consensus_indel_count != 0 ? consensus_indel_read_total_mapq/consensus_indel_count : 0.0 );
        }

        /** Returns average number of mismatches per read across all reads matching the ref (not containing any indel variants). */
        public double getAvRefMapq() {
            int coverage_ref = total_coverage-all_indel_count;
            return ( coverage_ref != 0 ? (all_read_total_mapq - all_indel_read_total_mapq )/coverage_ref : 0.0 );
        }

        /** Returns fraction of bases in NQS window around the indel that are mismatches, across all reads,
         * in which consensus indel is observed. NOTE: NQS window for indel containing reads is defined around
         * the indel itself (e.g. for a 10-base deletion spanning [X,X+9], the 5-NQS window is {[X-5,X-1],[X+10,X+15]}
         * */
        public double getNQSConsensusMMRate() {
            if ( consensus_indel_read_bases_in_nqs_window == 0 ) return 0;
            return ((double)consensus_indel_read_mismatches_in_nqs_window)/consensus_indel_read_bases_in_nqs_window;
        }

        /** Returns fraction of bases in NQS window around the indel start position that are mismatches, across all reads
         * that align to the ref (i.e. contain no indel observation at the current position). NOTE: NQS window for ref
         * reads is defined around the event start position, NOT around the actual consensus indel.
         * */
        public double getNQSRefMMRate() {
            int num_ref_bases = total_bases_in_nqs_window - indel_read_bases_in_nqs_window;
            if ( num_ref_bases == 0 ) return 0;
            return ((double)(total_mismatches_in_nqs_window - indel_read_mismatches_in_nqs_window))/num_ref_bases;
        }

        /** Returns average base quality in NQS window around the indel, across all reads,
         * in which consensus indel is observed. NOTE: NQS window for indel containing reads is defined around
         * the indel itself (e.g. for a 10-base deletion spanning [X,X+9], the 5-NQS window is {[X-5,X-1],[X+10,X+15]}
         * */
        public double getNQSConsensusAvQual() {
            if ( consensus_indel_read_bases_in_nqs_window == 0 ) return 0;
            return ((double)consensus_indel_read_base_qual_in_nqs_window)/consensus_indel_read_bases_in_nqs_window;
        }

        /** Returns fraction of bases in NQS window around the indel start position that are mismatches, across all reads
         * that align to the ref (i.e. contain no indel observation at the current position). NOTE: NQS window for ref
         * reads is defined around the event start position, NOT around the actual consensus indel.
         * */
        public double getNQSRefAvQual() {
            int num_ref_bases = total_bases_in_nqs_window - indel_read_bases_in_nqs_window;
            if ( num_ref_bases == 0 ) return 0;
            return ((double)(total_base_qual_in_nqs_window - indel_read_base_qual_in_nqs_window))/num_ref_bases;
        }

        public int getTotalNQSMismatches() { return total_mismatches_in_nqs_window; }

        public int getAllVariantCount() { return all_indel_count; }
        public int getConsensusVariantCount() { return consensus_indel_count; }

//        public boolean failsNQSMismatch() {
//            //TODO wrong fraction: mismatches are counted only in indel-containing reads, but total_coverage is used!
//            return ( indel_read_mismatches_in_nqs_window > NQS_MISMATCH_CUTOFF ) ||
//                    ( indel_read_mismatches_in_nqs_window > total_coverage * AV_MISMATCHES_PER_READ );
//        }

        public IndelVariant getVariant() { return consensus_indel; }

        public boolean isCall() {
            boolean ret =  ( consensus_indel_count >= minIndelCount &&
                    (double)consensus_indel_count > minFraction * total_coverage &&
                    (double) consensus_indel_count > minConsensusFraction*all_indel_count && total_coverage >= minCoverage);
            if ( DEBUG && ! ret ) System.out.println("DEBUG>>  NOT a call: count="+consensus_indel_count+
                        " total_count="+all_indel_count+" cov="+total_coverage+
                " minConsensusF="+((double)consensus_indel_count)/all_indel_count+
                    " minF="+((double)consensus_indel_count)/total_coverage);
            return ret;

        }

        /** Utility method: finds the indel variant with the largest count (ie consensus) among all the observed
         * variants, and sets the counts of consensus observations and all observations of any indels (including non-consensus)
         * @param variants
         * @return
         */
        private void findConsensus(List<IndelVariant> variants) {
            for ( IndelVariant var : variants ) {
                if ( DEBUG ) System.out.println("DEBUG>> Variant "+var.getBases()+" (cnt="+var.getCount()+")");
                int cnt = var.getCount();
                all_indel_count +=cnt;
                if ( cnt > consensus_indel_count ) {
                    consensus_indel = var;
                    consensus_indel_count = cnt;
                }
            }
            if ( DEBUG && consensus_indel != null ) System.out.println("DEBUG>> Returning: "+consensus_indel.getBases()+
                    " (cnt="+consensus_indel.getCount()+") with total count of "+all_indel_count);
        }



        public void printBedLine(Writer bed) {
            int event_length;
            if ( consensus_indel == null ) event_length = 0;
            else {
                event_length = consensus_indel.lengthOnRef();
                if ( event_length < 0 ) event_length = 0;
            }

            StringBuffer message = new StringBuffer();
            message.append(refName+"\t"+(pos-1)+"\t");
            message.append((pos-1+event_length)+"\t");
            if ( consensus_indel != null ) {
                message.append((event_length>0? "-":"+")+consensus_indel.getBases());
            } else {
                message.append('.');
            }
            message.append(":"+all_indel_count+"/"+total_coverage);
            try {
                bed.write(message.toString()+"\n");
            } catch (IOException e) {
               throw new UserException.CouldNotCreateOutputFile(bedOutput, "Error encountered while writing into output BED file", e);
            }
        }

        public String makeEventString() {
            int event_length;
            if ( consensus_indel == null ) event_length = 0;
            else {
                event_length = consensus_indel.lengthOnRef();
                if ( event_length < 0 ) event_length = 0;
            }
            StringBuffer message = new StringBuffer();
            message.append(refName);
            message.append('\t');
            message.append(pos-1);
            message.append('\t');
            message.append(pos-1+event_length);
            message.append('\t');
            if ( consensus_indel != null ) {
                message.append((event_length>0?'-':'+'));
                message.append(consensus_indel.getBases());
            } else {
                message.append('.');
            }
            return message.toString();
        }

        public String makeStatsString(String prefix) {
             StringBuilder message = new StringBuilder();
             message.append(prefix+"OBS_COUNTS[C/A/T]:"+getConsensusVariantCount()+"/"+getAllVariantCount()+"/"+getCoverage());
             message.append('\t');
             message.append(prefix+"AV_MM[C/R]:"+String.format("%.2f/%.2f",getAvConsensusMismatches(),
                                 getAvRefMismatches()));
             message.append('\t');
             message.append(prefix+"AV_MAPQ[C/R]:"+String.format("%.2f/%.2f",getAvConsensusMapq(),
                                getAvRefMapq()));
             message.append('\t');
             message.append(prefix+"NQS_MM_RATE[C/R]:"+String.format("%.2f/%.2f",getNQSConsensusMMRate(),getNQSRefMMRate()));
             message.append('\t');
             message.append(prefix+"NQS_AV_QUAL[C/R]:"+String.format("%.2f/%.2f",getNQSConsensusAvQual(),getNQSRefAvQual()));

             PrimitivePair.Int strand_cons = getConsensusStrandCounts();
             PrimitivePair.Int strand_ref = getRefStrandCounts();
             message.append('\t');
             message.append(prefix+"STRAND_COUNTS[C/C/R/R]:"+strand_cons.first+"/"+strand_cons.second+"/"+strand_ref.first+"/"+strand_ref.second);

             message.append('\t');
             message.append(prefix+"OFFSET_RSTART:"+from_start_median+"/"+from_start_mad);
             message.append('\t');
             message.append(prefix+"OFFSET_REND:"+from_end_median+"/"+from_end_mad);

             return message.toString();

         }

        /**
         * Places alignment statistics into attribute map and returns the map. If attr parameter is null,
         * a new map is allocated, filled and returned. If attr is not null, new attributes are added to that
         * preexisting map, and the same instance of the (updated) map is returned.
         *
         * @param attr
         * @return
         */
        public Map<String,Object> makeStatsAttributes(Map<String,Object> attr) {
             if ( attr == null ) attr = new HashMap<String, Object>();

             VCFIndelAttributes.recordDepth(getConsensusVariantCount(),getAllVariantCount(),getCoverage(),attr);

             VCFIndelAttributes.recordAvMM(getAvConsensusMismatches(),getAvRefMismatches(),attr);

             VCFIndelAttributes.recordAvMapQ(getAvConsensusMapq(),getAvRefMapq(),attr);

             VCFIndelAttributes.recordNQSMMRate(getNQSConsensusMMRate(),getNQSRefMMRate(),attr);

             VCFIndelAttributes.recordNQSAvQ(getNQSConsensusAvQual(),getNQSRefAvQual(),attr);

             VCFIndelAttributes.recordOffsetFromStart(from_start_median,from_start_mad,attr);

             VCFIndelAttributes.recordOffsetFromEnd(from_end_median,from_end_mad,attr);

             PrimitivePair.Int strand_cons = getConsensusStrandCounts();
             PrimitivePair.Int strand_ref = getRefStrandCounts();

             VCFIndelAttributes.recordStrandCounts(strand_cons.first,strand_cons.second,strand_ref.first,strand_ref.second,attr);
             return attr;
         }
    }

    interface IndelListener {
        public void addObservation(int pos, IndelVariant.Type t, String bases, int fromStart, int fromEnd, ExpandedSAMRecord r);
    }

    class WindowContext implements IndelListener {
            private Set<ExpandedSAMRecord> reads;
            private int start=0; // where the window starts on the ref, 1-based
            private CircularArray< List< IndelVariant > > indels;

            private List<IndelVariant> emptyIndelList = new ArrayList<IndelVariant>();


            public WindowContext(int start, int length) {
                this.start = start;
                indels = new CircularArray< List<IndelVariant> >(length);
//                reads = new LinkedList<SAMRecord>();
                reads = new HashSet<ExpandedSAMRecord>();
            }

            /** Returns 1-based reference start position of the interval this object keeps context for.
             *
             * @return
             */
            public int getStart() { return start; }

            /** Returns 1-based reference stop position (inclusive) of the interval this object keeps context for.
             *
             * @return
             */
            public int getStop() { return start + indels.length() - 1; }

            /** Resets reference start position to 0 and clears the context.
             *
             */
            public void clear() {
                start = 0;
                reads.clear();
                indels.clear();
            }

        /**
         * Returns true if any indel observations are present in the specified interval
         * [begin,end] (1-based, inclusive). Interval can be partially of fully outside of the
         * current context window: positions outside of the window will be ignored.
         * @param begin
         * @param end
         */
            public boolean hasIndelsInInterval(long begin, long end) {
                for ( long k = Math.max(start,begin); k < Math.min(getStop(),end); k++ ) {
                    if ( indelsAt(k) != emptyIndelList ) return true;
                }
                return false;               
            }

            public Set<ExpandedSAMRecord> getReads() { return reads; }

            /** Returns the number of reads spanning over the specified reference position                                                                                                       
             * (regardless of whether they have a base or indel at that specific location).
             * The second argument controls whether to count with indels in mind (this is relevant for insertions only,
             * deletions do not require any special treatment since they occupy non-zero length on the ref and since
             * alignment can not start or end with a deletion). For insertions, note that, internally, we assign insertions
             * to the reference position right after the actual event, and we count all events assigned to a given position.
             * This count (reads with indels) should be contrasted to reads without indels, or more rigorously, reads
             * that support the ref rather than the indel. Few special cases may occur here:
             * 1) an alignment that ends (as per getAlignmentEnd()) right before the current position but has I as its
             * last element: we have to count that read into the "coverage" at the current position for the purposes of indel
             * assessment, as the indel in that read <i>will</i> be counted at the current position, so the total coverage
             * should be consistent with that.
             */
             /* NOT IMPLEMENTED: 2) alsignments that start exactly at the current position do <i>not</i> count for the purpose of insertion
             * assessment since they do not contribute any evidence to either Ref or Alt=insertion hypothesis, <i>unless</i>
             * the alignment starts with I (so that we do have evidence for an indel assigned to the current position and
             * read should be counted). For deletions, reads starting at the current position should always be counted (as they
             * show no deletion=ref).
             * @param refPos position on the reference; must be within the bounds of the window
             */
            public int coverageAt(final long refPos, boolean countForIndels) {
                int cov = 0;
                for ( ExpandedSAMRecord read : reads ) {
                    if ( read.getSAMRecord().getAlignmentStart() > refPos || read.getSAMRecord().getAlignmentEnd() < refPos ) {
                        if ( countForIndels && read.getSAMRecord().getAlignmentEnd() == refPos - 1) {
                            Cigar c = read.getSAMRecord().getCigar();
                            if ( c.getCigarElement(c.numCigarElements()-1).getOperator() == CigarOperator.I ) cov++;
                        }
                        continue;
                    }
                    cov++;
                } 
                return cov;
            }


            /** Shifts current window to the right along the reference contig by the specified number of bases.
             * The context will be updated accordingly (indels and reads that go out of scope will be dropped).
             * @param offset
             */
            public void shift(int offset) {
                start += offset;

                indels.shiftData(offset);
                if ( indels.get(0) != null && indels.get(0).size() != 0 ) {
                    IndelVariant indel =  indels.get(0).get(0);

                    System.out.println("WARNING: Indel(s) at first position in the window ("+refName+":"+start+"): currently not supported: "+
                    (indel.getType()==IndelVariant.Type.I?"+":"-")+indel.getBases()+"; read: "+indel.getReadSet().iterator().next().getSAMRecord().getReadName()+"; site ignored");
                    indels.get(0).clear();
//                    throw new StingException("Indel found at the first position ("+start+") after a shift was performed: currently not supported: "+
//                    (indel.getType()==IndelVariant.Type.I?"+":"-")+indel.getBases()+"; reads: "+indel.getReadSet().iterator().next().getSAMRecord().getReadName());
                }
                
                Iterator<ExpandedSAMRecord> read_iter = reads.iterator();

                while ( read_iter.hasNext() ) {
                    ExpandedSAMRecord r = read_iter.next();
                    if ( r.getSAMRecord().getAlignmentEnd() < start ) { // discard reads and associated data that went out of scope
                        read_iter.remove();
                    }
                }
            }

            public void add(SAMRecord read, byte [] ref) {

                if ( read.getAlignmentStart() < start ) return; // silently ignore reads starting before the window start

                ExpandedSAMRecord er = new ExpandedSAMRecord(read,ref,read.getAlignmentStart()-start,this);
                //TODO duplicate records may actually indicate a problem with input bam file; throw an exception when the bug in CleanedReadInjector is fixed
                if ( reads.contains(er)) return; // ignore duplicate records
                reads.add(er);
            }

            public void addObservation(int pos, IndelVariant.Type type, String bases, int fromStart, int fromEnd, ExpandedSAMRecord rec) {
                List<IndelVariant> indelsAtSite;
                try {
                    indelsAtSite = indels.get(pos);
                } catch (IndexOutOfBoundsException e) {
                    SAMRecord r = rec.getSAMRecord();
                    System.out.println("Failed to add indel observation, probably out of coverage window bounds (trailing indel?):\nRead "+
                            r.getReadName()+": "+
                        "read length="+r.getReadLength()+"; cigar="+r.getCigarString()+"; start="+
                        r.getAlignmentStart()+"; end="+r.getAlignmentEnd()+"; window start="+getStart()+
                        "; window end="+getStop());
                    throw e;
                }

                if ( indelsAtSite == null ) {
                    indelsAtSite = new ArrayList<IndelVariant>();
                    indels.set(pos, indelsAtSite);
                }

                IndelVariant indel = null;
                for ( IndelVariant v : indelsAtSite ) {
                    if ( ! v.equals(type, bases) ) continue;

                    indel = v;
                    indel.addObservation(rec);
                    break;
                }
                
                if ( indel == null ) {  // not found:
                    indel = new IndelVariant(rec, type, bases);
                    indelsAtSite.add(indel);
                }
                indel.addReadPositions(fromStart,fromEnd);
            }

            public List<IndelVariant> indelsAt( final long refPos ) {
                List<IndelVariant> l = indels.get((int)( refPos - start ));
                if ( l == null ) return emptyIndelList;
                else return l;
            }


        }


    class ExpandedSAMRecord {
        private SAMRecord read;
        private byte[] mismatch_flags;
        private byte[] expanded_quals;
        private int mms;

        public ExpandedSAMRecord(SAMRecord r, byte [] ref, long offset, IndelListener l) {

            read = r;
            final long rStart = read.getAlignmentStart();
            final long rStop = read.getAlignmentEnd();
            final byte[] readBases = read.getReadString().toUpperCase().getBytes();

            ref = new String(ref).toUpperCase().getBytes();

            mismatch_flags = new byte[(int)(rStop-rStart+1)];
            expanded_quals = new byte[(int)(rStop-rStart+1)];

            // now let's extract indels:

            Cigar c = read.getCigar();
            final int nCigarElems = c.numCigarElements();


            int readLength = 0; // length of the aligned part of the read NOT counting clipped bases
            for ( CigarElement cel : c.getCigarElements() ) {

                switch(cel.getOperator()) {
                case H:
                case S:
                case D:
                case N:
                case P:
                    break; // do not count gaps or clipped bases
                case I:
                case M:
                    readLength += cel.getLength();
                    break; // advance along the gapless block in the alignment
                default :
                    throw new IllegalArgumentException("Unexpected operator in cigar string: "+cel.getOperator());
                }
            }

            int fromStart = 0;
            int posOnRead = 0;
            int posOnRef = 0; // the chunk of reference ref[] that we have access to is aligned with the read:
                                  // its start on the actual full reference contig is r.getAlignmentStart()
            for ( int i = 0 ; i < nCigarElems ; i++ ) {

                final CigarElement ce = c.getCigarElement(i);
                IndelVariant.Type type = null;
                String indel_bases = null;
                int eventPosition = posOnRef;

                switch(ce.getOperator()) {
                case H: break; // hard clipped reads do not have clipped indel_bases in their sequence, so we just ignore the H element...
                case I:
                    type = IndelVariant.Type.I;
                    indel_bases = read.getReadString().substring(posOnRead,posOnRead+ce.getLength());
                    // will increment position on the read below, there's no 'break' statement yet...
                case S:
                    // here we also skip soft-clipped indel_bases on the read; according to SAM format specification,
                    // alignment start position on the reference points to where the actually aligned
                    // (not clipped) indel_bases go, so we do not need to increment reference position here
                    posOnRead += ce.getLength();
                    break;
                case D:
                    type = IndelVariant.Type.D;
                    indel_bases = new String( ref, posOnRef, ce.getLength() );
                    for( int k = 0 ; k < ce.getLength(); k++, posOnRef++ ) mismatch_flags[posOnRef] = expanded_quals[posOnRef] = -1;

                    break;
                case M:
                    for ( int k = 0; k < ce.getLength(); k++, posOnRef++, posOnRead++ ) {
                        if ( readBases[posOnRead] != ref[posOnRef] )  { // mismatch!
                            mms++;
                            mismatch_flags[posOnRef] = 1;
                        }
                        expanded_quals[posOnRef] = read.getBaseQualities()[posOnRead];
                    }
                    fromStart += ce.getLength();
                    break; // advance along the gapless block in the alignment
                default :
                    throw new IllegalArgumentException("Unexpected operator in cigar string: "+ce.getOperator());
                }

                if ( type == null ) continue; // element was not an indel, go grab next element...

                // we got an indel if we are here...
                if ( i == 0 ) logger.debug("Indel at the start of the read "+read.getReadName());
                if ( i == nCigarElems - 1) logger.debug("Indel at the end of the read "+read.getReadName());

                // note that here we will be assigning indels to the first deleted base or to the first
                // base after insertion, not to the last base before the event!
                int fromEnd = readLength - fromStart;
                if ( type == IndelVariant.Type.I ) fromEnd -= ce.getLength();

                l.addObservation((int)(offset+eventPosition), type, indel_bases, fromStart, fromEnd, this);

                if ( type == IndelVariant.Type.I ) fromStart += ce.getLength();

            }
        }

        public SAMRecord getSAMRecord() { return read; }

        public byte [] getExpandedMMFlags() { return mismatch_flags; }

        public byte [] getExpandedQuals() { return expanded_quals; }

        public int getMMCount() { return mms; }

        public boolean equals(Object o) {
            if ( this == o ) return true;
            if ( read == null ) return false;
            if ( o instanceof SAMRecord ) return read.equals(o);
            if ( o instanceof ExpandedSAMRecord ) return read.equals(((ExpandedSAMRecord)o).read);
            return false;
        }


    }

}


class VCFIndelAttributes {
    public static String ALLELIC_DEPTH_KEY = "AD";
    public static String DEPTH_TOTAL_KEY = VCFConstants.DEPTH_KEY;

    public static String MAPQ_KEY = "MQS";

    public static String MM_KEY = "MM";

    public static String NQS_MMRATE_KEY = "NQSMM";

    public static String NQS_AVQ_KEY = "NQSBQ";

    public static String STRAND_COUNT_KEY = "SC";
    public static String RSTART_OFFSET_KEY = "RStart";
    public static String REND_OFFSET_KEY = "REnd";

    public static Set<? extends VCFHeaderLine> getAttributeHeaderLines() {
        Set<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();

        lines.add(new VCFFormatHeaderLine(ALLELIC_DEPTH_KEY, 2, VCFHeaderLineType.Integer, "# of reads supporting consensus indel/reference at the site"));
        lines.add(new VCFFormatHeaderLine(DEPTH_TOTAL_KEY, 1, VCFHeaderLineType.Integer, "Total coverage at the site"));

        lines.add(new VCFFormatHeaderLine(MAPQ_KEY, 2, VCFHeaderLineType.Float, "Average mapping qualities of consensus indel-supporting reads/reference-supporting reads"));

        lines.add(new VCFFormatHeaderLine(MM_KEY, 2, VCFHeaderLineType.Float, "Average # of mismatches per consensus indel-supporting read/per reference-supporting read"));

        lines.add(new VCFFormatHeaderLine(NQS_MMRATE_KEY, 2, VCFHeaderLineType.Float, "Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads"));

        lines.add(new VCFFormatHeaderLine(NQS_AVQ_KEY, 2, VCFHeaderLineType.Float, "Within NQS window: average quality of bases from consensus indel-supporting reads/from reference-supporting reads"));

        lines.add(new VCFFormatHeaderLine(STRAND_COUNT_KEY, 4, VCFHeaderLineType.Integer, "Strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads"));

        lines.add(new VCFFormatHeaderLine(RSTART_OFFSET_KEY, 2, VCFHeaderLineType.Integer, "Median/mad of indel offsets from the starts of the reads"));
        lines.add(new VCFFormatHeaderLine(REND_OFFSET_KEY, 2, VCFHeaderLineType.Integer, "Median/mad of indel offsets from the ends of the reads"));

        return lines;
    }

    public static Map<String,Object> recordStrandCounts(int cnt_cons_fwd, int cnt_cons_rev, int cnt_ref_fwd, int cnt_ref_rev, Map<String,Object> attrs) {
        attrs.put(STRAND_COUNT_KEY, new Integer[] {cnt_cons_fwd, cnt_cons_rev, cnt_ref_fwd, cnt_ref_rev} );
        return attrs;
    }

    public static Map<String,Object> recordDepth(int cnt_cons, int cnt_indel, int cnt_total, Map<String,Object> attrs) {
        attrs.put(ALLELIC_DEPTH_KEY, new Integer[] {cnt_cons, cnt_indel} );
        attrs.put(DEPTH_TOTAL_KEY, cnt_total);
        return attrs;
    }

    public static Map<String,Object> recordAvMapQ(double cons, double ref, Map<String,Object> attrs) {
        attrs.put(MAPQ_KEY, new Float[] {(float)cons, (float)ref} );
        return attrs;
    }

    public static Map<String,Object> recordAvMM(double cons, double ref, Map<String,Object> attrs) {
        attrs.put(MM_KEY, new Float[] {(float)cons, (float)ref} );
        return attrs;
    }

    public static Map<String,Object> recordNQSMMRate(double cons, double ref, Map<String,Object> attrs) {
        attrs.put(NQS_MMRATE_KEY, new Float[] {(float)cons, (float)ref} );
        return attrs;
    }

    public static Map<String,Object> recordNQSAvQ(double cons, double ref, Map<String,Object> attrs) {
        attrs.put(NQS_AVQ_KEY, new Float[] {(float)cons, (float)ref} );
        return attrs;
    }

    public static Map<String,Object> recordOffsetFromStart(int median, int mad, Map<String,Object> attrs) {
        attrs.put(RSTART_OFFSET_KEY, new Integer[] {median, mad} );
        return attrs;
    }

    public static Map<String,Object> recordOffsetFromEnd(int median, int mad, Map<String,Object> attrs) {
        attrs.put(REND_OFFSET_KEY, new Integer[] {median, mad} );
        return attrs;
    }
}
