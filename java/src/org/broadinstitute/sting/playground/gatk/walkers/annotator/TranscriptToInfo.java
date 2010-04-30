package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.AnnotatorROD;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;


/**
 *  Takes a UCSC refGene table ( -B arg must be: -B refgene,AnnotatorInfoTable,/path/to/refgene_file.txt) and generates the big table of nucleotides containing
 *  annotations for each possible variant at each transcript position (eg. 4 variants for each position).
 *
 *  The map & reduce types are both TreeMap<String, String>
 *  TreeMap entries represent lines in the output file. The TreeMap key is a combination of a given output line's position (so that this key can be used to sort all output lines
 *  by reference order), as well as allele and gene name (so that its unique across all output lines). The String value is the output line itself.
 */
@Reference(window=@Window(start=-4,stop=4))
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE}, referenceMetaData={ @RMD(name="refgene",type=AnnotatorROD.class) })
public class TranscriptToInfo extends RodWalker<TreeMap<String, String>, TreeMap<String, String>> implements TreeReducible<TreeMap<String, String>>
{

/**
   RefGene column names:
    # bin           182                     smallint(5) unsigned    range    Indexing field to speed chromosome range queries.
    # name          NM_021649               varchar(255)            values   Name of gene (usually transcript_id from GTF)
    # chrom         chr5                    varchar(255)            values   Reference sequence chromosome or scaffold
    # strand        -                       char(1)                 values   + or - for strand
    # txStart       114942238               int(10) unsigned        range    Transcription start position
    # txEnd         114966008               int(10) unsigned        range    Transcription end position
    # cdsStart      114944144               int(10) unsigned        range    Coding region start
    # cdsEnd        114944852               int(10) unsigned        range    Coding region end
    # exonCount     2                       int(10) unsigned        range    Number of exons
    # exonStarts    114942238,114965692,    longblob                         Exon start positions
    # exonEnds      114944911,114966008,    longblob                         Exon end positions
    # id            0                       int(10) unsigned        range    Unique identifier
    # name2         TICAM2                  varchar(255)            values   Alternate name (e.g. gene_id from GTF)
    # cdsStartStat  cmpl                    enum('none','unk','incmpl','cmpl')    values    enum('none','unk','incmpl','cmpl')
    # cdsEndStat    cmpl                    enum('none','unk','incmpl','cmpl')    values    enum('none','unk','incmpl','cmpl')
    # exonFrames    0,-1,                   longblob                         Exon frame {0,1,2}, or -1 if no frame for exon
*/


    /** Output columns */
    public static final String OUTPUT_CHRPOS = AnnotatorROD.CHRPOS_COLUMN;
    public static final String OUTPUT_HAPLOTYPE_REFERENCE = AnnotatorROD.HAPLOTYPE_REFERENCE_COLUMN;
    public static final String OUTPUT_HAPLOTYPE_ALTERNATE = AnnotatorROD.HAPLOTYPE_ALTERNATE_COLUMN;
    public static final String OUTPUT_HAPLOTYPE_STRAND = AnnotatorROD.HAPLOTYPE_STRAND_COLUMN;

    public static final String OUTPUT_GENE_NAME = "geneName"; //eg. NDUFS2
    public static final String OUTPUT_TRANSCRIPT_ACCESSION = "transcriptId"; //eg. NM_001212414.4

    //public static final String OUTPUT_IN_TRANSCRIPT = "inTranscript"; //eg. true
    public static final String OUTPUT_IN_CODING_REGION = "inCodingRegion"; //eg. true

    public static final String OUTPUT_FRAME = "frame"; //eg. 0,1,2
    public static final String OUTPUT_POSITION_TYPE = "positionType"; //eg. utr5, cds, utr3, intron, intergenic

    public static final String OUTPUT_MRNA_COORD = "mrnaCoord"; //1-based offset within the transcript
    public static final String OUTPUT_CODING_COORD = "codingCoord"; //1-based offset within the cds region

    public static final String OUTPUT_SPLICE_DISTANCE = "spliceDist"; //eg. integer, bp to nearest exon/intron boundary

    public static final String OUTPUT_CODON_NUMBER = "codonCoord"; //eg. 20
    public static final String OUTPUT_REFERENCE_CODON = "referenceCodon";
    public static final String OUTPUT_REFERENCE_AA = "referenceAA";
    public static final String OUTPUT_VARIANT_CODON = "variantCodon";
    public static final String OUTPUT_VARIANT_AA = "variantAA";

    public static final String OUTPUT_CHANGES_AMINO_ACID = "changesAA"; //eg. true
    public static final String OUTPUT_FUNCTIONAL_CLASS = "functionalClass"; //eg. missense

    public static final String OUTPUT_CODING_COORD_STR = "codingCoordStr";
    public static final String OUTPUT_PROTEIN_COORD_STR = "proteinCoordStr";

    public static final String OUTPUT_SPLICE_INFO = "spliceInfo"; //(eg "splice-donor -4", or "splice-acceptor 3") for the 10bp surrounding each exon/intron boundary
    public static final String OUTPUT_UORF_CHANGE = "uorfChange"; // (eg +1 or -1, indicating the addition or interruption of an ATG trinucleotide in the annotated utr5)



    public static final String[] OUTPUT_COLUMN_NAMES = new String[] {

        OUTPUT_CHRPOS,
        OUTPUT_HAPLOTYPE_REFERENCE,
        OUTPUT_HAPLOTYPE_ALTERNATE,
        OUTPUT_HAPLOTYPE_STRAND,

        OUTPUT_GENE_NAME,
        OUTPUT_TRANSCRIPT_ACCESSION,

        OUTPUT_POSITION_TYPE,


        OUTPUT_FRAME,
        OUTPUT_MRNA_COORD,
        OUTPUT_CODING_COORD,
        OUTPUT_CODON_NUMBER,
        OUTPUT_SPLICE_DISTANCE,

        OUTPUT_REFERENCE_CODON,
        OUTPUT_REFERENCE_AA,
        OUTPUT_VARIANT_CODON,
        OUTPUT_VARIANT_AA,
        OUTPUT_CHANGES_AMINO_ACID,
        OUTPUT_FUNCTIONAL_CLASS,

        OUTPUT_CODING_COORD_STR,
        OUTPUT_PROTEIN_COORD_STR,

        OUTPUT_IN_CODING_REGION,
        OUTPUT_SPLICE_INFO,
        OUTPUT_UORF_CHANGE,
    };


    /**
     * Container for all data fields from a single refGene.txt row.
     */
    protected static class RefGeneTranscriptRecord
    {
        public static final String REFGENE_NAME = "name";  //eg. NM_021649
        public static final String REFGENE_STRAND = "strand";  //eg. +
        public static final String REFGENE_TXSTART = "txStart";
        public static final String REFGENE_TXEND = "txEnd";
        public static final String REFGENE_CDS_START = "cdsStart";
        public static final String REFGENE_CDS_END = "cdsEnd";
        public static final String REFGENE_EXON_COUNT = "exonCount";
        public static final String REFGENE_EXON_STARTS = "exonStarts";
        public static final String REFGENE_EXON_ENDS = "exonEnds";
        public static final String REFGENE_EXON_FRAMES = "exonFrames";
        public static final String REFGENE_NAME2 = "name2"; //eg. TICAM2


        /**
         * This StringBuffer accumulates the entire transcript sequence.
         * This buffer is used instead of using the GATK window mechanism
         * because arbitrary-length look-aheads and look-behinds are needed to deal
         * with codons that span splice-junctions. The window mechanism
         * requires hard-coding the window size, which would translate into a
         * limit on maximum supported intron size. To avoid this, the
         * sequence is accumulated as the transcript is scanned left-to-right.
         * Then, all calculations are performed at the end.
         */
        public StringBuilder txSequence;  //the sequence of the entire transcript in order from 5' to 3'
        public StringBuilder utr5Sequence; //the protein coding sequence (with introns removed) in order from 5' to 3'
        public StringBuilder cdsSequence; //the protein coding sequence (with introns removed) in order from 5' to 3'

        public boolean positiveStrand;
        public String name; //eg. NM_021649
        public String name2; //eg. TICAM2

        public String txChrom; //The chromosome
        public long txStart;
        public long txEnd;

        public long cdsStart;
        public long cdsEnd;

        public long[] exonStarts;
        public long[] exonEnds;
        public int[] exonFrames;

        private AnnotatorROD rod;


        public RefGeneTranscriptRecord(final AnnotatorROD refGeneRod) {
            this.rod = refGeneRod;
            //String binStr = rod.get("bin");
            //String idStr = refGeneRod.get("id"); //int(10) unsigned range Unique identifier ( usually 0 for some reason - even for translated )
            String strandStr = refGeneRod.get("strand");
            if(strandStr == null) {
                throw new IllegalArgumentException("refGene record doesn't contain a 'strand' column. Make sure refGene input file has a header and the usual columns: \"" + strandStr + "\"");
            } else if(strandStr.equals("+")) {
                positiveStrand = true;
            } else if(strandStr.equals("-")) {
                positiveStrand = false;
            } else {
                throw new IllegalArgumentException("refGene record contains unexpected value for 'strand' column: \"" + strandStr + "\"");
            }

            name = refGeneRod.get("name");
            name2 = refGeneRod.get("name2");

            //String txStartStr = refGeneRod.get(REFGENE_TXSTART);  //These fields are used to generate column 1 of the ROD file (eg. they get turned into chr:txStart-txStop)
            //String txEndStr = refGeneRod.get(REFGENE_TXEND);

            final GenomeLoc loc = refGeneRod.getLocation();
            txChrom = loc.getContig();
            txStart = loc.getStart();
            txEnd = loc.getStop();

            String cdsStartStr = refGeneRod.get(REFGENE_CDS_START);
            String cdsEndStr = refGeneRod.get(REFGENE_CDS_END);

            cdsStart = Long.parseLong(cdsStartStr);
            cdsEnd = Long.parseLong(cdsEndStr);

            txSequence = new StringBuilder( (int) (txEnd - txStart + 1) );  //the sequence of the entire transcript in order from 5' to 3'
            if(isProteinCodingTranscript()) {
                utr5Sequence = new StringBuilder( positiveStrand ? (int) (cdsStart - txStart + 1) : (int) (txEnd - cdsEnd + 1) ); //TODO reduce init size by size of introns
                cdsSequence = new StringBuilder( (int) (cdsEnd - cdsStart + 1) ); //TODO reduce init size by size of introns
            }

            String exonCountStr = refGeneRod.get(REFGENE_EXON_COUNT);
            String exonStartsStr = refGeneRod.get(REFGENE_EXON_STARTS);
            String exonEndsStr = refGeneRod.get(REFGENE_EXON_ENDS);
            String exonFramesStr = refGeneRod.get(REFGENE_EXON_FRAMES);

            String[] exonStartStrs = exonStartsStr.split(",");
            String[] exonEndStrs = exonEndsStr.split(",");
            String[] exonFrameStrs = exonFramesStr.split(",");

            int exonCount = Integer.parseInt(exonCountStr);
            if(exonCount != exonStartStrs.length || exonCount != exonEndStrs.length || exonCount != exonFrameStrs.length)
            {
                throw new RuntimeException("exonCount != exonStarts.length || exonCount != exonEnds.length || exonCount != exonFrames.length. Exon starts: " + exonStartsStr + ", Exon ends: " + exonEndsStr + ",  Exon frames: " + exonFramesStr + ", Exon count: " + exonCountStr +". refGeneRod = " + refGeneRod);
            }

            exonStarts = new long[exonCount];
            exonEnds = new long[exonCount];
            exonFrames = new int[exonCount];
            for(int i = 0; i < exonCount; i++) {
                exonStarts[i] = Long.parseLong(exonStartStrs[i]);
                  exonEnds[i] = Long.parseLong(exonEndStrs[i]);
                exonFrames[i] = Integer.parseInt(exonFrameStrs[i]);
            }
        }


        /**
         * Takes a genomic position on the same contig as the transcript, and
         * returns true if this position falls within an exon.
         */
        public boolean isWithinExon(final long genomPosition) {
            for(int i = 0; i < exonStarts.length; i++) {
                final long curStart = exonStarts[i];
                if(genomPosition < curStart) {
                    return false;
                }
                final long curStop = exonEnds[i];
                if(genomPosition <= curStop) {
                    return true;
                }
            }

            return false;
        }

        /**
         * Computes the distance to the nearest splice-site.
         * The returned value is negative its on the 5' side (eg. upstream) of the juntion, and
         * positive if its on the 3' side.
         */
        public int computeDistanceToNearestSpliceSite(final long genomPosition) {
            int prevDistance = Integer.MAX_VALUE;
            for(int i = 0; i < exonStarts.length; i++) {
                final long curStart = exonStarts[i];
                int curDistance = (int) (curStart-genomPosition);
                if(genomPosition < curStart) {
                    //position is within the current intron
                    if(prevDistance < curDistance) {
                        return positiveStrand ? prevDistance : -prevDistance;
                    } else {
                        return positiveStrand ? -curDistance : curDistance;
                    }
                } else {
                    prevDistance = (int) (genomPosition - curStart) + 1;
                }

                final long curStop = exonEnds[i];
                curDistance = (int) (curStop-genomPosition) + 1;
                if(genomPosition <= curStop) {
                    //position is within an exon
                    if(prevDistance < curDistance) {
                        return positiveStrand ? prevDistance : -prevDistance;
                    } else {
                        return positiveStrand ? -curDistance : curDistance;
                    }
                } else {
                    prevDistance = (int) (genomPosition - curStop);
                }
            }

            throw new IllegalArgumentException("Genomic position: [" + genomPosition +"] not found within transcript: " + this +". " +
                    "This method should not have been called for this position. NOTE: this method assumes that all transcripts start " +
                    "with an exon and end with an exon (rather than an intron). Is this wrong?");
            //return prevDistance; //out of exons. return genomPosition-curStop
        }


        /**
         * Returns true if this is a coding transcript (eg. is translated
         * into proteins). Returns false for non-coding RNA.
         */
        public boolean isProteinCodingTranscript() {
            return cdsStart < cdsEnd;
        }

        @Override
        public String toString() {
            return rod.toString();
        }



        /**
         * Computes the coding coord of the 1st nucleotide in the transcript.
         * If the 1st nucleotide is in the 5'utr, the returned value will be negative.
         * Otherwise (if the 1st nucleotide is CDS), the returned value is 1.
         */
        public int computeInitialCodingCoord() {
            if(!isProteinCodingTranscript()) {
                throw new StingException("This method should only be called for protein-coding transcripts");
            }

            if(positiveStrand)
            {
                if( cdsStart == exonStarts[0] ) {
                    //the 1st nucleotide of the transcript is CDS.
                    return 1;
                }

                int result = 0;
                for(int i = 0; i < exonStarts.length; i++)
                {
                    final long exonStart = exonStarts[i];
                    final long exonEnd = exonEnds[i];
                    if(cdsStart <= exonEnd) { //eg. exonEnd is now on the 3' side of cdsStart
                        //this means cdsStart is within the current exon
                        result += (cdsStart - exonStart) + 1;
                        break;
                    } else {
                        //cdsStart is downstream of the current exon
                        result += (exonEnd - exonStart) + 1;
                    }
                }
                return -result; //negate because 5' UTR coding coord is negative
            }
            else //(negative strand)
            {
                final long cdsStart_5prime = cdsEnd;
                if(cdsStart_5prime == exonEnds[exonEnds.length - 1]) {
                    //the 1st nucleotide of the transcript is CDS.
                    return 1;
                }

                int result = 0;
                for(int i = exonEnds.length - 1; i >= 0; i--)
                {
                    final long exonStart = exonEnds[i]; //when its the negative strand, the 5' coord of the 1st exon is exonEnds[i]
                    final long exonEnd = exonStarts[i];
                    if( exonEnd <= cdsStart_5prime ) { //eg. exonEnd is now on the 3' side of cdsStart
                        //this means cdsStart is within the current exon
                        result += -(cdsStart_5prime - exonStart) + 1;
                        break;
                    } else {
                        //cdsStart is downstream of the current exon
                        result += -(exonEnd - exonStart) + 1;
                    }
                }
                return -result; //negate because 5' UTR coding coord is negative
            }
        }
    }


    protected String outputFilename = null;

    protected int transcriptsProcessedCounter = 0;


    /** Possible values for the "POSITION_TYPE" output column. */
    enum PositionType {
        intergenic, intron, utr5, CDS, utr3, non_coding_exon, non_coding_intron
    }


    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        //init the output file
        outputFilename = getToolkit().getArguments().outFileName;
        if(outputFilename == null) {
            throw new StingException("Output file not specified. (Used -o command-line-arg)" );
        }

        if(!new File(outputFilename).canWrite()) {
            throw new StingException("Unable to write to output file: " + outputFilename );
        }

    }

    /**
     * Initialize the output filename.
     *
     * @return 0
     */
    public TreeMap<String, String> reduceInit()
    {
        return new TreeMap<String, String>() {

            public String toString() {
                return size() + " total entries";
            }
        };
    }


    /**
     * For each site of interest, generate the appropriate fields.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public TreeMap<String, String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( tracker == null )
            return null;


        final Collection<RODRecordList> rods = tracker.getBoundRodTracks();
        if(rods.size() == 0) {
            //if there's nothing overlapping this locus, skip it.
            return null;
        }

        final TreeMap<String, String> result = new TreeMap<String, String>();


        final List<Object> refGeneRODs = (List<Object>) tracker.getReferenceMetaData("refgene");

        for(Object refGeneRodObject : refGeneRODs) {
            if(! (refGeneRodObject instanceof AnnotatorROD) ) {
                throw new RuntimeException("The 'refgene' -B arg must have the type: 'AnnotatorInputTable'.");
            }

            final AnnotatorROD refGeneRod = (AnnotatorROD) refGeneRodObject;
            RefGeneTranscriptRecord parsedRefGeneRod = (RefGeneTranscriptRecord) refGeneRod.getTemporaryAttribute("parsedRefGeneRod");
            if( parsedRefGeneRod == null ) {
                parsedRefGeneRod = new RefGeneTranscriptRecord(refGeneRod);
                refGeneRod.setTemporaryAttribute("parsedRefGeneRod", parsedRefGeneRod);
            }

            //populate the parsedRefGeneRod.positiveStrand data structure
            if(parsedRefGeneRod.positiveStrand) {
                parsedRefGeneRod.txSequence.append(ref.getBase());
            } else {
                parsedRefGeneRod.txSequence.insert(0, ref.getBase());
            }

            //continue populating the parsedRefGeneRod.utr5Sequence and  parsedRefGeneRod.cdsSequence data structures
            final long position = ref.getLocus().getStart();
            if(parsedRefGeneRod.isProteinCodingTranscript() && parsedRefGeneRod.isWithinExon(position) )
            {
                //we're within an exon of a proteinCodingTranscript

                if(parsedRefGeneRod.positiveStrand)
                {
                    if(position < parsedRefGeneRod.cdsStart)
                    {
                        parsedRefGeneRod.utr5Sequence.append(ref.getBase()); //within utr5
                    }
                    else if(position >= parsedRefGeneRod.cdsStart && position <= parsedRefGeneRod.cdsEnd)
                    {
                        parsedRefGeneRod.cdsSequence.append(ref.getBase()); //within CDS
                    }
                }
                else
                {
                    if(position > parsedRefGeneRod.cdsEnd)
                    {
                        //As we more left-to-right (aka. 3' to 5'), we do insert(0,..) to reverse the sequence so that it become 5' to 3' in parsedRefGeneRod.utr5Sequence.
                        parsedRefGeneRod.utr5Sequence.insert(0,ref.getBase()); //within utr5.
                    }
                    else if(position >= parsedRefGeneRod.cdsStart && position <= parsedRefGeneRod.cdsEnd)
                    {
                        parsedRefGeneRod.cdsSequence.insert(0,ref.getBase()); //within CDS
                    }
                }
            }

            if(position == parsedRefGeneRod.txEnd)
            {
                //we've reached the end of the transcript - compute all data and write it out.
                try
                {
                    generateOutputRecordsForROD(parsedRefGeneRod, result);
                }
                catch(Exception e)
                {
                    logger.error(Thread.currentThread().getName() + " - Unexpected error occurred at position: [" + parsedRefGeneRod.txChrom + ":" + position + "] in transcript: " + parsedRefGeneRod, e);
                }

                transcriptsProcessedCounter++;
                if(transcriptsProcessedCounter % 100 == 0) {
                    logger.info(new Date() + ": " +  transcriptsProcessedCounter + " transcripts processed");
                }
            }
        }

        return result;
    }


    private static final char[] ALLELES = {'A','C','G','T'};

    private static final String OUTPUT_HEADER_LINE;

    static {
        StringBuilder outputHeaderLine = new StringBuilder();
        for( final String column : OUTPUT_COLUMN_NAMES ) {
            if(outputHeaderLine.length() != 0) {
                outputHeaderLine.append( AnnotatorROD.DELIMITER );
            }
            outputHeaderLine.append(column);
        }

        OUTPUT_HEADER_LINE = outputHeaderLine.toString();
    }


    /**
     * Returns the filename.
     *
     * @param parsedRefGeneRod
     * @param result
     * @return
     */
    private void generateOutputRecordsForROD(RefGeneTranscriptRecord parsedRefGeneRod, TreeMap<String, String> result) throws IOException
    {
        //Transcripts that don't produce proteins are indicated in refGene by cdsStart == cdsEnd
        //These will be handled by generating only one record, with haplotypeAlternate == "*".
        final boolean isProteinCodingTranscript = parsedRefGeneRod.isProteinCodingTranscript();

        final boolean positiveStrand = parsedRefGeneRod.positiveStrand; //alias

        if(isProteinCodingTranscript && parsedRefGeneRod.cdsSequence.length() % 3 != 0) {
            logger.warn("WARNING: Transcript " + parsedRefGeneRod.name +" at position ["+ parsedRefGeneRod.txChrom + ":" +parsedRefGeneRod.txStart + "-" + parsedRefGeneRod.txEnd + "] has " + parsedRefGeneRod.cdsSequence.length() + " nucleotides in its CDS region, which is not divisible by 3. Skipping...");
            //discard transcripts where CDS length is not a multiple of 3
            return;
        }

        final long txStart_5prime = positiveStrand ? parsedRefGeneRod.txStart : parsedRefGeneRod.txEnd; //1-based, inclusive
        final long txEnd_3prime = positiveStrand ? parsedRefGeneRod.txEnd : parsedRefGeneRod.txStart; //1-based, inclusive
        final int increment_5to3 = positiveStrand ? 1 : -1; //whether to increment or decrement
        final int strandSign = increment_5to3; //alias

        final long cdsStart_5prime = positiveStrand ? parsedRefGeneRod.cdsStart : parsedRefGeneRod.cdsEnd; //1-based, inclusive
        final long cdsEnd_3prime = positiveStrand ? parsedRefGeneRod.cdsEnd : parsedRefGeneRod.cdsStart ; //1-based, inclusive

        int frame = 0; //the frame of the current position
        int txOffset_from5 = 1; //goes from txStart 5' to txEnd 3' for both + and - strand
        int utr5Count_from5 = 0;
        char[] utr5NucBuffer_5to3 = null; //used to find uORFs - size = 5 because to hold the 3 codons that overlap any given position: [-2,-1,0], [-1,0,1], and [0,1,2]

        int codonCount_from5 = 1; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand - counts the number of codons
        int codingCoord_from5 = isProteinCodingTranscript ? parsedRefGeneRod.computeInitialCodingCoord() : -1; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand
        boolean codingCoordResetForCDS = false;
        boolean codingCoordResetForUtr3 = false;
        final char[] currentCodon_5to3 = isProteinCodingTranscript ? new char[3] : null; //holds the current RNA codon - 5' to 3'

        final boolean mitochondrial = parsedRefGeneRod.txChrom.toLowerCase().startsWith("chrm");

        PositionType positionType = null;
        boolean isWithinIntronAndFarFromSpliceJunction = false;
        long intronStart_5prime = -1;
        long intronEnd_5prime = -1;

        final Map<String, String> outputLineFields = new HashMap<String, String>();

        for(long txCoord_5to3 = txStart_5prime; txCoord_5to3 != txEnd_3prime + increment_5to3; txCoord_5to3 += increment_5to3)
        {

            //compute certain attributes of the current position
            final boolean isWithinExon = parsedRefGeneRod.isWithinExon(txCoord_5to3); //TODO if necessary, this can be sped up by keeping track of current exon/intron

            final int distanceToNearestSpliceSite = parsedRefGeneRod.computeDistanceToNearestSpliceSite(txCoord_5to3);
            final boolean isWithin10bpOfSpliceJunction = Math.abs(distanceToNearestSpliceSite) <= 10;


            //increment coding coord is necessary
            if(isWithinExon) {
                codingCoord_from5++;
            }

            //figure out the current positionType
            final PositionType prevPositionType = positionType; //save the position before it is updated
            if(isProteinCodingTranscript)
            {
                if(isWithinExon)
                {
                    if( strandSign*(txCoord_5to3 - cdsStart_5prime) < 0 ) {  //utr5  (multiplying by strandSign is like doing absolute value.)
                        positionType = PositionType.utr5;
                    } else if( strandSign*(txCoord_5to3 - cdsEnd_3prime) > 0 ) {  //utr3  (multiplying by strandSign is like doing absolute value.)
                        positionType = PositionType.utr3;
                    } else {
                        positionType = PositionType.CDS;
                    }
                } else {
                    positionType = PositionType.intron;
                }
            } else {
                if(isWithinExon) {
                    positionType = PositionType.non_coding_exon;
                } else {
                    positionType = PositionType.non_coding_intron;
                }
            }

            //handle transitions
            if(positionType == PositionType.CDS && prevPositionType != PositionType.CDS && !codingCoordResetForCDS) {
                //transitioning from utr5 to CDS, reset the coding coord from -1 to 1.
                codingCoord_from5 = 1;
                codingCoordResetForCDS = true;
            } else if(positionType == PositionType.utr3 && prevPositionType != PositionType.utr3 && !codingCoordResetForUtr3) {
                //transitioning from CDS to utr3, reset the coding coord to 1.
                codingCoord_from5 = 1;
                codingCoordResetForUtr3 = true;
            }


            try
            {

                //handle introns
                boolean wasWithinIntronAndFarFromSpliceJunction = isWithinIntronAndFarFromSpliceJunction;
                isWithinIntronAndFarFromSpliceJunction = !isWithinExon && !isWithin10bpOfSpliceJunction;

                if(!wasWithinIntronAndFarFromSpliceJunction && isWithinIntronAndFarFromSpliceJunction) {
                    //save intron start
                    intronStart_5prime = txCoord_5to3;

                } else if(wasWithinIntronAndFarFromSpliceJunction && !isWithinIntronAndFarFromSpliceJunction) {
                    //output intron record
                    intronEnd_5prime = txCoord_5to3 - increment_5to3;

                    final long intronStart = (intronStart_5prime < intronEnd_5prime ? intronStart_5prime : intronEnd_5prime) ;
                    final long intronEnd = (intronEnd_5prime > intronStart_5prime ? intronEnd_5prime : intronStart_5prime);
                    outputLineFields.clear();
                    outputLineFields.put(OUTPUT_CHRPOS, parsedRefGeneRod.txChrom + ":" + intronStart +  "-" + intronEnd);
                    outputLineFields.put(OUTPUT_HAPLOTYPE_REFERENCE, Character.toString( '*' ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_ALTERNATE, Character.toString( '*' ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_STRAND, positiveStrand ? "+" : "-");
                    outputLineFields.put(OUTPUT_GENE_NAME, parsedRefGeneRod.name2 );
                    outputLineFields.put(OUTPUT_TRANSCRIPT_ACCESSION, parsedRefGeneRod.name );

                    outputLineFields.put(OUTPUT_POSITION_TYPE, positionType.toString() );

                    if(isProteinCodingTranscript)
                    {
                        outputLineFields.put(OUTPUT_IN_CODING_REGION, Boolean.toString(positionType == PositionType.CDS) );
                    }

                    final String sortKey = computeSortKey(parsedRefGeneRod.txChrom, intronStart, '*', parsedRefGeneRod);

                    addThisLineToResult(sortKey, outputLineFields, result);
                }


                //when in utr5, compute the utr5NucBuffer_5to3 which is later used to compute the OUTPUT_UORF_CHANGE field
                if(positionType == PositionType.utr5)
                {
                    if(utr5Count_from5 < parsedRefGeneRod.utr5Sequence.length())
                    {
                        if(utr5NucBuffer_5to3 == null) {
                            //initialize
                            utr5NucBuffer_5to3 = new char[5];
                            utr5NucBuffer_5to3[3] = parsedRefGeneRod.utr5Sequence.charAt( utr5Count_from5 );
                            if(!positiveStrand) {
                                utr5NucBuffer_5to3[3] = BaseUtils.simpleComplement(utr5NucBuffer_5to3[3]);
                            }

                            if(utr5Count_from5 + 1 < parsedRefGeneRod.utr5Sequence.length() ) {
                                utr5NucBuffer_5to3[4] = parsedRefGeneRod.utr5Sequence.charAt( utr5Count_from5 + 1 );
                                if(!positiveStrand) {
                                    utr5NucBuffer_5to3[4] = BaseUtils.simpleComplement(utr5NucBuffer_5to3[4]);
                                }
                            }
                        }

                        //as we move 5' to 3', shift nucleotides down to the 5' end, making room for the new 3' nucleotide:
                        utr5NucBuffer_5to3[0] = utr5NucBuffer_5to3[1];
                        utr5NucBuffer_5to3[1] = utr5NucBuffer_5to3[2];
                        utr5NucBuffer_5to3[2] = utr5NucBuffer_5to3[3];
                        utr5NucBuffer_5to3[3] = utr5NucBuffer_5to3[4];

                        char nextRefBase = 0;
                        if( utr5Count_from5 + 2 < parsedRefGeneRod.utr5Sequence.length() )
                        {
                            nextRefBase = parsedRefGeneRod.utr5Sequence.charAt( utr5Count_from5 + 2 );

                            //reverse complement if necessary
                            if(!positiveStrand) {
                                nextRefBase = BaseUtils.simpleComplement(nextRefBase);
                            }
                        }
                        utr5NucBuffer_5to3[4] = nextRefBase;

                        //check for bad bases
                        if( (utr5NucBuffer_5to3[0] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[0])) ||
                            (utr5NucBuffer_5to3[1] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[1])) ||
                            (utr5NucBuffer_5to3[2] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[2])) ||
                            (utr5NucBuffer_5to3[3] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[3])) ||
                            (utr5NucBuffer_5to3[4] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[4])))
                        {
                            logger.error("Exception: Skipping current position [" + parsedRefGeneRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedRefGeneRod.name +". utr5NucBuffer_5to3 contains irregular base:" + utr5NucBuffer_5to3[0] + utr5NucBuffer_5to3[1] + utr5NucBuffer_5to3[2] + utr5NucBuffer_5to3[3] + utr5NucBuffer_5to3[4]);// +". Transcript is: " + parsedRefGeneRod);
                            continue;
                        }

                    } else { // if(utr5Count_from5 >= parsedRefGeneRod.utr5Sequence.length())
                        //defensive programming
                        logger.error("Exception: Skipping current position [" + parsedRefGeneRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedRefGeneRod.name +". utr5Count_from5 is now " + utr5Count_from5 + ", while parsedRefGeneRod.utr5Sequence.length() == " + parsedRefGeneRod.utr5Sequence.length() + ". This means parsedRefGeneRod.utr5Sequence isn't as long as it should be. This is a bug in handling this record: " + parsedRefGeneRod);
                        continue;
                    }
                }


                //when in CDS, compute current codon
                if(positionType == PositionType.CDS)
                {
                    if(frame == 0)
                    {
                        currentCodon_5to3[0] = parsedRefGeneRod.cdsSequence.charAt( codingCoord_from5 - 1 ); //subtract 1 to go to zero-based coords
                        currentCodon_5to3[1] = parsedRefGeneRod.cdsSequence.charAt( codingCoord_from5 );
                        currentCodon_5to3[2] = parsedRefGeneRod.cdsSequence.charAt( codingCoord_from5 + 1);

                        if(!positiveStrand) {
                            currentCodon_5to3[0] = BaseUtils.simpleComplement(currentCodon_5to3[0]);
                            currentCodon_5to3[1] = BaseUtils.simpleComplement(currentCodon_5to3[1]);
                            currentCodon_5to3[2] = BaseUtils.simpleComplement(currentCodon_5to3[2]);
                        }
                    }

                    //check for bad bases
                    if(!BaseUtils.isRegularBase(currentCodon_5to3[0]) || !BaseUtils.isRegularBase(currentCodon_5to3[1]) || !BaseUtils.isRegularBase(currentCodon_5to3[2])) {
                        logger.error("Exception: Skipping current position [" + parsedRefGeneRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedRefGeneRod.name +". CDS codon contains irregular base:" + currentCodon_5to3[0] + currentCodon_5to3[1] + currentCodon_5to3[2]);// +". Transcript is: " + parsedRefGeneRod);
                        continue;
                    }

                }

                char haplotypeReference = parsedRefGeneRod.txSequence.charAt( txOffset_from5 - 1 );
                if(!BaseUtils.isRegularBase(haplotypeReference)) {
                    //check for bad bases
                    logger.error("Exception: Current position [" + parsedRefGeneRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedRefGeneRod.name + " contains an irregular base:" + haplotypeReference); // +". Transcript is: " + parsedRefGeneRod);
                    continue;
                }

                for(char haplotypeAlternate : ALLELES )
                {
                    outputLineFields.clear();

                    if(!isProteinCodingTranscript || isWithinIntronAndFarFromSpliceJunction) {
                        haplotypeReference = '*';
                        haplotypeAlternate = '*';
                    }

                    //compute simple OUTPUT fields.
                    outputLineFields.put(OUTPUT_CHRPOS, parsedRefGeneRod.txChrom + ":" + txCoord_5to3);
                    outputLineFields.put(OUTPUT_HAPLOTYPE_REFERENCE, Character.toString( haplotypeReference ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_ALTERNATE, Character.toString( haplotypeAlternate ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_STRAND, positiveStrand ? "+" : "-");
                    outputLineFields.put(OUTPUT_GENE_NAME, parsedRefGeneRod.name2 );
                    outputLineFields.put(OUTPUT_TRANSCRIPT_ACCESSION, parsedRefGeneRod.name );

                    outputLineFields.put(OUTPUT_POSITION_TYPE, positionType.toString() );
                    outputLineFields.put(OUTPUT_MRNA_COORD, Integer.toString(txOffset_from5) );
                    outputLineFields.put(OUTPUT_SPLICE_DISTANCE, Integer.toString(distanceToNearestSpliceSite) );

                    //compute OUTPUT_SPLICE_INFO
                    final String spliceInfoString;
                    if(isWithin10bpOfSpliceJunction) {
                        if(distanceToNearestSpliceSite < 0) {
                            //is on the 5' side of the splice junction
                            if(isWithinExon) {
                                spliceInfoString = "splice-donor " + distanceToNearestSpliceSite;
                            } else {
                                spliceInfoString = "splice-acceptor " + distanceToNearestSpliceSite;
                            }
                        } else {
                            if(isWithinExon) {
                                spliceInfoString = "splice-acceptor " + distanceToNearestSpliceSite;
                            } else {
                                spliceInfoString = "splice-donor " + distanceToNearestSpliceSite;
                            }
                        }
                        outputLineFields.put(OUTPUT_SPLICE_INFO, spliceInfoString);
                    }

                    //compute OUTPUT_IN_CODING_REGION
                    if(isProteinCodingTranscript)
                    {
                        outputLineFields.put(OUTPUT_IN_CODING_REGION, Boolean.toString(positionType == PositionType.CDS) );
                    }


                    //compute OUTPUT_UORF_CHANGE
                    if(positionType == PositionType.utr5)
                    {
                        String refCodon1 = (Character.toString(utr5NucBuffer_5to3[0]) + Character.toString(utr5NucBuffer_5to3[1]) + utr5NucBuffer_5to3[2]).toLowerCase();
                        String refCodon2 = (Character.toString(utr5NucBuffer_5to3[1]) + Character.toString(utr5NucBuffer_5to3[2]) + utr5NucBuffer_5to3[3]).toLowerCase();
                        String refCodon3 = (Character.toString(utr5NucBuffer_5to3[2]) + Character.toString(utr5NucBuffer_5to3[3]) + utr5NucBuffer_5to3[4]).toLowerCase();

                        String varCodon1 = (Character.toString(utr5NucBuffer_5to3[0]) + Character.toString(utr5NucBuffer_5to3[1]) + haplotypeAlternate).toLowerCase();
                        String varCodon2 = (Character.toString(utr5NucBuffer_5to3[1]) + Character.toString(haplotypeAlternate) + utr5NucBuffer_5to3[3]).toLowerCase();
                        String varCodon3 = (Character.toString(haplotypeAlternate) + Character.toString(utr5NucBuffer_5to3[3]) + utr5NucBuffer_5to3[4]).toLowerCase();

                        //check for +1 (eg. addition of new ATG uORF) and -1 (eg. disruption of existing ATG uORF)
                        String uORFChangeStr = null;
                        if( (refCodon1.equals("atg") && !varCodon1.equals("atg")) ||
                            (refCodon2.equals("atg") && !varCodon2.equals("atg")) ||
                            (refCodon3.equals("atg") && !varCodon3.equals("atg")))
                        {
                            uORFChangeStr = "-1";
                        }
                        else if((varCodon1.equals("atg") && !refCodon1.equals("atg")) ||
                                (varCodon2.equals("atg") && !refCodon2.equals("atg")) ||
                                (varCodon3.equals("atg") && !refCodon3.equals("atg")))
                        {
                            uORFChangeStr = "+1";
                        }

                        outputLineFields.put(OUTPUT_UORF_CHANGE, uORFChangeStr );
                    }
                    //compute CDS-specific fields
                    else if(positionType == PositionType.CDS)
                    {
                        final String referenceCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];
                        outputLineFields.put(OUTPUT_FRAME, Integer.toString(frame) );
                        outputLineFields.put(OUTPUT_CODON_NUMBER, Integer.toString(codonCount_from5) );

                        final AminoAcid refAA = AminoAcidTable.getAA(referenceCodon, mitochondrial);
                        outputLineFields.put(OUTPUT_REFERENCE_CODON, referenceCodon );
                        outputLineFields.put(OUTPUT_REFERENCE_AA, refAA.getCode());

                        final char temp = currentCodon_5to3[frame];
                        currentCodon_5to3[frame] = haplotypeAlternate;
                        final String variantCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];
                        final AminoAcid variantAA = AminoAcidTable.getAA(variantCodon, mitochondrial);
                        outputLineFields.put(OUTPUT_VARIANT_CODON, variantCodon );
                        outputLineFields.put(OUTPUT_VARIANT_AA, variantAA.getCode());

                        outputLineFields.put(OUTPUT_PROTEIN_COORD_STR, "p."+ refAA.getLetter() + Integer.toString(codonCount_from5) + variantAA.getLetter() ); //for example: "p.K76A"

                        currentCodon_5to3[frame] = temp;

                        boolean changesAA = !refAA.equals(variantAA);
                        outputLineFields.put(OUTPUT_CHANGES_AMINO_ACID, Boolean.toString(changesAA));
                        final String functionalClass;
                        if(changesAA) {
                            if(variantAA.isStop()) {
                                functionalClass ="nonsense";
                            } else {
                                functionalClass ="missense";
                            }
                        } else {
                            functionalClass ="silent";
                        }
                        outputLineFields.put(OUTPUT_FUNCTIONAL_CLASS, functionalClass);
                    }


                    //compute OUTPUT_CODING_COORD_STR
                    if(isProteinCodingTranscript)
                    {
                        //compute coding coord
                        final StringBuilder codingCoordStr = new StringBuilder();
                        codingCoordStr.append( mitochondrial ? "m." : "c." );
                        if(positionType == PositionType.utr3) {
                            codingCoordStr.append( '*' );
                        }

                        if(isWithinExon) {
                            codingCoordStr.append( Integer.toString(codingCoord_from5) );
                            codingCoordStr.append( haplotypeReference + ">" + haplotypeAlternate);
                        } else {
                            //intronic coordinates
                            if(distanceToNearestSpliceSite < 0) {
                                codingCoordStr.append( Integer.toString(codingCoord_from5 + 1) );
                            } else {
                                codingCoordStr.append( Integer.toString(codingCoord_from5 ) );
                                codingCoordStr.append( "+" );
                            }

                            codingCoordStr.append( Integer.toString( distanceToNearestSpliceSite ) );
                        }

                        outputLineFields.put(OUTPUT_CODING_COORD_STR, codingCoordStr.toString());
                    }


                    //generate the output line and add it to 'result' map.
                    if(!isWithinIntronAndFarFromSpliceJunction) {
                        final String sortKey = computeSortKey(parsedRefGeneRod.txChrom, txCoord_5to3, haplotypeAlternate, parsedRefGeneRod);

                        addThisLineToResult(sortKey, outputLineFields, result);
                    }

                    if( haplotypeAlternate == '*' ) {
                        //need only one record for this position with "*" for haplotypeAlternate, instead of the 4 individual alleles
                        break;
                    }

                } //ALLELE for-loop
            }
            finally
            {
                //increment coords
                txOffset_from5++;

                if(positionType == PositionType.utr5) {
                    utr5Count_from5++;
                } else if(positionType == PositionType.CDS) {
                    frame = (frame + 1) % 3;
                    if(frame == 0) {
                        codonCount_from5++;
                    }
                }
            }
        } // l for-loop

    } //method close


    /**
     * Utility method. Creates a line containing the outputLineFields, and adds it to result, hashed by the sortKey.
     *
     * @param sortKey  Reference-order sort key
     * @param outputLineFields Column-name to value pairs.
     * @param result  The hashmap that accumulates all results.
     */
    private void addThisLineToResult(final String sortKey, final Map<String, String> outputLineFields, TreeMap<String, String> result) {
        final StringBuilder outputLine = new StringBuilder();
        for( final String column : OUTPUT_COLUMN_NAMES ) {
            if(outputLine.length() != 0) {
                outputLine.append( AnnotatorROD.DELIMITER );
            }
            final String value = outputLineFields.get(column);
            if(value != null) {
                outputLine.append(value);
            }
        }

        result.put(sortKey, outputLine.toString());
    }



    /**
     * Merge the "value" file with the "sum" file and return the file name
     * of the merged data set.
     */
    @Override
    public TreeMap<String, String> reduce(TreeMap<String, String> value, TreeMap<String, String> sum)
    {
        if(value != null) {
            sum.putAll(value);
        }

        return sum;
    }

    @Override
    public TreeMap<String, String> treeReduce(TreeMap<String, String> lhs, TreeMap<String, String> rhs)
    {
        lhs.putAll(rhs);
        return lhs;
    }



    /**
     * Computes a reference-order sort key for a record with the given
     * UCSC chromosome name, position, and allele at that position.
     *
     * @param txChrom
     * @param position
     * @param allele
     * @param parsedRefGeneRod
     * @return
     */
    protected String computeSortKey(String txChrom, long position, char allele, RefGeneTranscriptRecord parsedRefGeneRod)
    {
        //compute chromosome sort key (this becomes the 1st column in the file and is later used to sort the output file using the GNU command-line sort utility)
        long key;
        final char lastChromChar = txChrom.charAt(txChrom.length() -1);
        switch( Character.toLowerCase(lastChromChar) ) {
        case 'm':
            key = 0;
            break;
        case 'x':
            key = 23;
            break;
        case 'y':
            key = 24;
            break;
        default:
            key = Integer.parseInt(txChrom.substring(3));
            break;
        }

        key++; //shift so there's no zero (otherwise, multiplication is screwed up in the next step)
        key += 100; //add a leading 1 to avoid having to create leading zeros
        key *= 10*1000L*1000L*1000L*1000L; //this way position doesn't overlap with chromosome keys
        key += position*1000L;
        key +=   allele;

        String result = Long.toString(key);
        if(parsedRefGeneRod != null) {
            result += (parsedRefGeneRod.name + parsedRefGeneRod.name2);
        }
        return result;
    }

    /**
     * Moves the file to the destination directory.
     *
     * @param result The file of the result
     */
    public void onTraversalDone(TreeMap<String, String> result)
    {
        if(result == null) {
            logger.info("Result is: null. No output file was generated.\n");
            return;
        }

        //move the fully merged temp file to the output file.
        logger.info("Result - writing " + result.size() + " lines to: " + outputFilename);
        BufferedWriter fileWriter  = null;
        try {
            fileWriter = new BufferedWriter(new FileWriter(outputFilename));
            fileWriter.write(OUTPUT_HEADER_LINE + '\n');
            for(String value : result.values() ) {
                fileWriter.write(value + "\n");
            }
        } catch(IOException e) {
            logger.info("Unable to write to file: " + outputFilename);
        } finally {
            try {
                if (fileWriter != null)
                    fileWriter.close();
            } catch(IOException e) {
                logger.info("Unable to write to file: " + outputFilename);
            }
        }

        logger.info("Deleting temp files..");
    }

}

