package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;

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
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;


/**
 *  Takes a refGene table ( -B arg must be: -B refgene,AnnotatorInfoTable,/path/to/refgene_file.txt) and generates the big table of nucleotides containing
 *  annotations for each possible variant at each transcript position (eg. 4 variants for each position).
 */
@Reference(window=@Window(start=-4,stop=4))
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE}, referenceMetaData={ @RMD(name="refgene",type=AnnotatorROD.class) })
public class TranscriptToInfo extends RodWalker<Integer, Integer>
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



    public static final int OUTPUT_STREAM_BUFFER_SIZE = 8*1024*1024;

    /** Output columns */
    public static final String OUTPUT_SORT_KEY = "00000_SORT_KEY";

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

    public static final String OUTPUT_SPLICE_DISTANCE = "minSpliceDistance"; //eg. 40

    public static final String OUTPUT_CODON_NUMBER = "codonCoord"; //eg. 20
    public static final String OUTPUT_REFERENCE_CODON = "referenceCodon";
    public static final String OUTPUT_REFERENCE_AA = "referenceAA";
    public static final String OUTPUT_VARIANT_CODON = "variantCodon";
    public static final String OUTPUT_VARIANT_AA = "variantAA";

    public static final String OUTPUT_CHANGES_AMINO_ACID = "changesAA"; //eg. true

    public static final String OUTPUT_CODING_COORD_STR = "coding_coord_str";
    public static final String OUTPUT_PROTEIN_COORD_STR = "protein_coord_str";
    public static final String OUTPUT_INTRON_COORD_STR = "intron_coord_str";

    public static final String[] OUTPUT_COLUMN_NAMES = new String[] {
        OUTPUT_SORT_KEY,

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

        OUTPUT_CODING_COORD_STR,
        OUTPUT_PROTEIN_COORD_STR,
        OUTPUT_INTRON_COORD_STR,

        OUTPUT_IN_CODING_REGION,
    };


    /**
     * Container for all data fields from a single refGene.txt row.
     */
    static class RefGeneTranscriptRecord
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
            if(strandStr.equals("+")) {
                positiveStrand = true;
            } else if(strandStr.equals("-")) {
                positiveStrand = false;
            } else {
                throw new IllegalArgumentException("refGene record contains unexpected strand value: \"" + strandStr + "\"");
            }

            name = refGeneRod.get("name");
            name2 = refGeneRod.get("name2");

            //String txStartStr = refGeneRod.get(REFGENE_TXSTART);  //These fields are used to generate column 1 of the ROD file (eg. they get turned into chr:txStart-txStop)
            //String txEndStr = refGeneRod.get(REFGENE_TXEND);

            final GenomeLoc loc = refGeneRod.getLocation();
            txChrom = loc.getContig();
            txStart = loc.getStart();
            txEnd = loc.getStop();

            txSequence = new StringBuilder( (int) (txEnd - txStart + 1) );  //the sequence of the entire transcript in order from 5' to 3'
            cdsSequence = new StringBuilder( (int) (cdsEnd - cdsStart + 1) ); //TODO reduce init size by size of introns

            String cdsStartStr = refGeneRod.get(REFGENE_CDS_START);
            String cdsEndStr = refGeneRod.get(REFGENE_CDS_END);

            cdsStart = Long.parseLong(cdsStartStr);
            cdsEnd = Long.parseLong(cdsEndStr);

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
         * Returns 1 when down-stream of splice-juction, -1 when upstream
         */
        public int distanceToSpliceSite(final long genomPosition) {
            int prevDistance = Integer.MAX_VALUE;
            int curDistance;
            for(int i = 0; i < exonStarts.length; i++) {
                final long curStart = exonStarts[i];
                curDistance = (int) (curStart-genomPosition);
                if(genomPosition < curStart) {
                    return Math.min(prevDistance, curDistance);
                } else {
                    prevDistance = (int) (genomPosition - curStart) + 1;
                }

                final long curStop = exonEnds[i];
                curDistance = (int) (curStop-genomPosition) + 1;
                if(genomPosition <= curStop) {
                    return Math.min(prevDistance, curDistance);
                } else {
                    prevDistance = (int) (genomPosition - curStop);
                }
            }

            return prevDistance; //out of exons. return genomPosition-curStop
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

    }


    protected int transcriptsProcessedCounter = 0;


    /** Possible values for the "POSITION_TYPE" output column. */
    enum PositionType {
        intergenic, intron, utr5, CDS, utr3, non_coding_exon, non_coding_intron
    }

    /** Output stream for computed results */
    private PrintStream outputStream;


    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        //init the output file
        String outputFilename = getToolkit().getArguments().outFileName;
        if(outputFilename == null ) {
            outputStream = System.out;
        } else {
            try {
                outputStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFilename), OUTPUT_STREAM_BUFFER_SIZE));
            } catch(Exception e) {
                throw new StingException("Unable to open output file for writing: " + outputFilename);
            }
        }

        outputStream.println(Utils.join(AnnotatorROD.DELIMITER, OUTPUT_COLUMN_NAMES));
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }


    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }



    /**
     * For each site of interest, generate the appropriate fields.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;


        Collection<RODRecordList> rods = tracker.getBoundRodTracks();
        if(rods.size() == 0) {
            //if there's nothing overlapping this locus, skip it.
            return 0;
        }




        final List<Object> refGeneRODs = (List<Object>) tracker.getReferenceMetaData("refgene");

        for(Object refGeneRodObject : refGeneRODs) {
            if(! (refGeneRodObject instanceof AnnotatorROD) ) {
                throw new RuntimeException("The 'refgene' -B arg must have the type: 'AnnotatorInputTable'.");
            }

            AnnotatorROD refGeneRod = (AnnotatorROD) refGeneRodObject;
            RefGeneTranscriptRecord parsedRefGeneRod = (RefGeneTranscriptRecord) refGeneRod.getTemporaryAttribute("parsedRefGeneRod");
            if( parsedRefGeneRod == null ) {
                parsedRefGeneRod = new RefGeneTranscriptRecord(refGeneRod);
                refGeneRod.setTemporaryAttribute("parsedRefGeneRod", parsedRefGeneRod);
            }

            if(parsedRefGeneRod.positiveStrand) {
                parsedRefGeneRod.txSequence.append(ref.getBase());
            } else {
                parsedRefGeneRod.txSequence.insert(0, ref.getBase());
            }

            final long position = ref.getLocus().getStart();
            if(parsedRefGeneRod.isProteinCodingTranscript() && position >= parsedRefGeneRod.cdsStart && position <= parsedRefGeneRod.cdsEnd && parsedRefGeneRod.isWithinExon(position) ) {
                if(parsedRefGeneRod.positiveStrand) {
                    parsedRefGeneRod.cdsSequence.append(ref.getBase());
                } else {
                    parsedRefGeneRod.cdsSequence.insert(0, ref.getBase());
                }
            }

            if(position == parsedRefGeneRod.txEnd) {
                //we've reached the end of the transcript - compute all data and write it out.
                generateOutputRecordsForROD(parsedRefGeneRod);

                transcriptsProcessedCounter++;
                if(transcriptsProcessedCounter % 100 == 0) {
                    System.err.println(new Date() + ": " +  transcriptsProcessedCounter + " transcripts processed");
                }
            }
        }

        return 1;
    }


    private static final char[] ALLELES = {'A','C','G','T'};



    private void generateOutputRecordsForROD(RefGeneTranscriptRecord parsedRefGeneRod) {

        //Transcripts that don't produce proteins are indicated in refGene by cdsStart == cdsEnd
        //These will be handled by generating only one record, with haplotypeAlternate == "*".
        final boolean isProteinCodingTranscript = parsedRefGeneRod.isProteinCodingTranscript();

        final boolean positiveStrand = parsedRefGeneRod.positiveStrand; //shorten the name

        if(isProteinCodingTranscript && parsedRefGeneRod.cdsSequence.length() % 3 != 0) {
            System.err.println("WARNING: The following record has " + parsedRefGeneRod.cdsSequence.length() + " nucleotides in its CDS region, which is not divisible by 3. Skipping... \n" + parsedRefGeneRod.toString());
            //TODO figure out why all CDS lengths aren't in multiples of 3 nucleotides - Sarah says its ok to thow these out.
            return;
        }


        //here, 5' and 3' are according to the DNA (not on the RNA)
        int frame = 0; //the frame of the current position
        int txOffset_from5 = 1; //goes from txStart 5' to txEnd 3' for both + and - strand
        int proteinOffset_from5 = 1; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand
        int codonOffset_from5 = 0; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand - counts the number of codons (this count is 1-based even though its set to 0 here.)

        char[] currentCodon_5to3 = null; //holds the current RNA codon - 5' to 3'


        final long txStart_5prime = positiveStrand ? parsedRefGeneRod.txStart : parsedRefGeneRod.txEnd; //1-based, inclusive
        final long txEnd_3prime = positiveStrand ? parsedRefGeneRod.txEnd : parsedRefGeneRod.txStart; //1-based, inclusive
        final int increment_5to3 = positiveStrand ? 1 : -1;
        final int strandSign = increment_5to3; //alias

        final long cdsStart_5prime = positiveStrand ? parsedRefGeneRod.cdsStart : parsedRefGeneRod.cdsEnd; //1-based, inclusive
        final long cdsEnd_3prime = positiveStrand ? parsedRefGeneRod.cdsEnd : parsedRefGeneRod.cdsStart ; //1-based, inclusive

        //compute chromosome key (later used to sort the output file)
        boolean mitochondrial = false;
        long chromKey;
        final char lastChromChar = parsedRefGeneRod.txChrom.charAt(parsedRefGeneRod.txChrom.length() -1);
        switch( Character.toLowerCase(lastChromChar) ) {
        case 'm':
            chromKey = 0;
            mitochondrial = true;
            break;
        case 'x':
            chromKey = 23;
            break;
        case 'y':
            chromKey = 24;
            break;
        default:
            chromKey = Integer.parseInt(parsedRefGeneRod.txChrom.substring(3));
            break;
        }

        chromKey++; //shift so there's no zero (otherwise, multiplication is screwed up in the next step)
        chromKey *= 100000000000L;


        for(long txCoord_5to3 = txStart_5prime; txCoord_5to3 != txEnd_3prime + increment_5to3; txCoord_5to3 += increment_5to3)
        {
            //get the reference sequence
            final char haplotypeReference = parsedRefGeneRod.txSequence.charAt( txOffset_from5 - 1 );

            //whether the current txCoord_5to3 is within an exon
            final boolean isWithinExon = parsedRefGeneRod.isWithinExon(txCoord_5to3); //TODO if necessary, this can be sped up by keeping track of current exon/intron

            //figure out what region the current position is in.
            PositionType positionType;
            if(isProteinCodingTranscript)
            {
                if(isWithinExon)
                {
                    if( strandSign*(txCoord_5to3 - cdsStart_5prime) < 0 ) {
                        positionType = PositionType.utr5;
                    } else if( strandSign*(txCoord_5to3 - cdsEnd_3prime) > 0 ) {
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

            //compute current codon
            if(positionType == PositionType.CDS && frame == 0)
            {
                if(currentCodon_5to3 == null) {
                    currentCodon_5to3 = new char[3];
                }

                codonOffset_from5++;

                currentCodon_5to3[0] = parsedRefGeneRod.cdsSequence.charAt( proteinOffset_from5 - 1 ); //subtract 1 to go to zero-based coords
                currentCodon_5to3[1] = parsedRefGeneRod.cdsSequence.charAt( proteinOffset_from5 );
                currentCodon_5to3[2] = parsedRefGeneRod.cdsSequence.charAt( proteinOffset_from5 + 1);

                if(!positiveStrand) {
                    currentCodon_5to3[0] = BaseUtils.simpleComplement(currentCodon_5to3[0]);
                    currentCodon_5to3[1] = BaseUtils.simpleComplement(currentCodon_5to3[1]);
                    currentCodon_5to3[2] = BaseUtils.simpleComplement(currentCodon_5to3[2]);
                }
            }


            final String distanceToSpliceSiteString = Integer.toString(parsedRefGeneRod.distanceToSpliceSite(txCoord_5to3));

            for(final char haplotypeAlternate : ALLELES )
            {
                final LinkedHashMap<String, String> outputLineFields = new LinkedHashMap<String, String>();

                outputLineFields.put(OUTPUT_SORT_KEY,  Long.toString( chromKey + txCoord_5to3) );

                outputLineFields.put(OUTPUT_CHRPOS, parsedRefGeneRod.txChrom + ":" + txCoord_5to3);
                outputLineFields.put(OUTPUT_HAPLOTYPE_REFERENCE, Character.toString( haplotypeReference ) );
                outputLineFields.put(OUTPUT_HAPLOTYPE_ALTERNATE, isProteinCodingTranscript ? Character.toString( haplotypeAlternate ) : "*");
                outputLineFields.put(OUTPUT_HAPLOTYPE_STRAND, positiveStrand ? "+" : "-");
                outputLineFields.put(OUTPUT_GENE_NAME, parsedRefGeneRod.name2 );
                outputLineFields.put(OUTPUT_TRANSCRIPT_ACCESSION, parsedRefGeneRod.name );

                outputLineFields.put(OUTPUT_POSITION_TYPE, positionType.toString() );
                outputLineFields.put(OUTPUT_MRNA_COORD, Integer.toString(txOffset_from5) );
                outputLineFields.put(OUTPUT_SPLICE_DISTANCE, distanceToSpliceSiteString );

                if(isProteinCodingTranscript)
                {
                    outputLineFields.put(OUTPUT_IN_CODING_REGION, Boolean.toString(positionType == PositionType.CDS) );

                    if(isWithinExon) {
                        outputLineFields.put(OUTPUT_CODING_COORD, Integer.toString(proteinOffset_from5) );

                        if(positionType == PositionType.CDS) {
                            final String referenceCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];
                            outputLineFields.put(OUTPUT_FRAME, Integer.toString(frame) );
                            outputLineFields.put(OUTPUT_CODON_NUMBER, Integer.toString(codonOffset_from5) );
                            outputLineFields.put(OUTPUT_PROTEIN_COORD_STR, "p."+ Integer.toString(codonOffset_from5) );
                            outputLineFields.put(OUTPUT_CODING_COORD_STR, "c."+ Integer.toString(proteinOffset_from5) );

                            outputLineFields.put(OUTPUT_REFERENCE_CODON, referenceCodon );
                            final String refAA = AminoAcidTable.getAA3LetterCode(referenceCodon, mitochondrial);
                            outputLineFields.put(OUTPUT_REFERENCE_AA, refAA);

                            final char temp = currentCodon_5to3[frame];
                            currentCodon_5to3[frame] = haplotypeAlternate;
                            final String variantCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];
                            final String variantAA = AminoAcidTable.getAA3LetterCode(variantCodon, mitochondrial);
                            outputLineFields.put(OUTPUT_VARIANT_CODON, variantCodon );
                            outputLineFields.put(OUTPUT_VARIANT_AA, variantAA);

                            currentCodon_5to3[frame] = temp;

                            outputLineFields.put(OUTPUT_CHANGES_AMINO_ACID, Boolean.toString(!refAA.equals(variantAA)));
                        }
                    }
                }


                //write out the record
                StringBuilder outputLine = new StringBuilder();
                for( String column : OUTPUT_COLUMN_NAMES) {
                    if(outputLine.length() != 0) {
                        outputLine.append(AnnotatorROD.DELIMITER);
                    }
                    String value = outputLineFields.get(column);
                    outputLine.append(value != null ? value : "");
                }

                outputStream.println(outputLine);

                if( !isProteinCodingTranscript ) {
                    //need only one record for this position, instead of 4 (one for each allele)
                    break;
                }

            } //ALLELE for-loop


            txOffset_from5++;

            if(positionType == PositionType.CDS) {
                proteinOffset_from5++;
                frame = (frame + 1) % 3;
            }



        } // l for-loop

    } //method close


    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        out.printf("Processed %d loci.\n", result);

        outputStream.close();
    }



}

