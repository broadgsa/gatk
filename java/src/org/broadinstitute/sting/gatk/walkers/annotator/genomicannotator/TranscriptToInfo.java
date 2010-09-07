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

package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.util.Interval;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableCodec;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
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
import org.broadinstitute.sting.utils.StingException;


/**
 *  Takes a table of transcripts (eg. UCSC refGene, knownGene, and CCDS tables) and generates the big table which contains
 *  annotations for each possible variant at each transcript position (eg. 4 variants at each genomic position).
 *
 *  Required args:
 *  -B - specifies the input file (ex. -B transcripts,AnnotatorInputTable,/path/to/transcript_table_file.txt)
 *  -n - Specifies which column(s) from the transcript table contain the gene name(s). (ex. -n name,name2  (for the UCSC refGene table))
 *       WARNING: The gene names for each record, when taken together, should provide a unique id for that record relative to all other records in the file.
 *
 *
 *  The map & reduce types are both TreeMap<String, String>.
 *  Each TreeMap entry represents one line in the output file. The TreeMap key is a combination of a given output line's position (so that this key can be used to sort all output lines
 *  by reference order), as well as allele and gene names (so that its unique across all output lines). The String value is the output line itself.
 */
@Reference(window=@Window(start=-4,stop=4))
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE}, referenceMetaData={ @RMD(name="transcripts",type=AnnotatorInputTableFeature.class) } )
public class TranscriptToInfo extends RodWalker<TreeMap<String, String>, TreeMap<String, String>> implements TreeReducible<TreeMap<String, String>>
{
    private static final String ROD_NAME = "transcripts";

    //@Argument(fullName="pass-through", shortName="t", doc="Optionally specifies which columns from the transcript table should be copied verbatim (aka. passed-through) to the records in the output table. For example, -B transcripts,AnnotatorInputTable,/data/refGene.txt -t id will cause the refGene id column to be copied to the output table.", required=false)
    //protected String[] PASS_THROUGH_COLUMNS = {};

    @Output
    private PrintWriter out;

    @Argument(fullName="unique-gene-name-columns", shortName="n", doc="Specifies which column(s) from the transcript table contains the gene name(s). For example, -B transcripts,AnnotatorInputTable,/data/refGene.txt -n name,name2  specifies that the name and name2 columns are gene names. WARNING: the gene names for each record, when taken together, should provide a unique id for that record relative to all other records in the file. If this is not the case, an error will be thrown. ", required=true)
    private String[] GENE_NAME_COLUMNS = {};

    private final char[] ALLELES = {'A','C','G','T'};

    /** Output columns */
    private final String OUTPUT_CHRPOS = "chrpos";
    private final String OUTPUT_HAPLOTYPE_REFERENCE = GenomicAnnotation.HAPLOTYPE_REFERENCE_COLUMN;
    private final String OUTPUT_HAPLOTYPE_ALTERNATE = GenomicAnnotation.HAPLOTYPE_ALTERNATE_COLUMN;
    private final String OUTPUT_HAPLOTYPE_STRAND = GenomicAnnotation.HAPLOTYPE_STRAND_COLUMN;

    private final String OUTPUT_IN_CODING_REGION = "inCodingRegion"; //eg. true

    private final String OUTPUT_FRAME = "frame"; //eg. 0,1,2
    private final String OUTPUT_POSITION_TYPE = "positionType"; //eg. utr5, cds, utr3, intron, intergenic

    private final String OUTPUT_MRNA_COORD = "mrnaCoord"; //1-based offset within the transcript

    private final String OUTPUT_SPLICE_DISTANCE = "spliceDist"; //eg. integer, bp to nearest exon/intron boundary

    private final String OUTPUT_CODON_NUMBER = "codonCoord"; //eg. 20
    private final String OUTPUT_REFERENCE_CODON = "referenceCodon";
    private final String OUTPUT_REFERENCE_AA = "referenceAA";
    private final String OUTPUT_VARIANT_CODON = "variantCodon";
    private final String OUTPUT_VARIANT_AA = "variantAA";

    private final String OUTPUT_CHANGES_AMINO_ACID = "changesAA"; //eg. true
    private final String OUTPUT_FUNCTIONAL_CLASS = "functionalClass"; //eg. missense

    private final String OUTPUT_CODING_COORD_STR = "codingCoordStr";
    private final String OUTPUT_PROTEIN_COORD_STR = "proteinCoordStr";

    private final String OUTPUT_SPLICE_INFO = "spliceInfo"; //(eg "splice-donor -4", or "splice-acceptor 3") for the 10bp surrounding each exon/intron boundary
    private final String OUTPUT_UORF_CHANGE = "uorfChange"; // (eg +1 or -1, indicating the addition or interruption of an ATG trinucleotide in the annotated utr5)

    private String OUTPUT_FILE_HEADER_LINE;

    //This list specifies the order of output columns in the big table.
    private final List<String> outputColumnNames = new LinkedList<String>();

    //Used to verify that the gene name portion of the key computed by computeSortKey is unique across all records.
    //This is to avoid silently loosing data in the TreeMaps.
    private final Map<String, TranscriptTableRecord> keyChecker = new HashMap<String, TranscriptTableRecord>();

    private int transcriptsProcessedCounter = 0;

    private long transcriptsThatDontStartWithMethionineOrEndWithStopCodonCounter = 0;
    private long transcriptsThatDontStartWithMethionineCounter = 0;
    private long transcriptsThatDontEndWithStopCodonCounter = 0;
    private long skippedTranscriptCounter = 0;

    private long skippedPositionsCounter = 0;
    private long totalPositionsCounter = 0;

    /** Possible values for the "POSITION_TYPE" output column. */
    private enum PositionType {
        intergenic, intron, utr5, CDS, utr3, non_coding_exon, non_coding_intron
    }


    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        //parse the GENE_NAME_COLUMNS arg and validate the column names
        final List<String> parsedGeneNameColumns = new LinkedList<String>();
        for(String token : GENE_NAME_COLUMNS) {
            parsedGeneNameColumns.addAll(Arrays.asList(token.split(",")));
        }
        GENE_NAME_COLUMNS = parsedGeneNameColumns.toArray(GENE_NAME_COLUMNS);


        ReferenceOrderedDataSource transcriptsDataSource = null;
        for(ReferenceOrderedDataSource ds : getToolkit().getRodDataSources()) {
            if(ds.getName().equals(ROD_NAME)) {
                transcriptsDataSource = ds;
                break;
            }
        }

        final ArrayList<String> header;
        try {
            header = AnnotatorInputTableCodec.readHeader(transcriptsDataSource.getReferenceOrderedData().getFile());
        } catch(Exception e) {
            throw new StingException("Failed when attempting to read header from file: " + transcriptsDataSource.getReferenceOrderedData().getFile(), e);
        }

        for(String columnName : GENE_NAME_COLUMNS) {
            if(!header.contains(columnName)) {
                throw new StingException("The column name '" + columnName + "' provided to -n doesn't match any of the column names in: " + transcriptsDataSource.getReferenceOrderedData().getFile());
            }
        }

        //init outputColumnNames list
        outputColumnNames.addAll( Arrays.asList(new String[] { OUTPUT_CHRPOS,OUTPUT_HAPLOTYPE_REFERENCE,OUTPUT_HAPLOTYPE_ALTERNATE,OUTPUT_HAPLOTYPE_STRAND,}) );
        outputColumnNames.addAll( Arrays.asList(GENE_NAME_COLUMNS) );
        outputColumnNames.addAll( Arrays.asList(new String[] {
                OUTPUT_POSITION_TYPE,

                OUTPUT_FRAME,
                OUTPUT_MRNA_COORD,
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
        }) );

        //init OUTPUT_HEADER_LINE
        StringBuilder outputHeaderLine = new StringBuilder();
        for( final String column : outputColumnNames ) {
            if(outputHeaderLine.length() != 0) {
                outputHeaderLine.append( AnnotatorInputTableCodec.DELIMITER );
            }
            outputHeaderLine.append(column);
        }

        OUTPUT_FILE_HEADER_LINE = outputHeaderLine.toString();

    }

    /**
     * Create an empty TreeMap.
     *
     * @return A new empty TreeMap.
     */
    public TreeMap<String, String> reduceInit()
    {
        return new TreeMap<String, String>() {
            @Override
            public String toString() {
                return "TreeMap of size: " + size();
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


        final List<Object> transcriptRODs = (List<Object>) tracker.getReferenceMetaData(ROD_NAME);

        //there may be multiple transcriptRODs that overlap this locus
        for(Object transcriptRodObject : transcriptRODs)
        {
            //parse this ROD if it hasn't been already.
            final AnnotatorInputTableFeature transcriptRod = (AnnotatorInputTableFeature) transcriptRodObject;
            TranscriptTableRecord parsedTranscriptRod = (TranscriptTableRecord) transcriptRod.getTemporaryAttribute("parsedTranscriptRod");
            if( parsedTranscriptRod == null ) {
                parsedTranscriptRod = new TranscriptTableRecord(transcriptRod, GENE_NAME_COLUMNS);
                transcriptRod.setTemporaryAttribute("parsedTranscriptRod", parsedTranscriptRod);
            }

            //populate parsedTranscriptRod.txSequence
            if(parsedTranscriptRod.positiveStrand) {
                parsedTranscriptRod.txSequence.append(ref.getBaseAsChar());
            } else {
                final char complementBase = BaseUtils.simpleComplement(ref.getBaseAsChar());
                parsedTranscriptRod.txSequence.insert(0, complementBase);
            }

            //populate parsedTranscriptRod.utr5Sequence and parsedTranscriptRod.cdsSequence
            final int position = (int) ref.getLocus().getStart();
            if(parsedTranscriptRod.isProteinCodingTranscript() && parsedTranscriptRod.isWithinExon(position) )
            {
                //we're within an exon of a proteinCodingTranscript

                if(parsedTranscriptRod.positiveStrand)
                {
                    if(position < parsedTranscriptRod.cdsStart)
                    {
                        parsedTranscriptRod.utr5Sequence.append(ref.getBaseAsChar()); //within utr5
                    }
                    else if(position >= parsedTranscriptRod.cdsStart && position <= parsedTranscriptRod.cdsEnd)
                    {
                        parsedTranscriptRod.cdsSequence.append(ref.getBaseAsChar()); //within CDS
                    }
                }
                else
                {
                    final char complementBase = BaseUtils.simpleComplement(ref.getBaseAsChar());
                    if(position > parsedTranscriptRod.cdsEnd)
                    {
                        //As we move left to right (aka. 3' to 5'), we do insert(0,..) to reverse the sequence so that it become 5' to 3' in parsedTranscriptRod.utr5Sequence.
                        parsedTranscriptRod.utr5Sequence.insert(0,complementBase); //within utr5.
                    }
                    else if(position >= parsedTranscriptRod.cdsStart && position <= parsedTranscriptRod.cdsEnd)
                    {
                        parsedTranscriptRod.cdsSequence.insert(0,complementBase); //within CDS
                    }
                }
            }

            if(position == parsedTranscriptRod.txEnd)
            {

                checkGeneNamesForUniqueness(parsedTranscriptRod);

                //we've reached the end of the transcript - compute all data and write it out.
                try
                {
                    generateOutputRecordsForROD(parsedTranscriptRod, result);
                }
                catch(Exception e)
                {
                    throw new RuntimeException(Thread.currentThread().getName() + " - Unexpected error occurred at position: [" + parsedTranscriptRod.txChrom + ":" + position + "] in transcript: " + parsedTranscriptRod, e);
                }

                transcriptsProcessedCounter++;
                if(transcriptsProcessedCounter % 100 == 0) {
                    logger.info(new Date() + ": " +  transcriptsProcessedCounter + " transcripts processed");
                }
            }
        }

        return result;
    }



    /**
     * Verify that gene names, when taken together, provide a unique key for this record
     * this guarantees that the computeSortKey(..) method will work as expected.
     * If this is not the case, disambiguate the gene names.
     *
     * @param parsedTranscriptRod
     */
    private void checkGeneNamesForUniqueness(TranscriptTableRecord parsedTranscriptRod)
    {
        final String geneNamePortionOfSortKey = computeGeneNamePortionOfSortKey(parsedTranscriptRod.geneNames);
        final TranscriptTableRecord collisionRecord = keyChecker.get( geneNamePortionOfSortKey );
        if(collisionRecord != null /* && new Interval(     - don't allow collitions even if positions are different.
                collisionRecord.txChrom,
                (int) collisionRecord.txStart,
                (int) collisionRecord.txEnd).intersects( new Interval(
                        parsedTranscriptRod.txChrom,
                        (int) parsedTranscriptRod.txStart,
                        (int) parsedTranscriptRod.txEnd)) */) {

            //disambiguate the gene names of parsedTranscriptRod
            int counter = 1;
            String newGeneNamePortionOfSortKey = null;
            String[] newGeneNames = null;
            do {
                newGeneNames = parsedTranscriptRod.geneNames.clone();
                newGeneNames[0] += ("." +  ++counter); //append '.1' or similar until you find a key that doesn't collide
                newGeneNamePortionOfSortKey = computeGeneNamePortionOfSortKey(newGeneNames);
            } while(keyChecker.containsKey(newGeneNamePortionOfSortKey));

            parsedTranscriptRod.geneNames = newGeneNames;

            keyChecker.put(newGeneNamePortionOfSortKey, parsedTranscriptRod );

            logger.warn("WARNING: The gene names of the following two records were identical, so they have been disambiguated:\n 1:" + collisionRecord + "\n 2:" + parsedTranscriptRod +".\n Otherwise, output data would likely be lost due to collisions in the keys generated by computeSortKey(..). ");
        } else {
            keyChecker.put(geneNamePortionOfSortKey, parsedTranscriptRod );
        }
    }




    /**
     * Returns the filename.
     *
     * @param parsedTranscriptRod
     * @param result
     * @return
     */
    private void generateOutputRecordsForROD(TranscriptTableRecord parsedTranscriptRod, TreeMap<String, String> result) throws IOException
    {
        //Transcripts that don't produce proteins are indicated in transcript by cdsStart == cdsEnd
        //These will be handled by generating only one record, with haplotypeAlternate == "*".
        final boolean isProteinCodingTranscript = parsedTranscriptRod.isProteinCodingTranscript();
        final boolean isMitochondrialTranscript = parsedTranscriptRod.txChrom.toLowerCase().startsWith("chrm");

        final boolean positiveStrand = parsedTranscriptRod.positiveStrand; //alias


        if(isProteinCodingTranscript && parsedTranscriptRod.cdsSequence.length() % 3 != 0) {
            if (!isMitochondrialTranscript) {
                logger.error("ERROR: Transcript " + parsedTranscriptRod +" at position ["+ parsedTranscriptRod.txChrom + ":" +parsedTranscriptRod.txStart + "-" + parsedTranscriptRod.txEnd + "] has " + parsedTranscriptRod.cdsSequence.length() + " nucleotides in its CDS region, which is not divisible by 3. Skipping...");
                //discard transcripts where CDS length is not a multiple of 3
                skippedTranscriptCounter++;
                return;
            } else {

                //In mitochondrial genes, the polyA tail may complete the stop codon, allowing transcript . To check for this special case:
                //1. check that the CDS covers the entire transcript
                //2. add 1 or 2 A's to the 3' end of the transcript (as needed to make it divisible by 3)
                //3. check whether the last 3 letters now form a stop codon using the mitochondrial AA table
                //4. If not, skip this gene, else incorporate the A's and process it like any other gene.

                if( parsedTranscriptRod.txSequence.length() == parsedTranscriptRod.cdsSequence.length()) {
                    do { //append A's until sequence length is divisible by 3
                        parsedTranscriptRod.txSequence.append('*');
                        parsedTranscriptRod.cdsSequence.append('a');
                        if(positiveStrand) {
                            parsedTranscriptRod.txEnd++;
                            parsedTranscriptRod.cdsEnd++;
                            parsedTranscriptRod.exonEnds[0]++;
                        } else {
                            parsedTranscriptRod.txStart--;
                            parsedTranscriptRod.cdsStart--;
                            parsedTranscriptRod.exonStarts[0]--;
                        }
                    } while( parsedTranscriptRod.cdsSequence.length() % 3 != 0);

                } else {
                    logger.error("ERROR: Mitochnodrial transcript " + parsedTranscriptRod +" at position ["+ parsedTranscriptRod.txChrom + ":" +parsedTranscriptRod.txStart + "-" + parsedTranscriptRod.txEnd + "] has " + parsedTranscriptRod.cdsSequence.length() + " nucleotides in its CDS region, which is not divisible by 3. The CDS does not cover the entire transcript, so its not possible to use A's from the polyA tail. Skipping...");
                    skippedTranscriptCounter++;
                    return;
                }
            }
        }


        //warn if the first codon isn't Methionine and/or the last codon isn't a stop codon.
        if(isProteinCodingTranscript) {
            final int cdsSequenceLength = parsedTranscriptRod.cdsSequence.length();

            final String firstCodon = parsedTranscriptRod.cdsSequence.substring(0, 3);
            final AminoAcid firstAA = isMitochondrialTranscript ? AminoAcidTable.getMitochondrialAA( firstCodon, true ) : AminoAcidTable.getEukaryoticAA( firstCodon ) ;

            final String lastCodon = parsedTranscriptRod.cdsSequence.substring(cdsSequenceLength - 3, cdsSequenceLength);
            final AminoAcid lastAA = isMitochondrialTranscript ? AminoAcidTable.getMitochondrialAA( lastCodon, false ) : AminoAcidTable.getEukaryoticAA( lastCodon ) ;

            if( firstAA != AminoAcidTable.METHIONINE && !lastAA.isStop()) {
                transcriptsThatDontStartWithMethionineOrEndWithStopCodonCounter++;
                logger.warn("WARNING: The CDS of transcript " + parsedTranscriptRod.geneNames[0] +" at position ["+ parsedTranscriptRod.txChrom + ":" +parsedTranscriptRod.txStart + "-" + parsedTranscriptRod.txEnd + "] does not start with Methionine or end in a stop codon. The first codon is: " + firstCodon + " (" + firstAA + "). The last codon is: " + lastCodon + " (" + lastAA + "). NOTE: This is just a warning - the transcript will be included in the output.");
            } else if( firstAA != AminoAcidTable.METHIONINE) {
                transcriptsThatDontStartWithMethionineCounter++;
                logger.warn("WARNING: The CDS of transcript " + parsedTranscriptRod.geneNames[0] +" at position ["+ parsedTranscriptRod.txChrom + ":" +parsedTranscriptRod.txStart + "-" + parsedTranscriptRod.txEnd + "] does not start with Methionine. The first codon is: " + firstCodon + " (" + firstAA + "). NOTE: This is just a warning - the transcript will be included in the output.");
            } else if(!lastAA.isStop()) {
                transcriptsThatDontEndWithStopCodonCounter++;
                logger.warn("WARNING: The CDS of transcript " + parsedTranscriptRod.geneNames[0] +" at position ["+ parsedTranscriptRod.txChrom + ":" +parsedTranscriptRod.txStart + "-" + parsedTranscriptRod.txEnd + "] does not end in a stop codon. The last codon is: " + lastCodon + " (" + lastAA + "). NOTE: This is just a warning - the transcript will be included in the output.");
            }
        }

        final int txStart_5prime = positiveStrand ? parsedTranscriptRod.txStart : parsedTranscriptRod.txEnd; //1-based, inclusive
        final int txEnd_3prime = positiveStrand ? parsedTranscriptRod.txEnd : parsedTranscriptRod.txStart; //1-based, inclusive
        final int increment_5to3 = positiveStrand ? 1 : -1; //whether to increment or decrement
        final int strandSign = increment_5to3; //alias

        final int cdsStart_5prime = positiveStrand ? parsedTranscriptRod.cdsStart : parsedTranscriptRod.cdsEnd; //1-based, inclusive
        final int cdsEnd_3prime = positiveStrand ? parsedTranscriptRod.cdsEnd : parsedTranscriptRod.cdsStart ; //1-based, inclusive

        int frame = 0; //the frame of the current position
        int txOffset_from5 = 1; //goes from txStart 5' to txEnd 3' for both + and - strand
        int utr5Count_from5 = 0;
        int mrnaCoord_from5 = 1; //goes from txStart 5' to txEnd 3' for both + and - strand, but only counts bases within exons.
        char[] utr5NucBuffer_5to3 = null; //used to find uORFs - size = 5 because to hold the 3 codons that overlap any given position: [-2,-1,0], [-1,0,1], and [0,1,2]

        int codonCount_from5 = 1; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand - counts the number of codons - 1-based
        int codingCoord_from5 = isProteinCodingTranscript ? parsedTranscriptRod.computeInitialCodingCoord() : -1; //goes from cdsStart 5' to cdsEnd 3' for both + and - strand
        boolean codingCoordResetForCDS = false;
        boolean codingCoordResetForUtr3 = false;
        final char[] currentCodon_5to3 = isProteinCodingTranscript ? new char[3] : null; //holds the current RNA codon - 5' to 3'



        PositionType positionType = null;
        boolean isWithinIntronAndFarFromSpliceJunction = false;
        int intronStart_5prime = -1;
        int intronEnd_5prime = -1;

        final Map<String, String> outputLineFields = new HashMap<String, String>();

        for(int txCoord_5to3 = txStart_5prime; txCoord_5to3 != txEnd_3prime + increment_5to3; txCoord_5to3 += increment_5to3)
        {
            ++totalPositionsCounter;

            //compute certain attributes of the current position
            final boolean isWithinExon = parsedTranscriptRod.isWithinExon(txCoord_5to3); //TODO if necessary, this can be sped up by keeping track of current exon/intron

            final int distanceToNearestSpliceSite = parsedTranscriptRod.computeDistanceToNearestSpliceSite(txCoord_5to3);
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

                    final int intronStart = (intronStart_5prime < intronEnd_5prime ? intronStart_5prime : intronEnd_5prime) ;
                    final int intronEnd = (intronEnd_5prime > intronStart_5prime ? intronEnd_5prime : intronStart_5prime);
                    outputLineFields.clear();
                    outputLineFields.put(OUTPUT_CHRPOS, parsedTranscriptRod.txChrom + ":" + intronStart +  "-" + intronEnd);
                    outputLineFields.put(OUTPUT_HAPLOTYPE_REFERENCE, Character.toString( '*' ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_ALTERNATE, Character.toString( '*' ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_STRAND, positiveStrand ? "+" : "-");
                    for(int i = 0; i < GENE_NAME_COLUMNS.length; i++) {
                        outputLineFields.put(GENE_NAME_COLUMNS[i], parsedTranscriptRod.geneNames[i] );
                    }


                    outputLineFields.put(OUTPUT_POSITION_TYPE, positionType.toString() );

                    if(isProteinCodingTranscript)
                    {
                        outputLineFields.put(OUTPUT_IN_CODING_REGION, Boolean.toString(positionType == PositionType.CDS) );
                    }

                    final String sortKey = computeSortKey(parsedTranscriptRod.txChrom, intronStart, '*', parsedTranscriptRod);

                    addThisLineToResult(sortKey, outputLineFields, result);
                }


                //when in utr5, compute the utr5NucBuffer_5to3 which is later used to compute the OUTPUT_UORF_CHANGE field
                if(positionType == PositionType.utr5)
                {
                    if(utr5Count_from5 < parsedTranscriptRod.utr5Sequence.length())
                    {
                        if(utr5NucBuffer_5to3 == null) {
                            //initialize
                            utr5NucBuffer_5to3 = new char[5];
                            utr5NucBuffer_5to3[3] = parsedTranscriptRod.utr5Sequence.charAt( utr5Count_from5 );

                            if(utr5Count_from5 + 1 < parsedTranscriptRod.utr5Sequence.length() ) {
                                utr5NucBuffer_5to3[4] = parsedTranscriptRod.utr5Sequence.charAt( utr5Count_from5 + 1 );
                            }
                        }

                        //as we move 5' to 3', shift nucleotides down to the 5' end, making room for the new 3' nucleotide:
                        utr5NucBuffer_5to3[0] = utr5NucBuffer_5to3[1];
                        utr5NucBuffer_5to3[1] = utr5NucBuffer_5to3[2];
                        utr5NucBuffer_5to3[2] = utr5NucBuffer_5to3[3];
                        utr5NucBuffer_5to3[3] = utr5NucBuffer_5to3[4];

                        char nextRefBase = 0;
                        if( utr5Count_from5 + 2 < parsedTranscriptRod.utr5Sequence.length() )
                        {
                            nextRefBase = parsedTranscriptRod.utr5Sequence.charAt( utr5Count_from5 + 2 );
                        }
                        utr5NucBuffer_5to3[4] = nextRefBase;

                        //check for bad bases
                        if( (utr5NucBuffer_5to3[0] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[0])) ||
                                (utr5NucBuffer_5to3[1] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[1])) ||
                                (utr5NucBuffer_5to3[2] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[2])) ||
                                (utr5NucBuffer_5to3[3] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[3])) ||
                                (utr5NucBuffer_5to3[4] != 0 && !BaseUtils.isRegularBase(utr5NucBuffer_5to3[4])))
                        {
                            logger.debug("Skipping current position [" + parsedTranscriptRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedTranscriptRod.geneNames.toString() +". utr5NucBuffer_5to3 contains irregular base:" + utr5NucBuffer_5to3[0] + utr5NucBuffer_5to3[1] + utr5NucBuffer_5to3[2] + utr5NucBuffer_5to3[3] + utr5NucBuffer_5to3[4]);// +". Transcript is: " + parsedTranscriptRod);
                            ++skippedPositionsCounter;
                            continue;
                        }

                    } else { // if(utr5Count_from5 >= parsedTranscriptRod.utr5Sequence.length())
                        //defensive programming
                        throw new RuntimeException("Exception: Skipping current position [" + parsedTranscriptRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedTranscriptRod.geneNames.toString() +". utr5Count_from5 is now " + utr5Count_from5 + ", while parsedTranscriptRod.utr5Sequence.length() == " + parsedTranscriptRod.utr5Sequence.length() + ". This means parsedTranscriptRod.utr5Sequence isn't as long as it should be. This is a bug in handling this record: " + parsedTranscriptRod);

                    }
                }


                //when in CDS, compute current codon
                if(positionType == PositionType.CDS)
                {
                    if(frame == 0)
                    {
                        currentCodon_5to3[0] = parsedTranscriptRod.cdsSequence.charAt( codingCoord_from5 - 1 ); //subtract 1 to go to zero-based coords
                        currentCodon_5to3[1] = parsedTranscriptRod.cdsSequence.charAt( codingCoord_from5 );
                        currentCodon_5to3[2] = parsedTranscriptRod.cdsSequence.charAt( codingCoord_from5 + 1);
                    }

                    //check for bad bases
                    if(!BaseUtils.isRegularBase(currentCodon_5to3[0]) || !BaseUtils.isRegularBase(currentCodon_5to3[1]) || !BaseUtils.isRegularBase(currentCodon_5to3[2])) {
                        logger.debug("Skipping current position [" + parsedTranscriptRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedTranscriptRod.geneNames.toString() +". CDS codon contains irregular base:" + currentCodon_5to3[0] + currentCodon_5to3[1] + currentCodon_5to3[2]);// +". Transcript is: " + parsedTranscriptRod);
                        ++skippedPositionsCounter;
                        continue;
                    }

                }

                char haplotypeReference = parsedTranscriptRod.txSequence.charAt( txOffset_from5 - 1 );
                if(!positiveStrand) {
                    haplotypeReference = BaseUtils.simpleComplement(haplotypeReference); //txSequence contents depend on whether its +/- strand
                }

                if(!BaseUtils.isRegularBase(haplotypeReference) && haplotypeReference != '*') { //* is special case for mitochondrial genes where polyA tail completes the last codon
                    //check for bad bases
                    logger.debug("Skipping current position [" + parsedTranscriptRod.txChrom + ":" +txCoord_5to3 + "] in transcript " + parsedTranscriptRod.geneNames.toString() + ". The reference contains an irregular base:" + haplotypeReference); // +". Transcript is: " + parsedTranscriptRod);
                    ++skippedPositionsCounter;
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
                    outputLineFields.put(OUTPUT_CHRPOS, parsedTranscriptRod.txChrom + ":" + txCoord_5to3);
                    outputLineFields.put(OUTPUT_HAPLOTYPE_REFERENCE, Character.toString( haplotypeReference ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_ALTERNATE, Character.toString( haplotypeAlternate ) );
                    outputLineFields.put(OUTPUT_HAPLOTYPE_STRAND, positiveStrand ? "+" : "-");
                    for(int i = 0; i < GENE_NAME_COLUMNS.length; i++) {
                        outputLineFields.put(GENE_NAME_COLUMNS[i], parsedTranscriptRod.geneNames[i] );
                    }

                    outputLineFields.put(OUTPUT_POSITION_TYPE, positionType.toString() );
                    if(isWithinExon) {
                        outputLineFields.put(OUTPUT_MRNA_COORD, Integer.toString(mrnaCoord_from5) );
                    }
                    outputLineFields.put(OUTPUT_SPLICE_DISTANCE, Integer.toString(distanceToNearestSpliceSite) );

                    //compute OUTPUT_SPLICE_INFO
                    final String spliceInfoString;
                    if(isWithin10bpOfSpliceJunction) {
                        if(distanceToNearestSpliceSite < 0) {
                            //is on the 5' side of the splice junction
                            if(isWithinExon) {
                                spliceInfoString = "splice-donor_" + distanceToNearestSpliceSite;
                            } else {
                                spliceInfoString = "splice-acceptor_" + distanceToNearestSpliceSite;
                            }
                        } else {
                            if(isWithinExon) {
                                spliceInfoString = "splice-acceptor_" + distanceToNearestSpliceSite;
                            } else {
                                spliceInfoString = "splice-donor_" + distanceToNearestSpliceSite;
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
                        String refCodon1 = (Character.toString(utr5NucBuffer_5to3[0]) + Character.toString(utr5NucBuffer_5to3[1]) + utr5NucBuffer_5to3[2]).toUpperCase();
                        String refCodon2 = (Character.toString(utr5NucBuffer_5to3[1]) + Character.toString(utr5NucBuffer_5to3[2]) + utr5NucBuffer_5to3[3]).toUpperCase();
                        String refCodon3 = (Character.toString(utr5NucBuffer_5to3[2]) + Character.toString(utr5NucBuffer_5to3[3]) + utr5NucBuffer_5to3[4]).toUpperCase();

                        String varCodon1 = (Character.toString(utr5NucBuffer_5to3[0]) + Character.toString(utr5NucBuffer_5to3[1]) + haplotypeAlternate).toUpperCase();
                        String varCodon2 = (Character.toString(utr5NucBuffer_5to3[1]) + Character.toString(haplotypeAlternate) + utr5NucBuffer_5to3[3]).toUpperCase();
                        String varCodon3 = (Character.toString(haplotypeAlternate) + Character.toString(utr5NucBuffer_5to3[3]) + utr5NucBuffer_5to3[4]).toUpperCase();

                        //check for +1 (eg. addition of new ATG uORF) and -1 (eg. disruption of existing ATG uORF)
                        String uORFChangeStr = null;
                        if( (refCodon1.equals("ATG") && !varCodon1.equals("ATG")) ||
                                (refCodon2.equals("ATG") && !varCodon2.equals("ATG")) ||
                                (refCodon3.equals("ATG") && !varCodon3.equals("ATG")))
                        {
                            uORFChangeStr = "-1";
                        }
                        else if((varCodon1.equals("ATG") && !refCodon1.equals("ATG")) ||
                                (varCodon2.equals("ATG") && !refCodon2.equals("ATG")) ||
                                (varCodon3.equals("ATG") && !refCodon3.equals("ATG")))
                        {
                            uORFChangeStr = "+1";
                        }

                        outputLineFields.put(OUTPUT_UORF_CHANGE, uORFChangeStr );
                    }
                    //compute CDS-specific fields
                    else if(positionType == PositionType.CDS)
                    {
                        final String referenceCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];
                        final String variantCodon = Character.toString(currentCodon_5to3[0]) + Character.toString(currentCodon_5to3[1]) + currentCodon_5to3[2];

                        final AminoAcid refAA = isMitochondrialTranscript ? AminoAcidTable.getMitochondrialAA( referenceCodon, codonCount_from5 == 1 ) : AminoAcidTable.getEukaryoticAA( referenceCodon ) ;
                        final AminoAcid variantAA = isMitochondrialTranscript ? AminoAcidTable.getMitochondrialAA( variantCodon, codonCount_from5 == 1 ) : AminoAcidTable.getEukaryoticAA( variantCodon ) ;

                        if ( refAA.isUnknown() || variantAA.isUnknown() ) {
                            logger.warn("Illegal amino acid detected: refCodon=" + referenceCodon + " altCodon=" + variantCodon);
                        }

                        outputLineFields.put(OUTPUT_FRAME, Integer.toString(frame) );
                        outputLineFields.put(OUTPUT_CODON_NUMBER, Integer.toString(codonCount_from5) );
                        outputLineFields.put(OUTPUT_REFERENCE_CODON, referenceCodon );
                        outputLineFields.put(OUTPUT_REFERENCE_AA, refAA.getCode());

                        final char temp = currentCodon_5to3[frame];
                        currentCodon_5to3[frame] = haplotypeAlternate;
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
                        codingCoordStr.append( isMitochondrialTranscript ? "m." : "c." );
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
                        final String sortKey = computeSortKey(parsedTranscriptRod.txChrom, txCoord_5to3, haplotypeAlternate, parsedTranscriptRod);

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
                if(isWithinExon) {
                    mrnaCoord_from5++;
                }

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
        for( final String column : outputColumnNames ) {
            if(outputLine.length() != 0) {
                outputLine.append( AnnotatorInputTableCodec.DELIMITER );
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
     * @param parsedTranscriptRod
     * @return
     */
    protected String computeSortKey(String txChrom, int position, char allele, TranscriptTableRecord parsedTranscriptRod)
    {
        //compute chromosome sort key (this becomes the 1st column in the file and is later used to sort the output file using the GNU command-line sort utility)
        long key;
        final char lastChromChar = txChrom.charAt(txChrom.length() -1);
        switch( Character.toLowerCase(lastChromChar) ) {
            case 'm':
            case 't': // for hg19 or b36
                key = 0;
                break;
            case 'x':
                key = 23;
                break;
            case 'y':
                key = 24;
                break;
            default:
                if( txChrom.startsWith("chr") ) {
                    key = Integer.parseInt(txChrom.substring(3));
                } else {
                    key = Integer.parseInt(txChrom);
                }
                break;
        }

        key++; //shift so there's no zero (otherwise, multiplication is screwed up in the next step)
        key += 100; //add a leading 1 to avoid having to create leading zeros
        key *= 10*1000L*1000L*1000L*1000L; //this way position doesn't overlap with chromosome keys
        key += position*1000L;
        key += allele;

        String result = Long.toString(key);
        if(parsedTranscriptRod != null) {
            result += computeGeneNamePortionOfSortKey(parsedTranscriptRod.geneNames);  //append all gene names, so that the key is unique across all possible splicing variants, etc.
        }
        return result;
    }

    /**
     * Computes the gene-name portion of the output record sort key.
     *
     * @param geneNames
     *
     * @return
     */
    private String computeGeneNamePortionOfSortKey(String[] geneNames) {
        StringBuilder result = new StringBuilder();
        for(String geneName : geneNames)
            result.append(geneName);  //append all gene names, so that the key is unique across all possible splicing variants, etc.
        return result.toString();
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

        logger.info("Writing " + result.size() + " lines to: " + out + ". Average of " + (totalPositionsCounter == 0 ? 0 : (10*result.size()/totalPositionsCounter)/10.0f) + " lines per genomic position.");
        out.println(OUTPUT_FILE_HEADER_LINE);
        for(String value : result.values() ) {
            out.println(value);
        }

        logger.info("Skipped " + skippedPositionsCounter + " in-transcript genomic positions out of "+ totalPositionsCounter + " total (" + ( totalPositionsCounter == 0 ? 0 : (100*skippedPositionsCounter)/totalPositionsCounter) + "%)");
        logger.info("Skipped " + skippedTranscriptCounter + " transcripts out of "+ transcriptsProcessedCounter + " total (" + ( transcriptsProcessedCounter == 0 ? 0 : (100*skippedTranscriptCounter)/transcriptsProcessedCounter) + "%)");
        logger.info("Protein-coding transcripts (eg. with a CDS region) that don't start with Methionine or end in a stop codon: " + transcriptsThatDontStartWithMethionineOrEndWithStopCodonCounter + " transcripts out of "+ transcriptsProcessedCounter + " total (" + ( transcriptsProcessedCounter == 0 ? 0 : (100*transcriptsThatDontStartWithMethionineOrEndWithStopCodonCounter)/transcriptsProcessedCounter) + "%)");
        logger.info("Protein-coding transcripts (eg. with a CDS region) that don't start with Methionine: " + transcriptsThatDontStartWithMethionineCounter + " transcripts out of "+ transcriptsProcessedCounter + " total (" + ( transcriptsProcessedCounter == 0 ? 0 : (100*transcriptsThatDontStartWithMethionineCounter)/transcriptsProcessedCounter) + "%)");
        logger.info("Protein-coding transcripts (eg. with a CDS region) that don't end in a stop codon: " + transcriptsThatDontEndWithStopCodonCounter + " transcripts out of "+ transcriptsProcessedCounter + " total (" + ( transcriptsProcessedCounter == 0 ? 0 : (100*transcriptsThatDontEndWithStopCodonCounter)/transcriptsProcessedCounter) + "%)");
    }


    /**
     * Container for all data fields from a single row of the transcript table.
     */
    protected static class TranscriptTableRecord
    {
        public static final String STRAND_COLUMN = "strand";  //eg. +
        public static final String CDS_START_COLUMN = "cdsStart";
        public static final String CDS_END_COLUMN = "cdsEnd";
        public static final String EXON_COUNT_COLUMN = "exonCount";
        public static final String EXON_STARTS_COLUMN = "exonStarts";
        public static final String EXON_ENDS_COLUMN = "exonEnds";
        //public static final String EXON_FRAMES_COLUMN = "exonFrames";


        /**
         * This StringBuffer accumulates the entire transcript sequence.
         * This buffer is used instead of using the GATK window mechanism
         * because arbitrary-length look-aheads and look-behinds are needed to deal
         * with codons that span splice-junctions in + & - strand transcripts.
         * The window mechanism requires hard-coding the window size, which would
         * translate into a limit on maximum supported intron size. To avoid this, the
         * sequence is accumulated as the transcript is scanned left-to-right.
         * Then, all calculations are performed at the end.
         */
        public StringBuilder txSequence;  //the sequence of the entire transcript in order from 5' to 3'
        public StringBuilder utr5Sequence; //the protein coding sequence (with introns removed) in order from 5' to 3'
        public StringBuilder cdsSequence; //the protein coding sequence (with introns removed) in order from 5' to 3'

        public boolean positiveStrand; //whether the transcript is on the + or the - strand.
        public String[] geneNames; //eg. NM_021649

        public String txChrom; //The chromosome name
        public int txStart;
        public int txEnd;

        public int cdsStart;
        public int cdsEnd;

        public int[] exonStarts;
        public int[] exonEnds;
        //public int[] exonFrames; - not used for anything, frame is computed another way




        /**
         * Constructor.
         *
         * @param transcriptRod A rod representing a single record in the transcript table.
         */
        public TranscriptTableRecord(final AnnotatorInputTableFeature transcriptRod, String[] geneNameColumns) {

            //String binStr = transcriptRod.get("bin");
            //String idStr = transcriptRod.get("id"); //int(10) unsigned range Unique identifier ( usually 0 for some reason - even for translated )
            String strandStr = transcriptRod.getColumnValue(STRAND_COLUMN);
            if(strandStr == null) {
                throw new IllegalArgumentException("Transcript table record doesn't contain a 'strand' column. Make sure the transcripts input file has a header and the usual columns: \"" + strandStr + "\"");
            } else if(strandStr.equals("+")) {
                positiveStrand = true;
            } else if(strandStr.equals("-")) {
                positiveStrand = false;
            } else {
                throw new IllegalArgumentException("Transcript table record contains unexpected value for 'strand' column: \"" + strandStr + "\"");
            }

            geneNames = new String[geneNameColumns.length];
            for(int i = 0; i < geneNameColumns.length; i++) {
                geneNames[i] = transcriptRod.getColumnValue(geneNameColumns[i]);
            }

            //String txStartStr = transcriptRod.get(TXSTART_COLUMN);  //These fields were used to generate column 1 of the ROD file (eg. they got turned into chr:txStart-txStop)
            //String txEndStr = transcriptRod.get(TXEND_COLUMN);
            txChrom = transcriptRod.getChr();
            txStart = transcriptRod.getStart();
            txEnd = transcriptRod.getEnd();

            String cdsStartStr = transcriptRod.getColumnValue(CDS_START_COLUMN);
            String cdsEndStr = transcriptRod.getColumnValue(CDS_END_COLUMN);

            cdsStart = Integer.parseInt(cdsStartStr);
            cdsEnd = Integer.parseInt(cdsEndStr);

            txSequence = new StringBuilder( (int) (txEnd - txStart + 1) );  //the sequence of the entire transcript in order from 5' to 3'
            if(isProteinCodingTranscript()) {
                utr5Sequence = new StringBuilder( positiveStrand ? (int) (cdsStart - txStart + 1) : (int) (txEnd - cdsEnd + 1) ); //TODO reduce init size by size of introns
                cdsSequence = new StringBuilder( (int) (cdsEnd - cdsStart + 1) ); //TODO reduce init size by size of introns
            }

            String exonCountStr = transcriptRod.getColumnValue(EXON_COUNT_COLUMN);
            String exonStartsStr = transcriptRod.getColumnValue(EXON_STARTS_COLUMN);
            String exonEndsStr = transcriptRod.getColumnValue(EXON_ENDS_COLUMN);
            //String exonFramesStr = transcriptRod.get(EXON_FRAMES_COLUMN);

            String[] exonStartStrs = exonStartsStr.split(",");
            String[] exonEndStrs = exonEndsStr.split(",");
            //String[] exonFrameStrs = exonFramesStr.split(",");

            int exonCount = Integer.parseInt(exonCountStr);
            if(exonCount != exonStartStrs.length || exonCount != exonEndStrs.length /* || exonCount != exonFrameStrs.length */)
            {
                throw new RuntimeException("exonCount != exonStarts.length || exonCount != exonEnds.length || exonCount != exonFrames.length. Exon starts: " + exonStartsStr + ", Exon ends: " + exonEndsStr + /*",  Exon frames: " + exonFramesStr + */", Exon count: " + exonCountStr +". transcriptRod = " + transcriptRod);
            }

            exonStarts = new int[exonCount];
            exonEnds = new int[exonCount];
            //exonFrames = new int[exonCount];
            for(int i = 0; i < exonCount; i++) {
                exonStarts[i] = Integer.parseInt(exonStartStrs[i]);
                exonEnds[i] = Integer.parseInt(exonEndStrs[i]);
                //exonFrames[i] = Integer.parseInt(exonFrameStrs[i]);
            }
        }


        /**
         * Takes a genomic position on the same contig as the transcript, and
         * returns true if this position falls within an exon.
         */
        public boolean isWithinExon(final int genomPosition) {
            for(int i = 0; i < exonStarts.length; i++) {
                final int curStart = exonStarts[i];
                if(genomPosition < curStart) {
                    return false;
                }
                final int curStop = exonEnds[i];
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
        public int computeDistanceToNearestSpliceSite(final int genomPosition) {
            int prevDistance = Integer.MAX_VALUE;
            for(int i = 0; i < exonStarts.length; i++) {
                final int curStart = exonStarts[i];
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

                final int curStop = exonEnds[i];
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
            StringBuilder sb = new StringBuilder();
            sb.append("chrpos=" + txChrom + ':' + txStart + '-' + txEnd + ", strand=" + (positiveStrand ? '+':'-') + ", gene-names=" + Arrays.toString(geneNames) + ", cds="+ cdsStart + '-' + cdsEnd + ", exonStarts=" + Arrays.toString(exonStarts) + ", exonEnds=" + Arrays.toString(exonEnds));
            return sb.toString();
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
                    final int exonStart = exonStarts[i];
                    final int exonEnd = exonEnds[i];
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
                final int cdsStart_5prime = cdsEnd;
                if(cdsStart_5prime == exonEnds[exonEnds.length - 1]) {
                    //the 1st nucleotide of the transcript is CDS.
                    return 1;
                }

                int result = 0;
                for(int i = exonEnds.length - 1; i >= 0; i--)
                {
                    final int exonStart = exonEnds[i]; //when its the negative strand, the 5' coord of the 1st exon is exonEnds[i]
                    final int exonEnd = exonStarts[i];
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


}

