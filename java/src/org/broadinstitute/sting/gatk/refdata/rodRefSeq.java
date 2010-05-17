package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 22, 2009
 * Time: 3:19:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class rodRefSeq extends BasicReferenceOrderedDatum implements Transcript {

    private String transcript_id;
    private int strand;
    private GenomeLoc transcript_interval;
    private GenomeLoc transcript_coding_interval;
    private List<GenomeLoc> exons;
    private String gene_name;
    private List<Integer> exon_frames;


    public rodRefSeq(String name) {
       super(name);
    }

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId() { return transcript_id; }
    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand() { return strand; }
    /** Returns transcript's full genomic interval (includes all exons with UTRs) */
    public GenomeLoc getLocation() { return transcript_interval; }
    /** Returns genomic interval of the coding sequence (does not include UTRs, but still includes introns, since it's a single interval on the DNA) */
    public GenomeLoc getCodingLocation() { return transcript_coding_interval; }
    /** Name of the gene this transcript corresponds to (NOT gene id such as Entrez etc) */
    public String getGeneName() { return gene_name; }
    /** Number of exons in this transcript */
    public int getNumExons() { return exons.size(); }
    /** Genomic location of the n-th exon; throws an exception if n is out of bounds */
    public GenomeLoc getExonLocation(int n) {
        if ( n >= exons.size() || n < 0 ) throw new StingException("Index out-of-bounds. Transcript has " + exons.size() +" exons; requested: "+n);
        return exons.get(n);
    }
    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<GenomeLoc> getExons() { return exons; }

    /** Returns all exons falling ::entirely:: inside an interval **/
    public List<GenomeLoc> getExonsInInterval( GenomeLoc interval ) {
        List<GenomeLoc> relevantExons = new ArrayList<GenomeLoc>(exons.size());
        for ( GenomeLoc exon : getExons() ) {
            if ( interval.containsP(exon) ) {
                relevantExons.add(exon);
            }
        }

        return relevantExons;
    }

    /** convenience method; returns the numbers of the exons in the interval **/
    public List<Integer> getExonNumbersInInterval( GenomeLoc interval ) {
        List<Integer> numbers = new ArrayList<Integer>();
        int iNo = 0;
        for ( GenomeLoc exon : getExons() ) {
            if ( interval.containsP(exon) ) {
                numbers.add(iNo);
            }
            iNo++;
        }

        return numbers;
    }

    public String getTranscriptUniqueGeneName() {
        return String.format("%s(%s)",getGeneName(),getTranscriptId());
    }

    public String getOverlapString(GenomeLoc position) {
        boolean is_exon = false;
        StringBuilder overlapString = new StringBuilder();
        int exonNo = 1;

        for ( GenomeLoc exon : exons ) {
            if ( exon.containsP(position) ) {
                overlapString.append(String.format("exon_%d",exonNo));
                is_exon = true;
                break;
            }
            exonNo ++;
        }

        if ( ! is_exon ) {
            if ( overlapsCodingP(position) ) {
                overlapString.append("Intron");
            } else {
                overlapString.append("UTR");
            }
        }

        return overlapString.toString();
    }

    /** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */
    public boolean overlapsP (GenomeLoc that) {
        return transcript_interval.overlapsP(that);
    }

    /** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript.
     * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
     * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript.
     * @see #overlapsExonP
     */
    public boolean overlapsCodingP (GenomeLoc that) {
        return transcript_coding_interval.overlapsP(that);
    }

    /** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */
    public boolean overlapsExonP (GenomeLoc that) {
        for ( GenomeLoc e : exons ) {
            if ( e.overlapsP(that) ) return true;
        }
        return false;
    }

    /** Fills this object from a text line in RefSeq (UCSC) text dump file */
    @Override
    public boolean parseLine(final Object header, String[] fields) {
        transcript_id = fields[1];
        if ( fields[3].length()==1 && fields[3].charAt(0)=='+') strand = 1;
        else if ( fields[3].length()==1 && fields[3].charAt(0)=='-') strand = -1;
        else throw new StingException("Expected strand symbol (+/-), found: "+fields[3]);

        String contig_name = fields[2];
        transcript_interval = GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));
        transcript_coding_interval = GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[6])+1, Integer.parseInt(fields[7]));
        gene_name = fields[12];
        String[] exon_starts = fields[9].split(",");
        String[] exon_stops = fields[10].split(",");
        String[] eframes = fields[15].split(",");

        assert exon_starts.length == exon_stops.length : "Data format error: numbers of exon start and stop positions differ";
        assert exon_starts.length == eframes.length : "Data format error: numbers of exons and exon frameshifts differ";

        exons = new ArrayList<GenomeLoc>(exon_starts.length);
        exon_frames = new ArrayList<Integer>(eframes.length);

        for ( int i = 0 ; i < exon_starts.length  ; i++ ) {
            exons.add(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(exon_starts[i])+1, Integer.parseInt(exon_stops[i]) ) );
            exon_frames.add(Integer.decode(eframes[i]));
        }
        return true;
    }

    public String toString() {
        StringBuilder b = new StringBuilder("000\t"); // first field is unused but required in th ecurrent format; just set to something
        b.append(transcript_id);   // #1
        b.append('\t');
        b.append(transcript_interval.getContig()); // #2
        b.append('\t');
        b.append( (strand==1?'+':'-') ); // #3
        b.append('\t');
        b.append( (transcript_interval.getStart() - 1) ); // #4
        b.append('\t');
        b.append( transcript_interval.getStop());  // #5
        b.append('\t');
        b.append( (transcript_coding_interval.getStart() - 1) ); // #6
        b.append('\t');
        b.append( transcript_coding_interval.getStop());  // #7
        b.append('\t');
        b.append(exons.size()); // #8
        b.append('\t');
        for ( GenomeLoc loc : exons ) { b.append( (loc.getStart()-1) ); b.append(','); } // #9
        b.append('\t');
        for ( GenomeLoc loc : exons ) { b.append( loc.getStop() ); b.append(','); } // #10
        b.append("\t0\t"); // # 11 - unused?
        b.append(gene_name); // # 12
        b.append("\tcmpl\tcmpl\t"); // #13, #14 - unused?
        for ( Integer f : exon_frames ) { b.append( f ); b.append(','); } // #15


        return b.toString();
    }

    /** Convenience method, which is packaged here for a lack of better place; it is indeed closely related to
     * rodRefSeq though: takes list of rods (transcripts) overlapping with a given position and determines whether
     * this position is fully whithin an exon of <i>any</i> of those transcripts. Passing null is safe (will return false).
     * NOTE: position can be still within a UTR, see #isCoding
     * @return true if it's an exon
     */
    public static boolean isExon(RODRecordList l) {

        if ( l == null ) return false;

        GenomeLoc loc = l.getLocation();

        for ( GATKFeature t : l ) {
            if ( ((rodRefSeq)t.getUnderlyingObject()).overlapsExonP(loc) ) return true;
        }
        return false;

    }

    /** Convenience method, which is packaged here for a lack of better place; it is indeed closely related to
     * rodRefSeq though: takes list of rods (transcripts) overlapping with a given position and determines whether
     * this position is fully whithin a coding region of <i>any</i> of those transcripts.
     * Passing null is safe (will return false).
     * NOTE: "coding" interval is defined as a single genomic interval, so it
	 * does not include the UTRs of the outermost exons, but it includes introns between exons spliced into a
	 * transcript, or internal exons that are not spliced into a given transcript. To check that a position is
	 * indeed within an exon but not in UTR, use #isCodingExon().
     * @return
     */
    public static boolean isCoding(RODRecordList l) {

        if ( l == null ) return false;

        GenomeLoc loc = l.getLocation();

        for ( GATKFeature t : l ) {
            if ( ((rodRefSeq)t.getUnderlyingObject()).overlapsCodingP(loc) ) return true;
        }
        return false;

    }

    /** Convenience method, which is packaged here for a lack of better place; it is indeed closely related to
     * rodRefSeq though: takes list of rods (transcripts) overlapping with a given position and determines whether
     * this position is fully whithin a coding exon portion (i.e. true coding sequence) of <i>any</i> of those transcripts.
     * Passing null is safe (will return false). In other words, this method returns true if the list contains a transcript,
     * for which the current position is within an exon <i>and</i> within a coding interval simultaneously.
     * @return
     */
    public static boolean isCodingExon(RODRecordList l) {

        if ( l == null ) return false;

        GenomeLoc loc = l.getLocation();

        for ( GATKFeature t : l ) {
            if ( ((rodRefSeq)t.getUnderlyingObject()).overlapsCodingP(loc) && ((rodRefSeq)t.getUnderlyingObject()).overlapsExonP(loc) ) return true;
        }
        return false;

    }

}
