package org.broadinstitute.sting.gatk.refdata.features.refseq;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.ArrayList;

/**
 * the ref seq codec
 */
public class RefSeqCodec implements FeatureCodec {

    @Override
    public Feature decodeLoc(String line) {
        if (line.startsWith("#")) return null;
        String fields[] = line.split("\t");
        if (fields.length < 3) throw new TribbleException("RefSeq (decodeLoc) : Unable to parse line -> " + line + ", we expected at least 3 columns, we saw " + fields.length);
        String contig_name = fields[2];
        return new RefSeqFeature(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));
    }

    /** Fills this object from a text line in RefSeq (UCSC) text dump file */
    @Override
    public Feature decode(String line) {
        if (line.startsWith("#")) return null;
        String fields[] = line.split("\t");

        // we reference postion 15 in the split array below, make sure we have at least that many columns
        if (fields.length < 16) throw new TribbleException("RefSeq (decode) : Unable to parse line -> " + line + ", we expected at least 16 columns, we saw " + fields.length);
        String contig_name = fields[2];
        RefSeqFeature feature = new RefSeqFeature(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));

        feature.setTranscript_id(fields[1]);
        if ( fields[3].length()==1 && fields[3].charAt(0)=='+') feature.setStrand(1);
        else if ( fields[3].length()==1 && fields[3].charAt(0)=='-') feature.setStrand(-1);
        else throw new UserException.MalformedFile("Expected strand symbol (+/-), found: "+fields[3] + " for line=" + line);


        feature.setTranscript_interval(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5])));
        feature.setTranscript_coding_interval(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[6])+1, Integer.parseInt(fields[7])));
        feature.setGene_name(fields[12]);
        String[] exon_starts = fields[9].split(",");
        String[] exon_stops = fields[10].split(",");
        String[] eframes = fields[15].split(",");

        if ( exon_starts.length == exon_stops.length )
            throw new UserException.MalformedFile("Data format error: numbers of exon start and stop positions differ for line=" + line);
        if ( exon_starts.length == eframes.length )
            throw new UserException.MalformedFile("Data format error: numbers of exons and exon frameshifts differ for line=" + line);

        ArrayList<GenomeLoc> exons = new ArrayList<GenomeLoc>(exon_starts.length);
        ArrayList<Integer> exon_frames = new ArrayList<Integer>(eframes.length);

        for ( int i = 0 ; i < exon_starts.length  ; i++ ) {
            exons.add(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(exon_starts[i])+1, Integer.parseInt(exon_stops[i]) ) );
            exon_frames.add(Integer.decode(eframes[i]));
        }

        feature.setExons(exons);
        feature.setExon_frames(exon_frames);
        return feature;
    }

    @Override
    public Object readHeader(LineReader reader) {
        return null;
    }

    @Override
    public Class getFeatureType() {
        return RefSeqCodec.class;
    }
}
