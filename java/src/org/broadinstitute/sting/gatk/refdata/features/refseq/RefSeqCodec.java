package org.broadinstitute.sting.gatk.refdata.features.refseq;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.List;

/**
 * the ref seq codec
 */
public class RefSeqCodec implements FeatureCodec {

    @Override
    public Feature decodeLoc(String line) {
        String fields[] = line.split("\t");
        String contig_name = fields[2];
        return new RefSeqFeature(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));
    }

    /** Fills this object from a text line in RefSeq (UCSC) text dump file */
    @Override
    public Feature decode(String line) {
        String fields[] = line.split("\t");

        String contig_name = fields[2];
        RefSeqFeature feature = new RefSeqFeature(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));

        feature.setTranscript_id(fields[1]);
        if ( fields[3].length()==1 && fields[3].charAt(0)=='+') feature.setStrand(1);
        else if ( fields[3].length()==1 && fields[3].charAt(0)=='-') feature.setStrand(-1);
        else throw new StingException("Expected strand symbol (+/-), found: "+fields[3]);


        feature.setTranscript_interval(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5])));
        feature.setTranscript_coding_interval(GenomeLocParser.parseGenomeLoc(contig_name, Integer.parseInt(fields[6])+1, Integer.parseInt(fields[7])));
        feature.setGene_name(fields[12]);
        String[] exon_starts = fields[9].split(",");
        String[] exon_stops = fields[10].split(",");
        String[] eframes = fields[15].split(",");

        assert exon_starts.length == exon_stops.length : "Data format error: numbers of exon start and stop positions differ";
        assert exon_starts.length == eframes.length : "Data format error: numbers of exons and exon frameshifts differ";

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
