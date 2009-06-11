package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.junit.Assert;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 1:04:27 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test reading the reference for a given read.
 */

public class ReadReferenceViewTest extends ReferenceViewTemplate {
    /**
     * Compares the contents of the fasta and view at a specified location.
     * @param loc
     */
    protected void validateLocation( GenomeLoc loc ) {
        SAMRecord read = buildSAMRecord( loc.getContig(), (int)loc.getStart(), (int)loc.getStop() );

        ShardDataProvider dataProvider = new ShardDataProvider(null,null,sequenceFile,null);
        ReadReferenceView view = new ReadReferenceView(dataProvider);

        ReferenceSequence expectedAsSeq = sequenceFile.getSubsequenceAt(loc.getContig(),loc.getStart(),loc.getStop());
        char[] expected = StringUtil.bytesToString(expectedAsSeq.getBases()).toCharArray();
        char[] actual = view.getReferenceBases(read);

        Assert.assertArrayEquals(String.format("Base array at  in shard %s does not match expected",loc.toString()),
                                 expected,
                                 actual);
    }


    /**
     * Build a SAM record featuring the absolute minimum required dataset.
     * TODO: Blatantly copied from LocusViewTemplate.  Refactor these into a set of tools.
     * @param contig Contig to populate.
     * @param alignmentStart start of alignment
     * @param alignmentEnd end of alignment
     * @return New SAM Record
     */
    protected SAMRecord buildSAMRecord( String contig, int alignmentStart, int alignmentEnd ) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(sequenceFile.getSequenceDictionary());

        SAMRecord record = new SAMRecord(header);

        record.setReferenceIndex(sequenceFile.getSequenceDictionary().getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(alignmentEnd-alignmentStart+1, CigarOperator.M));
        record.setCigar(cigar);
        return record;
    }


}
