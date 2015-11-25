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

package org.broadinstitute.gatk.engine.datasources.providers;

import org.testng.Assert;
import org.broadinstitute.gatk.utils.GenomeLoc;

import org.testng.annotations.Test;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 1:04:27 PM
 *
 */

/**
 * Test reading the reference for a given read.
 */

public class ReadReferenceViewUnitTest extends ReferenceViewTemplate {


    /**
     * tests that the ReadReferenceView correctly generates X's when a read overhangs the
     * end of a contig
     */
    @Test
    public void testOverhangingRead() {
        testOverhangingGivenSize(25,0);
        testOverhangingGivenSize(25,12);
        testOverhangingGivenSize(25,24);
    }


    /**
     * a private method, that tests getting the read sequence for reads that overlap the end of the
     * contig
     * @param readLength the length of the read
     * @param overlap the amount of overlap
     */
    private void testOverhangingGivenSize(int readLength, int overlap) {
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(sequenceFile.getSequenceDictionary().getSequences().size()-1);
        final long contigStart = selectedContig.getSequenceLength() - (readLength - overlap - 1);
        final long contigStop = selectedContig.getSequenceLength() + overlap;

        ReadShardDataProvider dataProvider = new ReadShardDataProvider(null,genomeLocParser,null,sequenceFile,null);
        ReadReferenceView view = new ReadReferenceView(dataProvider);

        SAMRecord rec = buildSAMRecord(selectedContig.getSequenceName(),(int)contigStart,(int)contigStop);
        ReferenceSequence expectedAsSeq = sequenceFile.getSubsequenceAt(selectedContig.getSequenceName(),(int)contigStart,selectedContig.getSequenceLength());
        //char[] expected = StringUtil.bytesToString(expectedAsSeq.getBases()).toCharArray();
        byte[] expected = expectedAsSeq.getBases();
        byte[] actual = view.getReferenceBases(rec);

        Assert.assertEquals((readLength - overlap), expected.length);
        Assert.assertEquals(readLength, actual.length);
        int xRange = 0;
        for (; xRange < (readLength - overlap); xRange++) {
            Assert.assertTrue(actual[xRange] != 'X');
        }
        for (; xRange < actual.length; xRange++) {
            Assert.assertTrue(actual[xRange] == 'X');
        }
    }


    /**
     * Compares the contents of the fasta and view at a specified location.
     * @param loc the location to validate
     */
    protected void validateLocation( GenomeLoc loc ) {
        SAMRecord read = buildSAMRecord( loc.getContig(), (int)loc.getStart(), (int)loc.getStop() );

        ReadShardDataProvider dataProvider = new ReadShardDataProvider(null,genomeLocParser,null,sequenceFile,null);
        ReadReferenceView view = new ReadReferenceView(dataProvider);

        ReferenceSequence expectedAsSeq = sequenceFile.getSubsequenceAt(loc.getContig(),loc.getStart(),loc.getStop());
        byte[] expected = expectedAsSeq.getBases();
        byte[] actual = view.getReferenceBases(read);

        org.testng.Assert.assertEquals(actual,expected,String.format("Base array at  in shard %s does not match expected",loc.toString()));
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
