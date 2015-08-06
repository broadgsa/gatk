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

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.iterators.RNAReadTransformer;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

/**
 * Reduce NDN cigar elements to one N element.
 *
 * <p>This read transformer will refactor cigar strings that contain N-D-N elements to one N element (with total length
 * of the three refactored elements). The engine parameter that activate this read transformer is
 * `--refactor_NDN_cigar_string` / `-fixNDN`</p>
 *
 * <h3>Rationale</h3>
 * <p>Some RNAseq aligners that use a known transcriptome resource (such as TopHat2) produce NDN elements in read CIGARS
 * when a small exon is entirely deleted during transcription, which ends up looking like [exon1]NDN[exon3]. Currently
 * we consider that the internal N-D-N motif is illegal and we error out when we encounter it. By refactoring the cigar string of
 * those specific reads, this read transformer allows users of TopHat and other tools to circumvent this problem without
 * affecting the rest of their dataset. From the point of view of variant calling, there is no meaningful difference between
 * the two representations.</p>
 *
 * <h3>Developer notes</h3>
 * <ul>
 *     <li>Any walker that needs to apply this functionality should apply that read transformer in its map function, since it won't be activated by the GATK engine.</li>
 * </ul>
 *
 *
 * @author ami
 * @since 04/22/14
 */

public class NDNCigarReadTransformer extends RNAReadTransformer {

    private boolean refactorReads;

    @Override
    public ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker) {
        refactorReads = engine.getArguments().REFACTOR_NDN_CIGAR_READS;

        return ApplicationTime.HANDLED_IN_WALKER;   //  NOTE: any walker that need that functionality should apply that read transformer in its map function, since it won't be activated by the GATK engine.
    }

    @Override
    public GATKSAMRecord apply(final GATKSAMRecord read) {
        if(read == null)
            throw new UserException.BadInput("try to transform a null GATKSAMRecord");
        final Cigar originalCigar = read.getCigar();
        if (originalCigar.isValid(read.getReadName(),-1) != null)
            throw new UserException.BadInput("try to transform a read with non-valid cigar string: readName: "+read.getReadName()+" Cigar String: "+originalCigar);
        read.setCigar(refactorNDNtoN(originalCigar));
        return read;
    }

    @Override
    public boolean enabled() {
        return refactorReads;
    }



    protected Cigar refactorNDNtoN(final Cigar originalCigar) {
        final Cigar refactoredCigar = new Cigar();
        final int cigarLength = originalCigar.numCigarElements();
        for(int i = 0; i < cigarLength; i++){
            final CigarElement element = originalCigar.getCigarElement(i);
            if(element.getOperator() == CigarOperator.N && thereAreAtLeast2MoreElements(i,cigarLength)){
                final CigarElement nextElement = originalCigar.getCigarElement(i+1);
                final CigarElement nextNextElement = originalCigar.getCigarElement(i+2);

                // if it is N-D-N replace with N (with the total length) otherwise just add the first N.
                if(nextElement.getOperator() == CigarOperator.D && nextNextElement.getOperator() == CigarOperator.N){
                    final int threeElementsLength = element.getLength() + nextElement.getLength() + nextNextElement.getLength();
                    final CigarElement refactoredElement = new CigarElement(threeElementsLength,CigarOperator.N);
                    refactoredCigar.add(refactoredElement);
                    i += 2; //skip the elements that were refactored
                }
                else
                    refactoredCigar.add(element);  // add only the first N
            }
            else
                refactoredCigar.add(element);  // add any non-N element
        }
        return refactoredCigar;
    }

    private boolean thereAreAtLeast2MoreElements(final int index, final int cigarLength){
        return index < cigarLength - 2;
    }

}

