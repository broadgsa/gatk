/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.phasing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

// Streams in VariantContext objects and streams out VariantContexts produced by merging phased segregating polymorphisms into MNP VariantContexts
public class MergePhasedSegregatingPolymorphismsToMNPvcfWriter implements VCFWriter {
    private VCFWriter innerWriter;

    private ReferenceSequenceFile referenceFileForMNPmerging;
    private int maxGenomicDistanceForMNP;

    private VCFRecord vcfrWaitingToMerge;
    private List<VCFRecord> filteredVcfrList;
    private int numMergedRecords;

    private Logger logger;

    private static class VCFRecord {
        public VariantContext vc;
        public byte refBase;

        public VCFRecord(VariantContext vc, byte refBase) {
            this.vc = vc;
            this.refBase = refBase;
        }
    }

    public MergePhasedSegregatingPolymorphismsToMNPvcfWriter(VCFWriter innerWriter, File referenceFile, int maxGenomicDistanceForMNP, Logger logger) {
        this.innerWriter = innerWriter;
        this.referenceFileForMNPmerging = new IndexedFastaSequenceFile(referenceFile);
        this.maxGenomicDistanceForMNP = maxGenomicDistanceForMNP;
        this.vcfrWaitingToMerge = null;
        this.filteredVcfrList = new LinkedList<VCFRecord>();
        this.numMergedRecords = 0;
        this.logger = logger;
    }

    public void writeHeader(VCFHeader header) {
        innerWriter.writeHeader(header);
    }

    public void close() {
        flush();
        innerWriter.close();
    }

    public void flush() {
        stopWaitingToMerge();
        innerWriter.flush();        
    }

    public void add(VariantContext vc, byte refBase) {
        logger.debug("Next VC input = " + VariantContextUtils.getLocation(vc));
        boolean curVcIsNotFiltered = vc.isNotFiltered();

        if (vcfrWaitingToMerge == null) {
            logger.debug("NOT Waiting to merge...");

            if (!filteredVcfrList.isEmpty())
                throw new ReviewedStingException("filteredVcfrList should be empty if not waiting to merge a vc!");

            if (curVcIsNotFiltered) { // still need to wait before can release vc
                logger.debug("Waiting for new variant " + VariantContextUtils.getLocation(vc));
                vcfrWaitingToMerge = new VCFRecord(vc, refBase);
            }
            else {
                logger.debug("DIRECTLY output " + VariantContextUtils.getLocation(vc));
                innerWriter.add(vc, refBase);
            }
        }
        else { // waiting to merge vcfrWaitingToMerge
            logger.debug("Waiting to merge " + VariantContextUtils.getLocation(vcfrWaitingToMerge.vc));

            if (!curVcIsNotFiltered) {
                logger.debug("Caching unprocessed output " + VariantContextUtils.getLocation(vc));
                filteredVcfrList.add(new VCFRecord(vc, refBase));
            }
            else { // waiting to merge vcfrWaitingToMerge, and curVcIsNotFiltered. So, attempt to merge them:
                boolean mergedRecords = false;
                if (minDistance(vcfrWaitingToMerge.vc, vc) <= maxGenomicDistanceForMNP) {
                    VariantContext mergedVc = VariantContextUtils.mergeIntoMNP(vcfrWaitingToMerge.vc, vc, referenceFileForMNPmerging);
                    if (mergedVc != null) {
                        mergedRecords = true;
                        vcfrWaitingToMerge = new VCFRecord(mergedVc, vcfrWaitingToMerge.refBase);
                        numMergedRecords++;
                    }
                }
                if (!mergedRecords) {
                    stopWaitingToMerge();
                    vcfrWaitingToMerge = new VCFRecord(vc, refBase);
                }
                logger.debug("Merged? = " + mergedRecords);
            }
        }
    }

    private void stopWaitingToMerge() {
        if (vcfrWaitingToMerge == null) {
            if (!filteredVcfrList.isEmpty())
                throw new ReviewedStingException("filteredVcfrList should be empty if not waiting to merge a vc!");
            return;
        }

        innerWriter.add(vcfrWaitingToMerge.vc, vcfrWaitingToMerge.refBase);
        vcfrWaitingToMerge = null;

        for (VCFRecord vcfr : filteredVcfrList)
            innerWriter.add(vcfr.vc, vcfr.refBase);
        filteredVcfrList.clear();
    }

    public int getNumMergedRecords() {
        return numMergedRecords;
    }

    public static int minDistance(VariantContext vc1, VariantContext vc2) {
        return VariantContextUtils.getLocation(vc1).minDistance(VariantContextUtils.getLocation(vc2));
    }

    /**
     * Gets a string representation of this object.
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
    }
}