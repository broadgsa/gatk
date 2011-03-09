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

package org.broadinstitute.sting.oneoffprojects.walkers.CNV;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Walks along reference and calculates the genes (from "refseq" ROD) for each interval defined in "intervals" ROD.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = GeneNamesIntervalWalker.REFSEQ_ROD_NAME, type = AnnotatorInputTableFeature.class)})

public class GeneNamesIntervalWalker extends RodWalker<GeneNames, GeneNames> {
    @Output
    protected PrintStream out;

    public final static String REFSEQ_ROD_NAME = "refseq";

    public final static String REFSEQ_NAME2 = "name2";


    public boolean isReduceByInterval() {
        return true;
    }

    public void initialize() {
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public GeneNames reduceInit() {
        return new GeneNames();
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public GeneNames map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        return new GeneNames().addGenes(tracker.getReferenceMetaData(REFSEQ_ROD_NAME));
    }

    public GeneNames reduce(GeneNames add, GeneNames runningCount) {
        if (add == null)
            add = new GeneNames();

        return runningCount.addIn(add);
    }

    /**
     * @param results the genes found in each interval.
     */
    public void onTraversalDone(List<Pair<GenomeLoc, GeneNames>> results) {
        for (Pair<GenomeLoc, GeneNames> result : results ) {
            GenomeLoc loc = result.getFirst();
            GeneNames names = result.getSecond();

            out.println(loc + "\t" + names);
        }
    }
}

class GeneNames {
    public Set<String> geneNames;

    public GeneNames() {
        this.geneNames = new HashSet<String>();
    }

    public GeneNames addIn(GeneNames other) {
        this.geneNames.addAll(other.geneNames);

        return this;
    }

    public GeneNames addGenes(List<Object> refSeqRODs) {
        for (Object refSeqObject : refSeqRODs) {
            AnnotatorInputTableFeature refSeqAnnotation = (AnnotatorInputTableFeature) refSeqObject;
            if (refSeqAnnotation.containsColumnName(GeneNamesIntervalWalker.REFSEQ_NAME2))
                geneNames.add(refSeqAnnotation.getColumnValue(GeneNamesIntervalWalker.REFSEQ_NAME2));
        }

        return this;
    }

    public String toString() {
        if (geneNames.isEmpty())
            return ".";

        StringBuilder sb = new StringBuilder();

        for (String gene : geneNames)
            sb.append(gene).append(";");

        return sb.toString();
    }
}

