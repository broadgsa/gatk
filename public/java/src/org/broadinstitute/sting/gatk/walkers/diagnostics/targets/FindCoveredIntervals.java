/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;

import java.io.PrintStream;

@PartitionBy(PartitionType.CONTIG)
@ActiveRegionExtension(extension = 0, maxRegion = 50000)
public class FindCoveredIntervals extends ActiveRegionWalker<GenomeLoc, Long> {
    @Output(required = true)
    private PrintStream out;

    @Override
    // Look to see if the region has sufficient coverage
    public double isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {

        int depth = ThresHolder.DEFAULTS.getFilteredCoverage(context.getBasePileup());

        // note the linear probability scale
        int coverageThreshold = 20;
        return Math.min((double) depth / coverageThreshold, 1);

    }

    @Override
    public GenomeLoc map(final ActiveRegion activeRegion, final RefMetaDataTracker tracker) {
        if (activeRegion.isActive)
            return activeRegion.getLocation();
        else
            return null;
    }

    @Override
    public Long reduceInit() {
        return 0L;
    }

    @Override
    public Long reduce(final GenomeLoc value, Long reduce) {
        if (value != null) {
            out.println(value.toString());
            return reduce++;
        } else
            return reduce;
    }

    @Override
    public void onTraversalDone(Long reduce) {
        logger.info(String.format("Found %d intervals", reduce));
    }
}
