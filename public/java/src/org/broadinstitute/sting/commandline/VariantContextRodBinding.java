/*
* Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Collection;

/**
* A RodBinding representing a walker argument that gets bound to a ROD track containing VariantContexts
*/
public class VariantContextRodBinding extends RodBinding {
    /**
     * Create a new RodBinding specialized to provide VariantContexts.
     * @param variableName the name of the field in the walker that we will bind the ROD track too
     * @param sourceFile the data source from which we will read the VCs
     * @param parser the Engine parser used to obtain information about this argument, such as its underlying file type
     */
    protected VariantContextRodBinding(final String variableName, final String sourceFile, final ParsingEngine parser) {
        super(VariantContext.class, variableName, sourceFile, parser);
    }

    /**
     * Forwarding method to identical tracker method
     */
    public Collection<VariantContext> getVariantContexts(final RefMetaDataTracker tracker,
                                                         final GenomeLoc curLocation,
                                                         final boolean requireStartHere,
                                                         final boolean takeFirstOnly ) {
        return tracker.getVariantContexts(variableName, curLocation, requireStartHere, takeFirstOnly);
    }

    /**
     * Forwarding method to identical tracker method
     * @param tracker
     * @param curLocation
     * @param requireStartHere
     * @return
     */
    public VariantContext getVariantContext(final RefMetaDataTracker tracker,
                                            final GenomeLoc curLocation,
                                            final boolean requireStartHere ) {
        return tracker.getVariantContext(variableName, curLocation, requireStartHere);
    }

    /**
     * Forwarding method to identical tracker method
     */
    public VariantContext getVariantContext(final RefMetaDataTracker tracker,
                                            final GenomeLoc curLocation) {
        return tracker.getVariantContext(variableName, curLocation);
    }
}
