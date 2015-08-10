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

package org.broadinstitute.gatk.tools.walkers.annotator.interfaces;

import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;

import java.util.List;
import java.util.Set;

@DocumentedGATKFeature(enable = true, groupName = HelpConstants.DOCS_CAT_ANNOT, summary = "Annotations available to VariantAnnotator and the variant callers (some restrictions apply)")
public abstract class VariantAnnotatorAnnotation {
    // return the INFO keys
    public abstract List<String> getKeyNames();

    // initialization method (optional for subclasses, and therefore non-abstract)
    public void initialize ( AnnotatorCompatible walker, GenomeAnalysisEngine toolkit, Set<VCFHeaderLine> headerLines ) { }
}