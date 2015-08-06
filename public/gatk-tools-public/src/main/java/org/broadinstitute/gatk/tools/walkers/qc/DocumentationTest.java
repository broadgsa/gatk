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

package org.broadinstitute.gatk.tools.walkers.qc;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Summary test
 *
 * <p>Body test</p>
 */
@Hidden
public class DocumentationTest extends RodWalker<Integer, Integer> {
    // the docs for the arguments are in the collection
    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * dbSNP comparison VCF.  By default, the dbSNP file is used to specify the set of "known" variants.
     * Other sets can be specified with the -knownName (--known_names) argument.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * detailed documentation about the argument goes here.
     */
    @Input(fullName="listofRodBinding", shortName = "disc", doc="Output variants that were not called in this Feature comparison track", required=false)
    private List<RodBinding<VariantContext>> listOfRodBinding = Collections.emptyList();

    @Input(fullName="optionalRodBinding", shortName = "conc", doc="Output variants that were also called in this Feature comparison track", required=false)
    private RodBinding<VariantContext> concordanceTrack;

    @Input(fullName="optionalRodBindingWithoutDefault", shortName = "optionalRodBindingWithoutDefault", doc="Output variants that were also called in this Feature comparison track", required=false)
    private RodBinding<VariantContext> noDefaultOptionalRodBinding;

    @Input(fullName="optionalRodBindingWithoutDefaultNull", shortName = "shortTest", doc="Output variants that were also called in this Feature comparison track", required=false)
    private RodBinding<VariantContext> noDefaultOptionalRodBindingNull = null;

    @Input(fullName="featureArg", shortName = "featureArg", doc="A RodBinding of feature", required=false)
    private RodBinding<Feature> featureArg = null;

    @Output(doc="VCFWriter")
    protected VariantContextWriter vcfWriter = null;

    @Advanced
    @Argument(fullName="setString", shortName="sn", doc="Sample name to be included in the analysis. Can be specified multiple times.", required=false)
    public Set<String> sampleNames;

    @Argument(fullName="setStringInitialized", shortName="setStringInitialized", doc="Sample name to be included in the analysis. Can be specified multiple times.", required=false)
    public Set<String> setStringInitialized = new HashSet<String>();

    @Argument(shortName="optionalArgWithMissinglessDefault", doc="One or more criteria to use when selecting the data.  Evaluated *after* the specified samples are extracted and the INFO-field annotations are updated.", required=false)
    public ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<String>();

    @Argument(shortName="AAAAA", fullName = "AAAAA", doc="Should be the first argument", required=false)
    public boolean FIRST_ARG = false;

    @Advanced
    @Argument(fullName="booleanArg", shortName="env", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean EXCLUDE_NON_VARIANTS = false;

    @Advanced
    @Argument(fullName="booleanArray", shortName="booleanArray", doc="x", required=false)
    private boolean[] boolArray = null;

    @Argument(fullName="enumTest", shortName="enumTest", doc="Test enum", required=false)
    private TestEnum TestEnumArg = TestEnum.ENUM2;
    public enum TestEnum {
        /** Docs for enum1 */
        ENUM1,
        /** Docs for enum2 */
        ENUM2
    }

    @Hidden
    @Argument(fullName="hiddenArg", shortName="keepAF", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean KEEP_AF_SPECTRUM = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) { return value + sum; }
    public void onTraversalDone(Integer result) { }
}
