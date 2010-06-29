package org.broadinstitute.sting.oneoffprojects.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluator;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.analysis.AminoAcid;
import org.broadinstitute.sting.utils.analysis.AminoAcidTable;
import org.broadinstitute.sting.utils.analysis.AminoAcidUtils;

import java.util.ArrayList;

/*
 * Copyright (c) 2010 The Broad Institute
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
 * @author chartl
 * @since June 28, 2010
 */

@Analysis(name = "Amino Acid Transition", description = "Calculates the Transition Matrix for coding variants; entries are Total, Num. Ti, Num. Tv, Ratio")
public class AminoAcidTransition extends VariantEvaluator {

    ////////////////////////////////////////////////////////////
    ////            INTERNAL DATA POINT CLASSES
    ////////////////////////////////////////////////////////////

    // a mapping from amino acid transition score histogram bin to Ti/Tv ratio
    @DataPoint(name="Amino Acid Table", description = "TiTv counts by amino acid change")
    AminoAcidTiTvTable acidTable = null;

    class TiTvCount {
        public int ti;
        public int tv;

        public TiTvCount() {
            ti = 0;
            tv = 0;
        }

        public int getTotal() {
            return ti + tv;
        }

        public double getRatio() {
            return ( (double) ti )/tv;
        }

        public String toString() {
            return String.format("%d:%d:%d:%f",getTotal(),ti,tv,getRatio());
        }
    }

    class AminoAcidTiTvTable implements TableType {

        private TiTvCount[][] countsByAAChange;

        public AminoAcidTiTvTable() {
            countsByAAChange = new TiTvCount[AminoAcid.values().length][AminoAcid.values().length];
            for ( int i = 0; i < AminoAcid.values().length; i ++ ) {
                for ( int j = 0; j < AminoAcid.values().length; j++ ) {
                    countsByAAChange[i][j] = new TiTvCount();
                }
            }
        }

        public Object[] getRowKeys() {
            return AminoAcidUtils.getAminoAcidCodes();

        }

        public Object[] getColumnKeys() {
            return AminoAcidUtils.getAminoAcidCodes();
        }

        public TiTvCount getCell(int x, int y) {
            return countsByAAChange[x][y];
        }

        public String getName() {
            return "AminoAcidTransitionTable";
        }

        public void update(AminoAcid reference, AminoAcid alternate, boolean isTransition) {
            TiTvCount counter = countsByAAChange[reference.ordinal()][alternate.ordinal()];
            if ( isTransition ) {
                counter.ti++;
            } else {
                counter.tv++;
            }
        }
    }

    ////////////////////////////////////////////////////////////
    ////        CORE VARIANT EVALUATOR DATA AND METHODS
    ////////////////////////////////////////////////////////////

    private String infoKey;
    private String infoValueSplit;
    private boolean useCodons;
    private AminoAcidTable lookup;

    public AminoAcidTransition(VariantEvalWalker parent) {
        super(parent);
        getParsingInformation(parent);
        lookup = new AminoAcidTable();
    }

    private void getParsingInformation(VariantEvalWalker parent) {
        // todo -- allow me to be flexible
        infoKey = "something";
        infoValueSplit = "somethingElse";
        useCodons = false;
    }

    public String getName() {
        return "AminoAcidTransitionTable";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public boolean enabled() {
        return true;
    }

    public String toString() {
        return getName();
    }

    public String update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String interesting = null;
        if ( eval != null && eval.hasAttribute(infoKey) ) {
            String[] parsedNames = ( (String) eval.getAttribute(infoKey)).split(infoValueSplit);
            String first = parsedNames [0];
            String second = parsedNames [1];
            AminoAcid reference;
            AminoAcid alternate;
            if ( useCodons ) {
                reference = lookup.getEukaryoticAA(first);
                alternate = lookup.getEukaryoticAA(second);
            } else {
                reference = lookup.getAminoAcidByCode(first);
                alternate = lookup.getAminoAcidByCode(second);
            }

            if ( reference == null ) {
                interesting = "Unknown Reference Codon";
            } else if ( alternate == null ) {
                interesting = "Unknown Alternate Codon";
            } else {
                acidTable.update(reference,alternate,eval.isTransition());
            }

        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    //public void finalizeEvaluation() {
    //
    //}
}