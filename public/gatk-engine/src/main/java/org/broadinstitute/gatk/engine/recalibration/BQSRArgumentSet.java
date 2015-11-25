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

package org.broadinstitute.gatk.engine.recalibration;

import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;

import java.io.File;
import java.util.List;

public class BQSRArgumentSet {
    // declare public, STL-style for easier and more efficient access:
    private File BQSR_RECAL_FILE;
    private int quantizationLevels;
    private List<Integer> staticQuantizedQuals;
    private boolean roundDown;
    private boolean disableIndelQuals;
    private boolean emitOriginalQuals;
    private int PRESERVE_QSCORES_LESS_THAN;
    private double globalQScorePrior;

    public BQSRArgumentSet(final GATKArgumentCollection args) {
        this.BQSR_RECAL_FILE = args.BQSR_RECAL_FILE;
        this.quantizationLevels = args.quantizationLevels;
        this.staticQuantizedQuals = args.staticQuantizationQuals;
        this.roundDown = args.roundDown;
        this.disableIndelQuals = args.disableIndelQuals;
        this.emitOriginalQuals = args.emitOriginalQuals;
        this.PRESERVE_QSCORES_LESS_THAN = args.PRESERVE_QSCORES_LESS_THAN;
        this.globalQScorePrior = args.globalQScorePrior;
    }

    public File getRecalFile() { return BQSR_RECAL_FILE; }

    public int getQuantizationLevels() { return quantizationLevels; }

    public List<Integer> getStaticQuantizedQuals() {return staticQuantizedQuals; }

    public boolean getRoundDown() {return roundDown; }

    public boolean shouldDisableIndelQuals() { return disableIndelQuals; }

    public boolean shouldEmitOriginalQuals() { return emitOriginalQuals; }

    public int getPreserveQscoresLessThan() { return PRESERVE_QSCORES_LESS_THAN; }

    public double getGlobalQScorePrior() { return globalQScorePrior; }

    public void setRecalFile(final File BQSR_RECAL_FILE) {
        this.BQSR_RECAL_FILE = BQSR_RECAL_FILE;
    }

    public void setQuantizationLevels(final int quantizationLevels) {
        this.quantizationLevels = quantizationLevels;
    }

    public void setStaticQuantizedQuals(final List<Integer> staticQuantizedQuals) { this.staticQuantizedQuals = staticQuantizedQuals; }

    public void setRoundDown(final boolean roundDown) {
        this.roundDown = roundDown;
    }

    public void setDisableIndelQuals(final boolean disableIndelQuals) {
        this.disableIndelQuals = disableIndelQuals;
    }

    public void setEmitOriginalQuals(final boolean emitOriginalQuals) {
        this.emitOriginalQuals = emitOriginalQuals;
    }

    public void setPreserveQscoresLessThan(final int PRESERVE_QSCORES_LESS_THAN) {
        this.PRESERVE_QSCORES_LESS_THAN = PRESERVE_QSCORES_LESS_THAN;
    }

    public void setGlobalQScorePrior(final double globalQScorePrior) {
        this.globalQScorePrior = globalQScorePrior;
    }
}
