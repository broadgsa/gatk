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

package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.PrintStream;
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jun 4, 2009
 * Time: 4:38:19 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ViolationVariantAnalysis extends BasicVariantAnalysis {
    PrintStream violationsOut = null;

    public ViolationVariantAnalysis(final String name) {
        super(name);
    }

    public void initialize(VariantEvalWalker master, PrintStream out, String filename) {
        super.initialize(master, out, filename);

        try {
            violationsOut = filename == null ? out : new PrintStream(new File(filename + ".violations"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    public PrintStream getCallPrintStream() {
        return violationsOut;
    }
}