/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broad.tribble.util.variantcontext.VariantContext;

/**
 * Created by IntelliJ IDEA.
 * User: depristo, ebanks
 * Date: Jan 22, 2010
 * Time: 2:25:19 PM
 *
 * Useful helper class to communicate the results of calculateGenotype to framework
 */
public class VariantCallContext {
    public VariantContext vc = null;
    public byte refBase;

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP) {
        this.vc = vc;
        this.confidentlyCalled = confidentlyCalledP;
    }

    VariantCallContext(VariantContext vc, byte ref, boolean confidentlyCalledP) {
        this.vc = vc;
        this.refBase = ref;
        this.confidentlyCalled = confidentlyCalledP;
    }

    // blank variant context => we're a ref site
    VariantCallContext(boolean confidentlyCalledP) {
        this.confidentlyCalled = confidentlyCalledP;
    }

    public void setRefBase(byte ref) {
        this.refBase = ref;
    }
}