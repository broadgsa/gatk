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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.sting.gatk.walkers.annotator.SnpEff.SnpEffEffect;

public class SnpEffUnitTest {

    @Test
    public void testParseWellFormedEffect() {
        String effectName = "NON_SYNONYMOUS_CODING";
        String[] effectMetadata = { "MODERATE", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990", "exon_1_69037_69829" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertTrue( effect.isWellFormed() && effect.isCoding() );
    }

    @Test
    public void testParseInvalidEffectNameEffect() {
        String effectName = "MADE_UP_EFFECT";
        String[] effectMetadata = { "MODERATE", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990", "exon_1_69037_69829" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertFalse(effect.isWellFormed());
    }

    @Test
    public void testParseInvalidEffectImpactEffect() {
        String effectName = "NON_SYNONYMOUS_CODING";
        String[] effectMetadata = { "MEDIUM", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990", "exon_1_69037_69829" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertFalse(effect.isWellFormed());
    }

    @Test
    public void testParseWrongNumberOfMetadataFieldsEffect() {
        String effectName = "NON_SYNONYMOUS_CODING";
        String[] effectMetadata = { "MODERATE", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertFalse(effect.isWellFormed());
    }

    @Test
    public void testParseSnpEffWarningEffect() {
        String effectName = "NON_SYNONYMOUS_CODING";
        String[] effectMetadata = { "MODERATE", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990", "exon_1_69037_69829", "SNPEFF_WARNING" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertTrue( ! effect.isWellFormed() && effect.getParseError().equals("SnpEff issued the following warning: SNPEFF_WARNING") );
    }

    @Test
    public void testParseSnpEffErrorEffect() {
        String effectName = "NON_SYNONYMOUS_CODING";
        String[] effectMetadata = { "MODERATE", "Aca/Gca", "T/A", "OR4F5", "protein_coding", "CODING", "ENST00000534990", "exon_1_69037_69829", "", "SNPEFF_ERROR" };

        SnpEffEffect effect = new SnpEffEffect(effectName, effectMetadata);
        Assert.assertTrue( ! effect.isWellFormed() && effect.getParseError().equals("SnpEff issued the following error: SNPEFF_ERROR") );
    }
}
