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

package org.broadinstitute.sting.utils.codecs.snpEff;

/**
 * A set of constants associated with the SnpEff codec.
 *
 * @author David Roazen
 */
public class SnpEffConstants {

    // Possible SnpEff biological effects and their associated impacts:
    public enum EffectType {
        START_GAINED                      (EffectImpact.HIGH),
        START_LOST                        (EffectImpact.HIGH),
        EXON_DELETED                      (EffectImpact.HIGH),
        FRAME_SHIFT                       (EffectImpact.HIGH),
        STOP_GAINED                       (EffectImpact.HIGH),
        STOP_LOST                         (EffectImpact.HIGH),
        SPLICE_SITE_ACCEPTOR              (EffectImpact.HIGH),
        SPLICE_SITE_DONOR                 (EffectImpact.HIGH),

        NON_SYNONYMOUS_CODING             (EffectImpact.MODERATE),
        UTR_5_DELETED                     (EffectImpact.MODERATE),
        UTR_3_DELETED                     (EffectImpact.MODERATE),
        CODON_INSERTION                   (EffectImpact.MODERATE),
        CODON_CHANGE_PLUS_CODON_INSERTION (EffectImpact.MODERATE),
        CODON_DELETION                    (EffectImpact.MODERATE),
        CODON_CHANGE_PLUS_CODON_DELETION  (EffectImpact.MODERATE),

        NONE                              (EffectImpact.LOW),
        CHROMOSOME                        (EffectImpact.LOW),
        INTERGENIC                        (EffectImpact.LOW),
        UPSTREAM                          (EffectImpact.LOW),
        UTR_5_PRIME                       (EffectImpact.LOW),
        SYNONYMOUS_START                  (EffectImpact.LOW),
        NON_SYNONYMOUS_START              (EffectImpact.LOW),
        CDS                               (EffectImpact.LOW),
        GENE                              (EffectImpact.LOW),
        TRANSCRIPT                        (EffectImpact.LOW),
        EXON                              (EffectImpact.LOW),
        SYNONYMOUS_CODING                 (EffectImpact.LOW),
        CODON_CHANGE                      (EffectImpact.LOW),
        SYNONYMOUS_STOP                   (EffectImpact.LOW),
        NON_SYNONYMOUS_STOP               (EffectImpact.LOW),
        INTRON                            (EffectImpact.LOW),
        UTR_3_PRIME                       (EffectImpact.LOW),
        DOWNSTREAM                        (EffectImpact.LOW),
        INTRON_CONSERVED                  (EffectImpact.LOW),
        INTERGENIC_CONSERVED              (EffectImpact.LOW),
        CUSTOM                            (EffectImpact.LOW);

        private final EffectImpact impact;

        EffectType ( EffectImpact impact ) {
            this.impact = impact;
        }

        public EffectImpact getImpact() {
            return impact;
        }
    }

    public enum EffectImpact {
        LOW       (1),
        MODERATE  (2),
        HIGH      (3);

        private final int severityRating;

        EffectImpact ( int severityRating ) {
            this.severityRating = severityRating;
        }

        public boolean isHigherImpactThan ( EffectImpact other ) {
            return this.severityRating > other.severityRating;
        }
    }

    // The kinds of variants supported by the SnpEff output format:
    public enum ChangeType {
        SNP,
        MNP,
        INS,
        DEL
    }

    // Possible zygosities of SnpEff variants:
    public enum Zygosity {
        Hom,
        Het
    }
}
