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

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants;
import org.broadinstitute.sting.utils.codecs.snpEff.SnpEffFeature;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * A set of genomic annotations based on the output of the SnpEff variant effect predictor tool
 * (http://snpeff.sourceforge.net/).
 *
 * For each variant, chooses one of the effects of highest biological impact from the SnpEff
 * output file (which must be provided on the command line via --snpEffFile:SnpEff <filename>),
 * and adds annotations on that effect.
 *
 * The possible biological effects and their associated impacts are defined in the class:
 * org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants
 *
 * @author David Roazen
 */
public class SnpEff extends InfoFieldAnnotation implements ExperimentalAnnotation {

    // SnpEff annotation key names:
    public static final String GENE_ID_KEY = "GENE_ID";
    public static final String GENE_NAME_KEY = "GENE_NAME";
    public static final String TRANSCRIPT_ID_KEY = "TRANSCRIPT_ID";
    public static final String EXON_ID_KEY = "EXON_ID";
    public static final String EXON_RANK_KEY = "EXON_RANK";
    public static final String WITHIN_NON_CODING_GENE_KEY = "WITHIN_NON_CODING_GENE";
    public static final String EFFECT_KEY = "EFFECT";
    public static final String EFFECT_IMPACT_KEY = "EFFECT_IMPACT";
    public static final String EFFECT_EXTRA_INFORMATION_KEY = "EFFECT_EXTRA_INFORMATION";
    public static final String OLD_NEW_AA_KEY = "OLD_NEW_AA";
    public static final String OLD_NEW_CODON_KEY = "OLD_NEW_CODON";
    public static final String CODON_NUM_KEY = "CODON_NUM";
    public static final String CDS_SIZE_KEY = "CDS_SIZE";

    public Map<String, Object> annotate ( RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc ) {
        List<Feature> features = tracker.getValues(Feature.class);

        // Add only annotations for one of the most biologically-significant effects as defined in
        // the SnpEffConstants class:
        SnpEffFeature mostSignificantEffect = getMostSignificantEffect(features);

        if ( mostSignificantEffect == null ) {
            return null;
        }

        return generateAnnotations(mostSignificantEffect);
    }

    private SnpEffFeature getMostSignificantEffect ( List<Feature> features ) {
        SnpEffFeature mostSignificantEffect = null;

        for ( Feature feature : features ) {
            if ( feature instanceof SnpEffFeature ) {
                SnpEffFeature snpEffFeature = (SnpEffFeature)feature;

                if ( mostSignificantEffect == null ||
                     snpEffFeature.isHigherImpactThan(mostSignificantEffect) ) {

                    mostSignificantEffect = snpEffFeature;
                }
            }
        }

        return mostSignificantEffect;
    }

    private Map<String, Object> generateAnnotations ( SnpEffFeature mostSignificantEffect ) {
        Map<String, Object> annotations = new LinkedHashMap<String, Object>(Utils.optimumHashSize(getKeyNames().size()));

        if ( mostSignificantEffect.hasGeneID() )
            annotations.put(GENE_ID_KEY, mostSignificantEffect.getGeneID());
        if ( mostSignificantEffect.hasGeneName() )
            annotations.put(GENE_NAME_KEY, mostSignificantEffect.getGeneName());
        if ( mostSignificantEffect.hasTranscriptID() )
            annotations.put(TRANSCRIPT_ID_KEY, mostSignificantEffect.getTranscriptID());
        if ( mostSignificantEffect.hasExonID() )
            annotations.put(EXON_ID_KEY, mostSignificantEffect.getExonID());
        if ( mostSignificantEffect.hasExonRank() )
            annotations.put(EXON_RANK_KEY, Integer.toString(mostSignificantEffect.getExonRank()));
        if ( mostSignificantEffect.isNonCodingGene() )
            annotations.put(WITHIN_NON_CODING_GENE_KEY, null);

        annotations.put(EFFECT_KEY, mostSignificantEffect.getEffect().toString());
        annotations.put(EFFECT_IMPACT_KEY, mostSignificantEffect.getEffectImpact().toString());
        if ( mostSignificantEffect.hasEffectExtraInformation() )
            annotations.put(EFFECT_EXTRA_INFORMATION_KEY, mostSignificantEffect.getEffectExtraInformation());

        if ( mostSignificantEffect.hasOldAndNewAA() )
            annotations.put(OLD_NEW_AA_KEY, mostSignificantEffect.getOldAndNewAA());
        if ( mostSignificantEffect.hasOldAndNewCodon() )
            annotations.put(OLD_NEW_CODON_KEY, mostSignificantEffect.getOldAndNewCodon());
        if ( mostSignificantEffect.hasCodonNum() )
            annotations.put(CODON_NUM_KEY, Integer.toString(mostSignificantEffect.getCodonNum()));
        if ( mostSignificantEffect.hasCdsSize() )
            annotations.put(CDS_SIZE_KEY, Integer.toString(mostSignificantEffect.getCdsSize()));

        return annotations;
    }

    public List<String> getKeyNames() {
        return Arrays.asList( GENE_ID_KEY,
                              GENE_NAME_KEY,
                              TRANSCRIPT_ID_KEY,
                              EXON_ID_KEY,
                              EXON_RANK_KEY,
                              WITHIN_NON_CODING_GENE_KEY,
                              EFFECT_KEY,
                              EFFECT_IMPACT_KEY,
                              EFFECT_EXTRA_INFORMATION_KEY,
                              OLD_NEW_AA_KEY,
                              OLD_NEW_CODON_KEY,
                              CODON_NUM_KEY,
                              CDS_SIZE_KEY
                            );
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(
            new VCFInfoHeaderLine(GENE_ID_KEY,                  1, VCFHeaderLineType.String,  "Gene ID"),
            new VCFInfoHeaderLine(GENE_NAME_KEY,                1, VCFHeaderLineType.String,  "Gene name"),
            new VCFInfoHeaderLine(TRANSCRIPT_ID_KEY,            1, VCFHeaderLineType.String,  "Transcript ID"),
            new VCFInfoHeaderLine(EXON_ID_KEY,                  1, VCFHeaderLineType.String,  "Exon ID"),
            new VCFInfoHeaderLine(EXON_RANK_KEY,                1, VCFHeaderLineType.Integer, "Exon rank"),
            new VCFInfoHeaderLine(WITHIN_NON_CODING_GENE_KEY,   0, VCFHeaderLineType.Flag,    "If present, gene is non-coding"),
            new VCFInfoHeaderLine(EFFECT_KEY,                   1, VCFHeaderLineType.String,  "One of the most high-impact effects across all transcripts at this site"),
            new VCFInfoHeaderLine(EFFECT_IMPACT_KEY,            1, VCFHeaderLineType.String,  "Impact of the effect " + Arrays.toString(SnpEffConstants.EffectImpact.values())),
            new VCFInfoHeaderLine(EFFECT_EXTRA_INFORMATION_KEY, 1, VCFHeaderLineType.String,  "Additional information about the effect"),
            new VCFInfoHeaderLine(OLD_NEW_AA_KEY,               1, VCFHeaderLineType.String,  "Old/New amino acid"),
            new VCFInfoHeaderLine(OLD_NEW_CODON_KEY,            1, VCFHeaderLineType.String,  "Old/New codon"),
            new VCFInfoHeaderLine(CODON_NUM_KEY,                1, VCFHeaderLineType.Integer, "Codon number"),
            new VCFInfoHeaderLine(CDS_SIZE_KEY,                 1, VCFHeaderLineType.Integer, "CDS size")
        );
    }
}
