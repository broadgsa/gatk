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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff;
import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff.EffectType;
import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff.InfoFieldKey;
import org.broadinstitute.gatk.tools.walkers.annotator.SnpEffUtil;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies variants as genes or coding regions, according to the effect modifier, as indicated by snpEff.
 * The 'gene' category includes category 'coding region', and additionally includes introns. 'Coding regions'
 * includes transcripts and, implicitly, UTRs.
 */
public class SnpEffPositionModifier extends VariantStratifier {

	public enum PositionModifier {
		GENE,          // EXON
		CODING_REGION, // CDS
		SPLICE_SITE,   // not a straight translation -- see getRelevantStates
		STOP_GAINED,   // STOP_GAINED
		STOP_LOST      // STOP_LOST
	}

	@Override
	public void initialize() {
		for (final PositionModifier type : PositionModifier.values()) states.add(type.name());
	}

	@Override
	public List<Object> getRelevantStates(
			final ReferenceContext ref,
			final RefMetaDataTracker tracker,
			final VariantContext comp,
			final String compName,
			final VariantContext eval,
			final String evalName,
			final String sampleName,
            final String FamilyName)
	{
		final List<Object> relevantStates = new ArrayList<Object>();
		if (eval != null && eval.isVariant() && eval.hasAttribute(InfoFieldKey.EFFECT_KEY.getKeyName())) {
			final SnpEff.EffectType effectType = SnpEff.EffectType.valueOf(
					eval.getAttribute(InfoFieldKey.EFFECT_KEY.getKeyName()).toString());

			if (SnpEffUtil.isSubTypeOf(effectType, EffectType.EXON))        relevantStates.add(PositionModifier.GENE.name());
			if (SnpEffUtil.isSubTypeOf(effectType, EffectType.CDS))         relevantStates.add(PositionModifier.CODING_REGION.name());
			if (SnpEffUtil.isSubTypeOf(effectType, EffectType.STOP_GAINED)) relevantStates.add(PositionModifier.STOP_GAINED.name());
			if (SnpEffUtil.isSubTypeOf(effectType, EffectType.STOP_LOST))   relevantStates.add(PositionModifier.STOP_LOST.name());

			if (SnpEffUtil.isSubTypeOf(effectType, EffectType.SPLICE_SITE_ACCEPTOR) ||
				SnpEffUtil.isSubTypeOf(effectType, EffectType.SPLICE_SITE_DONOR))
					relevantStates.add(PositionModifier.SPLICE_SITE.name());
		}

		return relevantStates;
	}
}
