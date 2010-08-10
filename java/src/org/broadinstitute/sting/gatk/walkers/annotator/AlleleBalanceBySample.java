package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFFormatHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broad.tribble.util.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.*;

import java.util.*;


public class AlleleBalanceBySample implements GenotypeAnnotation, ExperimentalAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
				Double ratio = annotateSNP(stratifiedContext, vc, g);
				if (ratio == null)
						return null;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", ratio.doubleValue()));
				return map;

   }

    private Double annotateSNP(StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g) {

				double ratio = -1;

        if ( !vc.isSNP() ) 
						return null;

				if ( !vc.isBiallelic() )
						return null;

        if ( g == null || !g.isCalled() )
            return null;

				if (!g.isHet())
						return null;

        Set<Allele> altAlleles = vc.getAlternateAlleles();
        if ( altAlleles.size() == 0 )
            return null;

				final String bases = new String(stratifiedContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().getBases());
				if ( bases.length() == 0 )
						return null;
				char refChr = vc.getReference().toString().charAt(0);
				char altChr = vc.getAlternateAllele(0).toString().charAt(0);

				int refCount = MathUtils.countOccurrences(refChr, bases);
				int altCount = MathUtils.countOccurrences(altChr, bases);

				// sanity check
				if ( refCount + altCount == 0 )
						return null;

				ratio = ((double)refCount / (double)(refCount + altCount));
				return ratio;
    }

    public List<String> getKeyNames() { return Arrays.asList("AB"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), -1, VCFHeaderLineType.Float, "Allele balance for each het genotype")); }
}