package org.broadinstitute.sting.gatk.walkers.genotyper;
import org.broadinstitute.sting.utils.genotype.VariationCall;
import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jan 22, 2010
 * Time: 2:25:19 PM
 *
 * Useful helper class to communicate the results of calculateGenotype to framework
 */
public class VariantCallContext {
    public VariationCall variation = null;
    public List<Genotype> genotypes = null;

    /** Was the site called confidently, either reference or variant? */
    public boolean confidentlyCalled = false;

    VariantCallContext(VariationCall variation, List<Genotype> genotypes, boolean confidentlyCalledP) {
        this.variation = variation;
        this.genotypes = genotypes;
        this.confidentlyCalled = confidentlyCalledP;
    }

    /** blank variation and genotypes => we're a ref site */
    VariantCallContext(boolean confidentlyCalledP) {
        this.confidentlyCalled = confidentlyCalledP;
    }
}