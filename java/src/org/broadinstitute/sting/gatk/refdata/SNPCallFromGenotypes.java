package org.broadinstitute.sting.gatk.refdata;

import java.util.List;

/**
 * loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom
 * chr1:1104840 A N 0.000000 -85.341265 -85.341265 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104841 C N 0.000000 -69.937928 -69.937928 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104842 A N 0.000000 -84.816002 -84.816002 0.000000 0.000000 324.000000 162 0 0
 *
 */
public interface SNPCallFromGenotypes extends AllelicVariant {
    public int nIndividuals();
    public int nHomRefGenotypes();
    public int nHetGenotypes();
    public int nHomVarGenotypes();
    public List<Genotype> getGenotypes();
}