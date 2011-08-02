package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 6/1/11
 * Time: 10:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class MultiallelicGenotypeLikelihoods {
    private String sample;
    private double[] GLs;
    private ArrayList<Allele> alleleList;
    private int depth;

    public MultiallelicGenotypeLikelihoods(String sample,
                                         ArrayList<Allele> A,
                                         double[] log10Likelihoods, int depth) {
        /* Check for consistency between likelihood vector and number of alleles */
        int numAlleles = A.size();
        if (log10Likelihoods.length != numAlleles*(numAlleles+1)/2)
            throw new StingException(("BUG: Incorrect length of GL vector when creating MultiallelicGenotypeLikelihoods object!"));

         this.sample = sample;
         this.alleleList = A;
         this.GLs = log10Likelihoods;
         this.depth = depth;
     }

     public String getSample() {
         return sample;
     }

      public double[] getLikelihoods() {
         return GLs;
     }

     public ArrayList<Allele> getAlleles() {
         return alleleList;
     }

     public int getDepth() {
         return depth;
     }

}
