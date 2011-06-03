package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broad.tribble.util.variantcontext.Allele;

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
                                         double[] log10AALikelihoods, int depth) {
         this.sample = sample;
         this.alleleList = A;
         this.GLs = log10AALikelihoods;
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
