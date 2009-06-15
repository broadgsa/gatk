package org.broadinstitute.sting.playground.utils;

import java.util.HashMap;

public class IndelLikelihood {
    private String   type;
    private String[] alleles;
    private double   p;
    private double   lod;

    public IndelLikelihood(String type, String[] alleles, double p, double lod) {
        initialize(type, alleles, p, lod);
    }

    public IndelLikelihood(String[] indels) {
        HashMap<String,Integer> indel_allele_counts = new HashMap<String,Integer>();

        for (int i = 0; i < indels.length; i++) {
            if (! indel_allele_counts.containsKey(indels[i])) {
                indel_allele_counts.put(indels[i], 1);
            } else {
                indel_allele_counts.put(indels[i], indel_allele_counts.get(indels[i])+1);
            }
        }

        Object[] keys = indel_allele_counts.keySet().toArray();
        String[] alleles = new String[keys.length];
        int[] counts = new int[keys.length];
        //double likelihoods[] = new double[keys.length];
        int null_count = 0;
        String max_allele = null;
        int max_count = -1;

        if ((keys.length > 0) && (! ((keys.length == 1) && (((String)keys[0]).equals("null"))))) {
            for (int i = 0; i < keys.length; i++) {
                Integer count = (Integer)indel_allele_counts.get(keys[i]);
                alleles[i] = (String)keys[i];
                counts[i] = count;
                if (alleles[i].equals("null")) { null_count = counts[i]; }
                else if (counts[i] > max_count) { max_count = counts[i]; max_allele = alleles[i]; }
                //System.out.printf("%s[%d] ", keys[i], count);
            }
            //System.out.printf("\n");

            double eps = 1e-3;
            double pRef = null_count*Math.log10(1.0 - eps)   + max_count*Math.log10(eps) + Math.log10(0.999);
            double pHet = null_count*Math.log10(0.5 - eps/2) + max_count*Math.log10(0.5-eps/2) + Math.log10(1e-3);
            double pHom = null_count*Math.log10(eps)         + max_count*Math.log10(1.0 - eps) + Math.log10(1e-5);

            double lodRef = pRef - Math.max(pHet, pHom);
            double lodHet = pHet - pRef;
            double lodHom = pHom - pRef;

            //System.out.printf("%s/%s %f %f\n", "null", "null", pRef, lodRef);
            //System.out.printf("%s/%s %f %f\n", max_allele, "null", pHet, lodHet);
            //System.out.printf("%s/%s %f %f\n", max_allele, max_allele, pHom, lodHom);
            //System.out.printf("\n");

            if (lodRef > 0) {
                // reference call
                String[] genotype = new String[2];
                genotype[0] = "null";
                genotype[1] = "null";
                
                //return new IndelLikelihood("ref", genotype, pRef, lodRef);
                initialize("ref", genotype, pRef, lodRef);
            } else if (lodHet > lodHom) {
                // het call
                String[] genotype = new String[2];
                genotype[0] = "null";
                genotype[1] = max_allele;

                //return new IndelLikelihood("het", genotype, pHet, lodHet);
                initialize("het", genotype, pHet, lodHet);
            } else {
                // hom call
                String[] genotype = new String[2];
                genotype[0] = max_allele;
                genotype[1] = max_allele;

                //return new IndelLikelihood("hom", genotype, pHom, lodHom);
                initialize("hom", genotype, pHom, lodHom);
            }
        }
    }

    private void initialize(String type, String[] alleles, double p, double lod) {
        this.type = type;
        this.alleles = alleles;
        this.p = p;
        this.lod = lod;
    }

    public String getType() { return type; }
    public String[] getAlleles() { return alleles; }
    public double getPosteriorProbability() { return p; }
    public double getLOD() { return lod; }

    public String toString() {
        return String.format("%s/%s %f %f", alleles[0], alleles[1], p, lod);
    }
}
