package edu.mit.broad.sting.atk;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Mar 9, 2009
 * Time: 3:34:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenotypeEvidence {

    int[] nuc2num = new int[128];
    int[] nucs = new int[4];
    int a = nucs[0];
    int c = nucs[1];
    int t = nucs[2];
    int g = nucs[3];
    float[] nuc_pcnt = new float[4];
    char ref;
    public float q; // % non-reference alleles
    public int refbases;
    public int allbases;

    public GenotypeEvidence(String bases, char ref){
        this.ref = ref;
        nuc2num['A'] = 0;
        nuc2num['C'] = 1;
        nuc2num['T'] = 2;
        nuc2num['G'] = 3;
        nuc2num['a'] = 0;
        nuc2num['c'] = 1;
        nuc2num['t'] = 2;
        nuc2num['g'] = 3;

        for (char b : bases.toCharArray()) {
            nucs[nuc2num[b]] += 1;
            /*switch (b) {
                case 'A': nucs[0] += 1; break;
                case 'C': nucs[1] += 1; break;
                case 'T': nucs[2] += 1; break;
                case 'G': nucs[3] += 1; break;
            } */
        }

        // Calculate q = ref. bases / nonref. bases
        refbases = nucs[nuc2num[ref]];
        allbases = bases.length();
        q = 1 - ((float)refbases / allbases);

        /*for (int i=0; i<4; i++) {
            nuc_pcnt[i] = (float)nucs[i] / len;
            //if
        }*/
    }


    public boolean SigNonref(float cutoff_fraction) {
      /*  for (char nuc : nucs) {
            
        }*/


        return true;
    }

    public void print() {

        System.out.format("A %2d | ", nucs[0]);
        System.out.format("C %2d | ", nucs[1]);
        System.out.format("T %2d | ", nucs[2]);
        System.out.format("G %2d | ", nucs[3]);
        System.out.format("Ref %s | ", ref);

    }



}
