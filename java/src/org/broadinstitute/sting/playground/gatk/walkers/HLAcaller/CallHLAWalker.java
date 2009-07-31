/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.genotype.calls.SSGGenotypeCall;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.gatk.walkers.*;
import java.io.*;
import java.util.*;
/**
 *
 * @author shermanjia
 */
public class CallHLAWalker extends LocusWalker<Integer, Pair<Long, Long>>{
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
    public boolean suppressPrinting = false;

    //String HLAdatabaseFile = "/Users/shermanjia/Work/HLA.sam";
    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.sam";


    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    int[] HLAstartpos;
    int[] HLAstoppos;
    int numHLAlleles = 0;
    double[][] Prob; String[][] Alleles;
    boolean THREE_BASE_ERRORS = false;
    String PRIORS_2ND_ON  = "0.000,0.302,0.366,0.142,0.000,0.548,0.370,0.000,0.319,0.000";
    String PRIORS_2ND_OFF = "0.480,0.769,0.744,0.538,0.575,0.727,0.768,0.589,0.762,0.505";
    double[] p2ndon = priorsArray(PRIORS_2ND_ON);
    double[] p2ndoff = priorsArray(PRIORS_2ND_OFF);
    boolean keepQ0Bases = false;

    //HLAreads.add("Italian Riviera");

    public Pair<Long, Long> reduceInit() {
        try{
            FileInputStream fstream = new FileInputStream(HLAdatabaseFile);
            // Get the object of DataInputStream
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine; String [] s = null;
            //Read File Line By Line
            while ((strLine = br.readLine()) != null)   {
                s = strLine.split("\\t");
                if (s.length>=10){
                    HLAreads.add(s[9]);
                    HLAnames.add(s[0]);
                    HLApositions.add(s[3]);
                    //System.out.println (s[0] + "\t" + s[3]);
                }
            }
            in.close();
            int n = HLApositions.size();
            numHLAlleles = n;
            HLAstartpos = new int[n]; HLAstoppos = new int[n];
            Prob = new double[n][n];
            for (int i = 0; i < n; i++){
                HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                for (int j = 0; j <n; j++){
                    Prob[i][j]=1;
                }
                //System.out.println (HLAstartpos[i] + " to " + HLAstoppos[i]);
            }

        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }
        return new Pair<Long,Long>(0l,0l);
    }


    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        long loc = context.getPosition();
        if ( !suppressPrinting )
            out.printf("%s: %d\t", context.getLocation(),pileup.size() );
        out.printf("ref=%s\t", ref);
        if( context.getReads().size() > 0 ) {
            //out.printf("RG for first read: %s%n",context.getReads().get(0).getReadName());
            int numAs = 0;
            int numCs = 0;
            int numGs = 0;
            int numTs = 0;
            char c = 'a';
            
            for (int i = 0;i < context.getReads().size();i++){
                //out.printf("%s%s ",i,bases.charAt(i));
                char base = bases.charAt(i);
                if (base == 'A'){numAs++;}
                if (base == 'C'){numCs++;}
                if (base == 'T'){numTs++;}
                if (base == 'G'){numGs++;}

            }
            out.printf("%sAs\t%sCs\t%sTs\t%sGs\t",numAs,numCs,numTs,numGs);

            GenotypeLikelihoods G = new GenotypeLikelihoods(THREE_BASE_ERRORS,0.999,0.000333,0.000667, p2ndon, p2ndoff, keepQ0Bases);
            SSGGenotypeCall geno = (SSGGenotypeCall)G.callGenotypes(tracker, ref, pileup);

            
            double mLikelihoods[] = geno.getLikelihoods();
            List<Genotype> mGenotypes = geno.getGenotypes();
            for (int j =0; j < geno.getBases().length();j++){
                out.printf("%s %s\t",mGenotypes.get(j).getBases(),mLikelihoods[j]);
            }
                    

            for (int j = 0; j < numHLAlleles; j++){
                //out.print(loc + "," + HLAstartpos[j] + "," + HLAstoppos[j]);
                long i = loc - HLAstartpos[j];
                if (loc >= HLAstartpos[j] && loc <= HLAstoppos[j]){
                    c = HLAreads.get(j).charAt((int) i);
                    //out.printf("%s",c);
                }
            }
            out.print("\n");
            //out.printf("RG for first read: %s%n",context.getReads().get(0).getReadName());

        }
        return context.getReads().size();

    }

    private double[] priorsArray(String priorsString) {
        String[] pstrs = priorsString.split(",");
        double[] pdbls = new double[pstrs.length];

        for (int i = 0; i < pstrs.length; i++) {
            pdbls[i] = Double.valueOf(pstrs[i]);
        }

        return pdbls;
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",
                ((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }
}
