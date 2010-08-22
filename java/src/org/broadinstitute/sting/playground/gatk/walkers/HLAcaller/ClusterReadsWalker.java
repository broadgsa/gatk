package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.ArrayList;
import java.util.Hashtable;
import java.io.PrintStream;

/**
 * Compares reads to longest read at each exon. Usage: java -jar GenomeAnalysisTK.jar -T ClusterReads -I INPUT.bam -R /broad/1KG/reference/human_b36_both.fasta [-filter INPUT.filter] | grep -v INFO | sort -k1 > OUTPUT
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class ClusterReadsWalker extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName = "filter", shortName = "filter", doc = "file containing reads to exclude", required = false)
    public String filterFile = "";

    @Argument(fullName = "maxAllowedMismatches", shortName = "maxAllowedMismatches", doc = "Max number of mismatches tolerated per read (default 7)", required = false)
    public int MAXALLOWEDMISMATCHES = 7;

    @Argument(fullName = "minRequiredMatches", shortName = "minRequiredMatches", doc = "Min number of matches required per read (default 7)", required = false)
    public int MINREQUIREDMATCHES = 0;

    String UniqueAllelesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/UniqueAlleles";

    String PolymorphicSitesFile = "/humgen/gsa-scr1/GSA/sjia/Sting/HLA.polymorphic.sites";

    boolean DatabaseLoaded = false;
    boolean DEBUG = false;

    Integer[] HLAstartpos, HLAstoppos, PolymorphicSites,NonPolymorphicSites;
    double[] SingleAlleleFrequencies;
    ArrayList<String> ReadsToDiscard = new ArrayList<String>();
    ArrayList<SAMRecord> AlignedReads = new ArrayList<SAMRecord>();

    Hashtable MaxNumMatches, MaxConcordance, NumCompared = new Hashtable();

    double[] nummatched, concordance, numcompared;

    SAMRecord A2, A3, B2, B3, C2, C3;
    int MaxMatchesA2 = 0, MaxMatchesA3 = 0, MaxMatchesB2 = 0, MaxMatchesB3 = 0, MaxMatchesC2 = 0, MaxMatchesC3 = 0;
    int A2start = 30018513, A2stop = 30018781;
    int A3start = 30019024, A3stop = 30019300;

    int C2start = 31347355, C2stop = 31347623;
    int C3start = 31346829, C3stop = 31347104;

    int B2start = 31432444, B2stop = 31432714;
    int B3start = 31431923, B3stop = 31432198;
    
    int minstartpos = 0;
    int maxstoppos = 0;


    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;
    CigarParser formatter = new CigarParser();

    public Integer reduceInit() {
        if (!DatabaseLoaded){
            DatabaseLoaded = true;

            PolymorphicSitesFileReader siteFileReader = new PolymorphicSitesFileReader();
            siteFileReader.ReadFile(PolymorphicSitesFile);
            PolymorphicSites = siteFileReader.GetPolymorphicSites();
            NonPolymorphicSites = siteFileReader.GetNonPolymorphicSites();


            out.printf("INFO  %s polymorphic and %s non-polymorphic sites found in HLA dictionary\n",PolymorphicSites.length,NonPolymorphicSites.length);

            if (!filterFile.equals("")){
                out.printf("INFO  Reading properties file ... ");
                SimilarityFileReader similarityReader = new SimilarityFileReader();
                similarityReader.ReadFile(filterFile,MAXALLOWEDMISMATCHES,MINREQUIREDMATCHES);
                ReadsToDiscard = similarityReader.GetReadsToDiscard();
                MaxNumMatches = similarityReader.GetNumMatches();
                MaxConcordance = similarityReader.GetConcordance();

                out.printf("Done! Found %s misaligned reads to discard.\n",ReadsToDiscard.size());
                for (int i = 0; i < ReadsToDiscard.size(); i++){
                    out.printf("MISALIGNED %s\n", ReadsToDiscard.get(i).toString());
                }
            }

            out.printf("INFO  Comparing reads ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {
        //Calculate concordance for this read and all overlapping reads
        if (!ReadsToDiscard.contains(read.getReadName())){
            AlignedReads.add(read);
            
            int readstart = read.getAlignmentStart();
            int readstop = read.getAlignmentEnd();
            int length = readstop - readstart + 1;
            if (MaxNumMatches.containsKey(read.getReadName())){
                int maxMatches = Integer.parseInt(MaxNumMatches.get(read.getReadName()).toString());
                double concordance = Double.parseDouble((String) MaxConcordance.get(read.getReadName()).toString());

                if (readstart < A2stop && readstop > A2start){
                    if (maxMatches > MaxMatchesA2 && concordance > 0.95){
                        MaxMatchesA2 = maxMatches;
                        A2 = read;
                    }
                } else if (readstart < A3stop && readstop > A3start){
                    if (maxMatches > MaxMatchesA3){
                        MaxMatchesA3 = maxMatches;
                        A3 = read;
                    }
                } else if (readstart < B2stop && readstop > B2start){
                    if (maxMatches > MaxMatchesB2){
                        MaxMatchesB2 = maxMatches;
                        B2 = read;
                    }
                } else if (readstart < B3stop && readstop > B3start){
                    if (maxMatches > MaxMatchesB3){
                        MaxMatchesB3 = maxMatches;
                        B3 = read;
                    }
                } else if (readstart < C2stop && readstop > C2start){
                    if (maxMatches > MaxMatchesC2){
                        MaxMatchesC2 = maxMatches;
                        C2 = read;
                    }
                } else if (readstart < C3stop && readstop > C3start){
                    if (maxMatches > MaxMatchesC3){
                        MaxMatchesC3 = maxMatches;
                        C3 = read;
                    }
                }
            }else{
                out.printf("Data for %s not found\n",read.getReadName());
            }
        }
        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        
        return value + sum;
    }

    public void onTraversalDone(Integer numreads) {
        SAMRecord read; String name, name2; String locus = "";
        int A2a = 0, A2b = 0, A2c = 0;
        int A3a = 0, A3b = 0, A3c = 0;
        int B2a = 0, B2b = 0, B2c = 0;
        int B3a = 0, B3b = 0, B3c = 0;
        int C2a = 0, C2b = 0, C2c = 0;
        int C3a = 0, C3b = 0, C3c = 0;
        double minA2 = -1, minA3 = -1, minB2 = -1, minB3 = -1, minC2 = -1, minC3 = -1;
        double maxA2 = 0, maxA3 = 0, maxB2 = 0, maxB3 = 0, maxC2 = 0, maxC3 = 0;
        double ratioA2 = 0, ratioA3 = 0, ratioB2 = 0, ratioB3 = 0, ratioC2 = 0, ratioC3 = 0;
        double maxA2l = 0, maxA3l = 0, maxB2l = 0, maxB3l = 0, maxC2l = 0, maxC3l = 0;
        double maxA2d = 0, maxA3d = 0, maxB2d = 0, maxB3d = 0, maxC2d = 0, maxC3d = 0;
                double a2, a3, b2, b3, c2, c3, normalized = 0;
        int readstart, readstop;
        double matches, compared, concordance;
        STATS stats;

        for (int i = 0; i < AlignedReads.size(); i++){
            read = AlignedReads.get(i);
            readstart = read.getAlignmentStart();
            readstop = read.getAlignmentEnd();
            if (readstart < A2stop && readstop > A2start){
                stats = CalculateConcordance(read, A2, "A2", A2start, A2stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minA2 < 0 || minA2 > concordance){minA2 = concordance;}; if (concordance > maxA2){maxA2 = concordance;}; if (matches > maxA2l){maxA2l = matches;}; if (compared-matches > maxA2d){maxA2d = compared-matches;}}else{A2c++;}
            } else if (readstart < A3stop && readstop > A3start){
                stats = CalculateConcordance(read, A3, "A3", A3start, A3stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minA3 < 0 || minA3 > concordance){minA3 = concordance;}; if (concordance > maxA3){maxA3 = concordance;}; if (matches > maxA3l){maxA3l = matches;}; if (compared-matches > maxA3d){maxA3d = compared-matches;}}else{A3c++;}
            } else if (readstart < B2stop && readstop > B2start){
                stats = CalculateConcordance(read, B2, "B2", B2start, B2stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minB2 < 0 || minB2 > concordance){minB2 = concordance;}; if (concordance > maxB2){maxB2 = concordance;}; if (matches > maxB2l){maxB2l = matches;}; if (compared-matches > maxB2d){maxB2d = compared-matches;}}else{B2c++;}
            } else if (readstart < B3stop && readstop > B3start){
                stats = CalculateConcordance(read, B3, "B3", B3start, B3stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minB3 < 0 || minB3 > concordance){minB3 = concordance;}; if (concordance > maxB3){maxB3 = concordance;}; if (matches > maxB3l){maxB3l = matches;}; if (compared-matches > maxB3d){maxB3d = compared-matches;}}else{B3c++;}
            } else if (readstart < C2stop && readstop > C2start){
                stats = CalculateConcordance(read, C2, "C2", C2start, C2stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minC2 < 0 || minC2 > concordance){minC2 = concordance;}; if (concordance > maxC2){maxC2 = concordance;}; if (matches > maxC2l){maxC2l = matches;}; if (compared-matches > maxC2d){maxC2d = compared-matches;}}else{C2c++;}
            } else if (readstart < C3stop && readstop > C3start){
                stats = CalculateConcordance(read, C3, "C3", C3start, C3stop); concordance = stats.getConcordance(); matches = stats.getNumMatches(); compared = stats.getNumCompared();
                if (stats.getNumCompared() > 40){if (minC3 < 0 || minC3 > concordance){minC3 = concordance;}; if (concordance > maxC3){maxC3 = concordance;}; if (matches > maxC3l){maxC3l = matches;}; if (compared-matches > maxC3d){maxC3d = compared-matches;}}else{C3c++;}
            }
        }


        for (int i = 0; i < AlignedReads.size(); i++){
            read = AlignedReads.get(i);
            readstart = read.getAlignmentStart();
            readstop = read.getAlignmentEnd();
            name = read.getReadName(); name2 = "";
            if (NumCompared.containsKey(name)){
                compared = Double.parseDouble((String) NumCompared.get(name).toString());
                matches = Double.parseDouble((String) MaxNumMatches.get(name).toString());
                concordance = Double.parseDouble((String) MaxConcordance.get(name).toString());
                if (matches > 40){
                    if (readstart < A2stop && readstop > A2start){
                        locus = "A2"; name2 = A2.getReadName();
                        a2 = (concordance - minA2)/(maxA2-minA2); if (a2 >= .5){A2a++;}else{A2b++;}; normalized = a2;
                    } else if (readstart < A3stop && readstop > A3start){
                        locus = "A3"; name2 = A3.getReadName();
                        a3 = (concordance - minA3)/(maxA3-minA3); if (a3 >= .5){A3a++;}else{A3b++;}; normalized = a3;
                    } else if (readstart < B2stop && readstop > B2start){
                        locus = "B2"; name2 = B2.getReadName();
                        b2 = (concordance - minB2)/(maxB2-minB2); if (b2 >= .5){B2a++;}else{B2b++;}; normalized = b2;
                    } else if (readstart < B3stop && readstop > B3start){
                        locus = "B3"; name2 = B3.getReadName();
                        b3 = (concordance - minB3)/(maxB3-minB3); if (b3 >= .5){B3a++;}else{B3b++;}; normalized = b3;
                    } else if (readstart < C2stop && readstop > C2start){
                        locus = "C2"; name2 = C2.getReadName();
                        c2 = (concordance - minC2)/(maxC2-minC2); if (c2 >= .5){C2a++;}else{C2b++;}; normalized = c2;
                    } else if (readstart < C3stop && readstop > C3start){
                        locus = "C3"; name2 = C3.getReadName();
                        c3 = (concordance - minC3)/(maxC3-minC3); if (c3 >= .5){C3a++;}else{C3b++;}; normalized = c3;
                    }
                    out.printf("%s\t%s\t%s\t%.0f\t%.0f\t%.3f\t%.3f\n",locus,name,name2,matches,compared,concordance,normalized);
                }else{
                    out.printf("%s (compared at %s sites) is too short\n",name,matches);
                }
            }else{
                out.printf("%s [%s to %s] not found\n",name,readstart,readstop);
            }
        }

        if (A2a > 0 && A2b > 0){if (A2a > A2b){ratioA2 = (double)A2b/(A2a+A2b);}else{ratioA2 = (double)A2a/(A2a+A2b);}}else{ratioA2 = -1;}
        if (A3a > 0 && A3b > 0){if (A3a > A3b){ratioA3 = (double)A3b/(A3a+A3b);}else{ratioA3 = (double)A3a/(A3a+A3b);}}else{ratioA3 = -1;}
        if (B2a > 0 && B2b > 0){if (B2a > B2b){ratioB2 = (double)B2b/(B2a+B2b);}else{ratioB2 = (double)B2a/(B2a+B2b);}}else{ratioB2 = -1;}
        if (B3a > 0 && B3b > 0){if (B3a > B3b){ratioB3 = (double)B3b/(B3a+B3b);}else{ratioB3 = (double)B3a/(B3a+B3b);}}else{ratioB3 = -1;}
        if (C2a > 0 && C2b > 0){if (C2a > C2b){ratioC2 = (double)C2b/(C2a+C2b);}else{ratioC2 = (double)C2a/(C2a+C2b);}}else{ratioC2 = -1;}
        if (C3a > 0 && C3b > 0){if (C3a > C3b){ratioC3 = (double)C3b/(C3a+C3b);}else{ratioC3 = (double)C3a/(C3a+C3b);}}else{ratioC3 = -1;}

        out.printf("RATIO_A2\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioA2,A2a,A2b,A2c,maxA2-minA2,maxA2l,maxA2d);
        out.printf("RATIO_A3\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioA3,A3a,A3b,A3c,maxA3-minA3,maxA3l,maxA3d);
        out.printf("RATIO_B2\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioB2,B2a,B2b,B2c,maxB2-minB2,maxB2l,maxB2d);
        out.printf("RATIO_B3\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioB3,B3a,B3b,B3c,maxB3-minB3,maxB3l,maxB3d);
        out.printf("RATIO_C2\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioC2,C2a,C2b,C2c,maxC2-minC2,maxC2l,maxC2d);
        out.printf("RATIO_C3\t%.2f\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\n",ratioC3,C3a,C3b,C3c,maxC3-minC3,maxC3l,maxC3d);

    }

    public class STATS {
       protected double concordance = 0.0;
       protected double numcompared = 0;
       protected double nummatches = 0;

       public STATS(double d, double i, double m) {
          concordance = d;
          numcompared = i;
          nummatches = m;
       }

       public double getConcordance() {
          return concordance;
       }

       public double getNumCompared() {
          return numcompared;
       }

       public double getNumMatches() {
          return nummatches;
       }
    }

    private STATS CalculateConcordance(SAMRecord read1, SAMRecord read2, String locus, int start, int stop){
        int start1 = read1.getAlignmentStart(), stop1 = read1.getAlignmentEnd();
        int start2 = read2.getAlignmentStart(), stop2 = read2.getAlignmentEnd();

        int pos;
        double numcompared = 0, nummatched = 0, concordance;
        char c1, c2;
        String s1 = formatter.FormatRead(read1.getCigarString(), read1.getReadString());
        String s2 = formatter.FormatRead(read2.getCigarString(), read2.getReadString());


        //Polymorphic sites: always increment denominator, increment numerator when bases are concordant
        for (int j = 0; j < PolymorphicSites.length; j++){
            pos = PolymorphicSites[j];
            if (pos >= start1 && pos <= stop1 && pos >= start2 && pos <= stop2 && pos >= start && pos <= stop){
                c1 = s1.charAt(pos-start1);
                c2 = s2.charAt(pos-start2);
                if (c1 != 'D'){//allow for deletions (sequencing errors)
                    numcompared++;
                    if (c1 == c2){
                        nummatched++;
                    }
                }
            }
        }

        //Non-polymorphic sites: increment denominator only when bases are discordant
        if (false){
            for (int j = 0; j < NonPolymorphicSites.length; j++){
                pos = NonPolymorphicSites[j];
                if (pos >= start1 && pos <= stop1 && pos >= start2 && pos <= stop2){
                    c1 = s1.charAt(pos-start1);
                    c2 = s2.charAt(pos-start2);
                    if (c1 != c2 && c1 != 'D'){//allow for deletions (sequencing errors)
                        numcompared++;
                    }
                }
            }
        }
            
        //Update concordance array
        concordance=nummatched/numcompared;

        MaxNumMatches.put(read1.getReadName(), nummatched);
        NumCompared.put(read1.getReadName(), numcompared);
        MaxConcordance.put(read1.getReadName(), concordance);
        //out.printf("%s\t%s\t%s\t%.0f\t%.0f\t%.3f\n",locus,read1.getReadName(),read2.getReadName(),nummatched,numcompared,concordance);

        return new STATS(concordance, numcompared, nummatched);
    }

}


