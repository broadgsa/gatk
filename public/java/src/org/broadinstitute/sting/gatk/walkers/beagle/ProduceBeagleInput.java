/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Gender;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VQSRCalibrationCurve;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 *  Converts the input VCF into a format accepted by the Beagle imputation/analysis program.
 * <p>
 *
 * <h2>Input</h2>
 * <p>
 * A VCF with variants to convert to Beagle format
 * </p>
 *
 * <h2>Outputs</h2>
 * <p>
 * A single text file which can be fed to Beagle
 * </p>
 * <p>
 * Optional: A file with a list of markers
 * </p>
  *
 * <h2>Examples</h2>
 * <pre>
 *     java -Xmx2g -jar dist/GenomeAnalysisTK.jar -L 20 \
 *      -R reffile.fasta -T ProduceBeagleInput \
 *      -V path_to_input_vcf/inputvcf.vcf -o path_to_beagle_output/beagle_output
 * </pre>
 *
 */

@DocumentedGATKFeature( groupName = "Variant Discovery Tools", extraDocs = {CommandLineGATK.class} )
public class ProduceBeagleInput extends RodWalker<Integer, Integer> {

    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Hidden
    @Input(fullName="validation", shortName = "validation", doc="Validation VCF file", required=false)
    public RodBinding<VariantContext> validation;


    @Output(doc="File to which BEAGLE input should be written",required=true)
    protected PrintStream  beagleWriter = null;

    @Hidden
     @Output(doc="File to which BEAGLE markers should be written", shortName="markers", fullName = "markers", required = false)
    protected PrintStream  markers = null;
    int markerCounter = 1;

    @Hidden
    @Input(doc="VQSqual calibration file", shortName = "cc", required=false)
    protected File VQSRCalibrationFile = null;
    protected VQSRCalibrationCurve VQSRCalibrator = null;

    @Hidden
    @Argument(doc="VQSqual key", shortName = "vqskey", required=false)
    protected String VQSLOD_KEY = "VQSqual";

    @Hidden
     @Argument(fullName = "inserted_nocall_rate", shortName = "nc_rate", doc = "Rate (0-1) at which genotype no-calls will be randomly inserted, for testing", required = false)
    public double insertedNoCallRate  = 0;
    @Hidden
     @Argument(fullName = "validation_genotype_ptrue", shortName = "valp", doc = "Flat probability to assign to validation genotypes. Will override GL field.", required = false)
    public double validationPrior = -1.0;
    @Hidden
     @Argument(fullName = "validation_bootstrap", shortName = "bs", doc = "Proportion of records to be used in bootstrap set", required = false)
    public double bootstrap = 0.0;
    @Hidden
     @Argument(fullName = "bootstrap_vcf",shortName = "bvcf", doc = "Output a VCF with the records used for bootstrapping filtered out", required = false)
    VariantContextWriter bootstrapVCFOutput = null;

    /**
     * If sample gender is known, this flag should be set to true to ensure that Beagle treats male Chr X properly.
     */
    @Argument(fullName = "checkIsMaleOnChrX", shortName = "checkIsMaleOnChrX", doc = "Set to true when Beagle-ing chrX and want to ensure male samples don't have heterozygous calls.", required = false)
    public boolean CHECK_IS_MALE_ON_CHR_X = false;

    @Hidden
    @Argument(fullName = "variant_genotype_ptrue", shortName = "varp", doc = "Flat probability prior to assign to variant (not validation) genotypes. Does not override GL field.", required = false)
    public double variantPrior = 0.96;

    private Set<String> samples = null;
    private Set<String> BOOTSTRAP_FILTER = new HashSet<String>( Arrays.asList("bootstrap") );
    private int bootstrapSetSize = 0;
    private int testSetSize = 0;
    private CachingFormatter formatter = new CachingFormatter("%5.4f ", 100000);
    private int certainFPs = 0;

    public void initialize() {

        samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(variantCollection.variants.getName()));

        beagleWriter.print("marker alleleA alleleB");
        for ( String sample : samples )
            beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));

        beagleWriter.println();

        if ( bootstrapVCFOutput != null ) {
            initializeVcfWriter();
        }

        if ( VQSRCalibrationFile != null ) {
            VQSRCalibrator = VQSRCalibrationCurve.readFromFile(VQSRCalibrationFile);
            logger.info("Read calibration curve");
            VQSRCalibrator.printInfo(logger);
        }
    }

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        if( tracker != null ) {
            GenomeLoc loc = context.getLocation();
            VariantContext variant_eval = tracker.getFirstValue(variantCollection.variants, loc);
            VariantContext validation_eval = tracker.getFirstValue(validation, loc);

            if ( goodSite(variant_eval,validation_eval) ) {
                if ( useValidation(validation_eval, ref) ) {
                    writeBeagleOutput(validation_eval, variant_eval, true, validationPrior);
                    return 1;
                } else {
                    if ( goodSite(variant_eval) ) {
                        writeBeagleOutput(variant_eval,validation_eval,false,variantPrior);
                        return 1;
                    } else { // todo -- if the variant site is bad, validation is good, but not in bootstrap set -- what do?
                        return 0;
                    }
                }
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    public boolean goodSite(VariantContext a, VariantContext b) {
        return goodSite(a) || goodSite(b);
    }

    public boolean goodSite(VariantContext v) {
        if ( canBeOutputToBeagle(v) ) {
            if ( VQSRCalibrator != null && VQSRCalibrator.certainFalsePositive(VQSLOD_KEY, v) ) {
                certainFPs++;
                return false;
            } else {
                return true;
            }
        } else {
            return false;
        }
    }

    public static boolean canBeOutputToBeagle(VariantContext v) {
        return v != null && ! v.isFiltered() && v.isBiallelic() && v.hasGenotypes();
    }

    public boolean useValidation(VariantContext validation, ReferenceContext ref) {
        if( goodSite(validation) ) {
            // if using record keeps us below expected proportion, use it
            logger.debug(String.format("boot: %d, test: %d, total: %d", bootstrapSetSize, testSetSize, bootstrapSetSize+testSetSize+1));
            if ( (bootstrapSetSize+1.0)/(1.0+bootstrapSetSize+testSetSize) <= bootstrap ) {
                if ( bootstrapVCFOutput != null ) {
                    bootstrapVCFOutput.add(new VariantContextBuilder(validation).filters(BOOTSTRAP_FILTER).make());
                }
                bootstrapSetSize++;
                return true;
            } else {
                if ( bootstrapVCFOutput != null ) {
                    bootstrapVCFOutput.add(validation);
                }
                testSetSize++;
                return false;
            }
        } else {
            if ( validation != null && bootstrapVCFOutput != null ) {
                bootstrapVCFOutput.add(validation);
            }
            return false;
        }
    }

    private final static double[] HAPLOID_FLAT_LOG10_LIKELIHOODS = MathUtils.toLog10(new double[]{ 0.5, 0.0, 0.5 });
    private final static double[] DIPLOID_FLAT_LOG10_LIKELIHOODS = MathUtils.toLog10(new double[]{ 0.33, 0.33, 0.33 });

    public void writeBeagleOutput(VariantContext preferredVC, VariantContext otherVC, boolean isValidationSite, double prior) {
        GenomeLoc currentLoc = GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), preferredVC);
        StringBuffer beagleOut = new StringBuffer();

        String marker = String.format("%s:%d ",currentLoc.getContig(),currentLoc.getStart());
        beagleOut.append(marker);
        if ( markers != null ) markers.append(marker).append("\t").append(Integer.toString(markerCounter++)).append("\t");
        for ( Allele allele : preferredVC.getAlleles() ) {
            String bglPrintString;
            if (allele.isNoCall())
                bglPrintString = "-";
            else
                bglPrintString = allele.getBaseString();  // get rid of * in case of reference allele

            beagleOut.append(String.format("%s ", bglPrintString));
            if ( markers != null ) markers.append(bglPrintString).append("\t");
        }
        if ( markers != null ) markers.append("\n");

        GenotypesContext preferredGenotypes = preferredVC.getGenotypes();
        GenotypesContext otherGenotypes = goodSite(otherVC) ? otherVC.getGenotypes() : null;
        for ( String sample : samples ) {
            boolean isMaleOnChrX = CHECK_IS_MALE_ON_CHR_X && getSample(sample).getGender() == Gender.MALE;

            Genotype genotype;
            boolean isValidation;
            // use sample as key into genotypes structure
            if ( preferredGenotypes.containsSample(sample) ) {
                genotype = preferredGenotypes.get(sample);
                isValidation = isValidationSite;
            } else if ( otherGenotypes != null && otherGenotypes.containsSample(sample) ) {
                genotype = otherGenotypes.get(sample);
                isValidation = ! isValidationSite;
            } else {
                // there is magically no genotype for this sample.
                throw new StingException("Sample "+sample+" arose with no genotype in variant or validation VCF. This should never happen.");
            }

            /*
             * Use likelihoods if: is validation, prior is negative; or: is not validation, has genotype key
             */
            double [] log10Likelihoods = null;
            if ( (isValidation && prior < 0.0) || genotype.hasLikelihoods() ) {
                log10Likelihoods = genotype.getLikelihoods().getAsVector();

                // see if we need to randomly mask out genotype in this position.
                if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() <= insertedNoCallRate ) {
                    // we are masking out this genotype
                    log10Likelihoods = isMaleOnChrX ? HAPLOID_FLAT_LOG10_LIKELIHOODS : DIPLOID_FLAT_LOG10_LIKELIHOODS;
                }

                if( isMaleOnChrX ) {
                    log10Likelihoods[1] = -255;  // todo -- warning this is dangerous for multi-allele case
                }
            }
            /**
             * otherwise, use the prior uniformly
             */
            else if (! isValidation && genotype.isCalled() && ! genotype.hasLikelihoods() ) {
                // hack to deal with input VCFs with no genotype likelihoods.  Just assume the called genotype
                // is confident.  This is useful for Hapmap and 1KG release VCFs.
                double AA = (1.0-prior)/2.0;
                double AB = (1.0-prior)/2.0;
                double BB = (1.0-prior)/2.0;

                if (genotype.isHomRef()) { AA = prior; }
                else if (genotype.isHet()) { AB = prior; }
                else if (genotype.isHomVar()) { BB = prior; }

                log10Likelihoods = MathUtils.toLog10(new double[]{ AA, isMaleOnChrX ? 0.0 : AB, BB });
            }
            else  {
                log10Likelihoods = isMaleOnChrX ? HAPLOID_FLAT_LOG10_LIKELIHOODS : DIPLOID_FLAT_LOG10_LIKELIHOODS;
            }

            writeSampleLikelihoods(beagleOut, preferredVC, log10Likelihoods);
        }

        beagleWriter.println(beagleOut.toString());
    }

    private void writeSampleLikelihoods( StringBuffer out, VariantContext vc, double[] log10Likelihoods ) {
        if ( VQSRCalibrator != null ) {
            log10Likelihoods = VQSRCalibrator.includeErrorRateInLikelihoods(VQSLOD_KEY, vc, log10Likelihoods);
        }

        double[] normalizedLikelihoods = MathUtils.normalizeFromLog10(log10Likelihoods);
        // see if we need to randomly mask out genotype in this position.
        for (double likeVal: normalizedLikelihoods) {
            out.append(formatter.format(likeVal));
//            out.append(String.format("%5.4f ",likeVal));
        }
    }


    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return value + sum; // count up the sites
    }

    public void onTraversalDone( Integer includedSites ) {
        logger.info("Sites included in beagle likelihoods file             : " + includedSites);
        logger.info(String.format("Certain false positive found from recalibration curve : %d (%.2f%%)",
                certainFPs, (100.0 * certainFPs) / (Math.max(certainFPs + includedSites, 1))));
    }

    private void initializeVcfWriter() {
        final List<String> inputNames = Arrays.asList(validation.getName());

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));
        hInfo.add(new VCFFilterHeaderLine("bootstrap","This site used for genotype bootstrapping with ProduceBeagleInputWalker"));

        bootstrapVCFOutput.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));
    }

    public static class CachingFormatter {
        private String format;
        private LRUCache<Double, String> cache;

        public String getFormat() {
            return format;
        }

        public String format(double value) {
            String f = cache.get(value);
            if ( f == null ) {
                f = String.format(format, value);
                cache.put(value, f);
//                if ( cache.usedEntries() < maxCacheSize ) {
//                    System.out.printf("CACHE size %d%n", cache.usedEntries());
//                } else {
//                    System.out.printf("CACHE is full %f%n", value);
//                }
//            }
//            } else {
//                System.out.printf("CACHE hit %f%n", value);
//            }
            }

            return f;
        }

        public CachingFormatter(String format, int maxCacheSize) {
            this.format = format;
            this.cache = new LRUCache<Double, String>(maxCacheSize);
        }
    }

    /**
    * An LRU cache, based on <code>LinkedHashMap</code>.
    *
    * <p>
    * This cache has a fixed maximum number of elements (<code>cacheSize</code>).
    * If the cache is full and another entry is added, the LRU (least recently used) entry is dropped.
    *
    * <p>
    * This class is thread-safe. All methods of this class are synchronized.
    *
    * <p>
    * Author: Christian d'Heureuse, Inventec Informatik AG, Zurich, Switzerland<br>
    * Multi-licensed: EPL / LGPL / GPL / AL / BSD.
    */
    public static class LRUCache<K,V> {

    private static final float   hashTableLoadFactor = 0.75f;

    private LinkedHashMap<K,V>   map;
    private int                  cacheSize;

    /**
    * Creates a new LRU cache.
    * @param cacheSize the maximum number of entries that will be kept in this cache.
    */
    public LRUCache (int cacheSize) {
       this.cacheSize = cacheSize;
       int hashTableCapacity = (int)Math.ceil(cacheSize / hashTableLoadFactor) + 1;
       map = new LinkedHashMap<K,V>(hashTableCapacity, hashTableLoadFactor, true) {
          // (an anonymous inner class)
          private static final long serialVersionUID = 1;
          @Override protected boolean removeEldestEntry (Map.Entry<K,V> eldest) {
             return size() > LRUCache.this.cacheSize; }}; }

    /**
    * Retrieves an entry from the cache.<br>
    * The retrieved entry becomes the MRU (most recently used) entry.
    * @param key the key whose associated value is to be returned.
    * @return    the value associated to this key, or null if no value with this key exists in the cache.
    */
    public synchronized V get (K key) {
       return map.get(key); }

    /**
    * Adds an entry to this cache.
    * The new entry becomes the MRU (most recently used) entry.
    * If an entry with the specified key already exists in the cache, it is replaced by the new entry.
    * If the cache is full, the LRU (least recently used) entry is removed from the cache.
    * @param key    the key with which the specified value is to be associated.
    * @param value  a value to be associated with the specified key.
    */
    public synchronized void put (K key, V value) {
       map.put (key, value); }

    /**
    * Clears the cache.
    */
    public synchronized void clear() {
       map.clear(); }

    /**
    * Returns the number of used entries in the cache.
    * @return the number of entries currently in the cache.
    */
    public synchronized int usedEntries() {
       return map.size(); }

    /**
    * Returns a <code>Collection</code> that contains a copy of all cache entries.
    * @return a <code>Collection</code> with a copy of the cache content.
    */
    public synchronized Collection<Map.Entry<K,V>> getAll() {
       return new ArrayList<Map.Entry<K,V>>(map.entrySet()); }

    } // end class LRUCache
}
