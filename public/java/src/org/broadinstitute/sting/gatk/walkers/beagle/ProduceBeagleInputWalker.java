/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broadinstitute.sting.commandline.*;
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
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

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

public class ProduceBeagleInputWalker extends RodWalker<Integer, Integer> {

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
    VCFWriter bootstrapVCFOutput = null;

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
                    bootstrapVCFOutput.add(VariantContext.modifyFilters(validation, BOOTSTRAP_FILTER));
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
        GenomeLoc currentLoc = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),preferredVC);
        StringBuffer beagleOut = new StringBuffer();

        String marker = String.format("%s:%d ",currentLoc.getContig(),currentLoc.getStart());
        beagleOut.append(marker);
        if ( markers != null ) markers.append(marker).append("\t").append(Integer.toString(markerCounter++)).append("\t");
        for ( Allele allele : preferredVC.getAlleles() ) {
            String bglPrintString;
            if (allele.isNoCall() || allele.isNull())
                bglPrintString = "-";
            else
                bglPrintString = allele.getBaseString();  // get rid of * in case of reference allele

            beagleOut.append(String.format("%s ", bglPrintString));
            if ( markers != null ) markers.append(bglPrintString).append("\t");
        }
        if ( markers != null ) markers.append("\n");

        Map<String,Genotype> preferredGenotypes = preferredVC.getGenotypes();
        Map<String,Genotype> otherGenotypes = goodSite(otherVC) ? otherVC.getGenotypes() : null;
        for ( String sample : samples ) {
            boolean isMaleOnChrX = CHECK_IS_MALE_ON_CHR_X && getSample(sample).getGender() == Gender.MALE;

            Genotype genotype;
            boolean isValidation;
            // use sample as key into genotypes structure
            if ( preferredGenotypes.keySet().contains(sample) ) {
                genotype = preferredGenotypes.get(sample);
                isValidation = isValidationSite;
            } else if ( otherGenotypes != null && otherGenotypes.keySet().contains(sample) ) {
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
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), inputNames));
        hInfo.add(new VCFFilterHeaderLine("bootstrap","This site used for genotype bootstrapping with ProduceBeagleInputWalker"));

        bootstrapVCFOutput.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));
    }

    public static class CachingFormatter {
        private int maxCacheSize = 0;
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
            this.maxCacheSize = maxCacheSize;
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
