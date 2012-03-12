package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.commons.lang.NotImplementedException;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.security.cert.CertificateNotYetValidException;
import java.util.*;

import org.broadinstitute.sting.utils.codecs.vcf.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 8/30/11
 * Time: 10:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class UGBoundAF extends RodWalker<VariantContext,Integer> {

    @Output(shortName="vcf",fullName="VCF",doc="file to write to",required=true)
    VCFWriter writer;

    @Input(shortName="V",fullName="Variants",doc="variant tracks to use in calculation",required=true)
    List<RodBinding<VariantContext>> variants;

    private static double EPS_LOWER_LIMIT = Math.pow(10,-6.0);

    private HashMap<Integer,Pair<Double,Double>> epsilonPosteriorCache = new HashMap<Integer,Pair<Double,Double>>(8192);
    private HashMap<Integer,Double> logAC0Cache = new HashMap<Integer, Double>(8192);
    private int QUANTIZATION_FACTOR = 1000;


    public void initialize() {
        Set<VCFHeaderLine> allHeaderLines = new HashSet<VCFHeaderLine>(1024);
        for ( RodBinding<VariantContext> v : variants ) {
            String trackName = v.getName();
            Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));
            Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>(vcfHeaders.get(trackName).getMetaData());
        }
        allHeaderLines.add(new VCFInfoHeaderLine("AFB",2,VCFHeaderLineType.Float,"The 95% bounds on the allele "+
                "frequency. First value=95% probability AF>x. Second value=95% probability AF<x."));
        writer.writeHeader(new VCFHeader(allHeaderLines));
    }
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext unused ) {
        List<VariantContext> allVariants = tracker.getValues(variants);
        if ( allVariants.size() == 0 ) {
            return null;
        }

        List<Allele> alternateAlleles = getAllAlternateAlleles(allVariants);
        VariantContextBuilder builder = new VariantContextBuilder(allVariants.get(0).subContextFromSamples(new TreeSet<String>()));
        if ( alternateAlleles.size() > 1 ) {
            logger.warn("Multiple Segregating Variants at position "+ref.getLocus().toString());
            alternateAlleles.add(allVariants.get(0).getReference());
            builder.alleles(alternateAlleles);
            builder.filters(String.format("MULTIPLE_SEGREGATING[%s]", Utils.join(",",alternateAlleles)));
        } else {
            // get all the genotype likelihoods
            GenotypesContext context = GenotypesContext.create();
            int numNoCall = 0;
            for ( VariantContext v : allVariants ) {
                numNoCall += v.getNoCallCount();
                context.addAll(v.getGenotypes());
            }
            builder.attribute("AFB",boundAlleleFrequency(getACPosteriors(context)));
        }

        return builder.make();
    }

    private List<Allele> getAllAlternateAlleles(List<VariantContext> variants) {
        List<Allele> alleles = new ArrayList<Allele>(3); // some overhead
        for ( VariantContext v : variants ) {
            alleles.addAll(v.getAlternateAlleles());
        }
        return alleles;
    }

    @Override
    public Integer reduce(VariantContext value, Integer sum) {
        if ( value == null )
            return sum;
        writer.add(value);
        return ++sum;
    }

    private int N_ITERATIONS = 1;
    private double[] getACPosteriors(GenotypesContext gc) {
        // note this uses uniform priors (!)

        double[][] zeroPriors = new double[1][1+2*gc.size()];
        AlleleFrequencyCalculationResult result = new AlleleFrequencyCalculationResult(2,2*gc.size());
        // todo -- allow multiple alleles here
        for ( int i = 0; i < N_ITERATIONS; i ++ ) {
            ExactAFCalculationModel.linearExactMultiAllelic(gc, 2, zeroPriors, result, false);
        }

        return result.log10AlleleFrequencyPosteriors[0];
    }

    private String boundAlleleFrequency(double[] ACLikelihoods) {
        // note that no-calls are unnecessary: the ML likelihoods take nocalls into account as 0,0,0 GLs
        // thus, for sites with K 100,40,0 likelihoods and M no-calls, the likelihoods will be
        // agnostic between 2*K alleles through 2*(K+M) alleles - exactly what we want to marginalize over

        // want to pick a lower limit x and upper limit y such that
        // int_{f = x to y} sum_{c = 0 to 2*AN} P(AF=f | c, AN) df = 0.95
        // int_{f=x to y} calculateAFPosterior(f) df = 0.95
        // and that (y-x) is minimized

        // this is done by quantizing [0,1] into small bins and, since the distribution is
        // unimodal, greedily adding them until the probability is >= 0.95

        throw new ReviewedStingException("This walker is unsupported, and is not fully implemented", new NotImplementedException("bound allele frequency not implemented"));
    }

    private double calculateAFPosterior(double[] likelihoods, double af) {
        double[] probLiks = new double[likelihoods.length];
        for ( int c = 0; c < likelihoods.length; c++) {
            probLiks[c] = calculateAFPosterior(c,likelihoods.length,af);
        }

        return MathUtils.log10sumLog10(probLiks);
    }

    private double calculateAFPosterior(int ac, int n, double af) {
        // evaluate the allele frequency posterior distribution at AF given AC observations of N chromosomes
        switch ( ac ) {
            case 0:
                return logAC0Coef(n) + n*Math.log10(1 - af) - Math.log10(af);
            case 1:
                return Math.log10(n) + (n-1)*Math.log10(1-af) - n*Math.log10(1-EPS_LOWER_LIMIT);
            case 2:
                return Math.log10(n) + Math.log10(n-1) + Math.log10(af) + (n-2)*Math.log10(1-af) - Math.log10(1-(n-1)*EPS_LOWER_LIMIT) - (n-1)*Math.log10(EPS_LOWER_LIMIT);
            default:
                return  (ac-1)*Math.log10(af)+ac*Math.log10( (double) n-ac)-(n-ac)*af*Math.log10(Math.E) - MathUtils.log10Gamma(ac);
        }
    }

    private double logAC0Coef(int an) {
        if ( ! logAC0Cache.containsKey(an) ) {
            double coef = -Math.log10(EPS_LOWER_LIMIT);
            for ( int k = 1; k <= an; k++ ) {
                // note this should typically just be
                // term = ( 1 - Math.pow(EPS_LOWER_LIMIT,k) ) * MathUtils.binomialCoefficient(an,k) / k
                // but the 1-E term will just be 1, so we do the following to mitigate this problem
                double binom = MathUtils.binomialCoefficient(an,k);
                double eps_correction = EPS_LOWER_LIMIT*Math.pow(binom,1/k);
                double term = binom/k - Math.pow(eps_correction,k);
                if ( k % 2 == 0 ) {
                    coef += term;
                }  else {
                    coef -= term;
                }
            }

            logAC0Cache.put(an,coef);
        }

        return logAC0Cache.get(an);
    }

    private double adaptiveSimpson(double[] likelihoods, double start, double stop, double err, int cap) {
        double mid = (start + stop)/2;
        double size = stop-start;
        double fa = calculateAFPosterior(likelihoods,start);
        double fb = calculateAFPosterior(likelihoods,mid);
        double fc = calculateAFPosterior(likelihoods,stop);
        double s = (size/6)*(fa + 4*fc + fb);
        double h = simpAux(likelihoods,start,stop,err,s,fa,fb,fc,cap);
        return h;
    }

    private double simpAux(double[] likelihoods, double a,double b,double eps,double s,double fa,double fb,double fc,double cap){
        if ( s == 0 )
                return -300.0;
        double c = ( a + b )/2;
        double h = b-a;
        double d = (a + c)/2;
        double e = (c + b)/2;
        double fd = calculateAFPosterior(likelihoods, d);
        double fe = calculateAFPosterior(likelihoods, e);
        double s_l = (h/12)*(fa + 4*fd + fc);
        double s_r = (h/12)*(fc + 4*fe + fb);
        double s_2 = s_l + s_r;
        if ( cap <= 0 || Math.abs(s_2 - s) <= 15*eps ){
            return Math.log10(s_2 + (s_2 - s)/15.0);
        }

        return MathUtils.approximateLog10SumLog10(simpAux(likelihoods,a,c,eps/2,s_l,fa,fc,fd,cap-1),simpAux(likelihoods, c, b, eps / 2, s_r, fc, fb, fe, cap - 1));
    }
}
