package org.broadinstitute.sting.oneoffprojects.walkers.association;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.bed.BedParser;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/25/11
 * Time: 11:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class RegionalAssociationRecalibrator extends RodWalker<RegionalAssociationRecalibrator.LocHashingPair,Set<RegionalAssociationRecalibrator.LocHashingPair>> {

    @Argument(shortName="w",fullName="whiteningFile",doc="File containing the mean, max vectors and covariance matrix",required=false)
    File whitenFile = null;
    @Argument(shortName="p",fullName="normParameter",doc="Exponent for vector norm; p > 4 will induce a p=infinity approximation",required=false)
    double pNorm = 1.0;
    @Hidden
    @Argument(shortName="cg",fullName="capacityGuess",doc="this is hidden",required = false)
    int initialCapacity = 3500000;

    @Output
    public PrintStream out;

    private DataWhitener whitener;
    private double[] dataHolder;

    public void initialize() {
        if ( whitenFile != null ) {
            whitener = getWhitenerFromFile(whitenFile);
        }
    }

    public Set<LocHashingPair> reduceInit() { return new HashSet<LocHashingPair>(initialCapacity); }

    public LocHashingPair map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext ali ) {
        if ( tracker == null ) { return null; }

        ArrayList<TableFeature> features = new ArrayList<TableFeature>(7);

        for ( GATKFeature feature : tracker.getAllRods() ) {
            Object uo = feature.getUnderlyingObject();
            if ( uo instanceof TableFeature ) {
                features.add((TableFeature)uo);
            }
        }

        if ( dataHolder == null ) {
            dataHolder = new double[features.size()];
        }

        int idx = 0;
        for ( TableFeature f : features ) {
            dataHolder[idx] = Double.parseDouble(f.getValue(3));
            idx++;
        }

        double normVal;

        if ( whitener != null ) {
            normVal = calculateNorm(whitener.whiten(dataHolder),pNorm);
        } else {
            normVal = calculateNorm(dataHolder,pNorm);
        }

        return new LocHashingPair(features.get(0).getLocation(),normVal);
    }

    public Set<LocHashingPair> reduce(LocHashingPair map, Set<LocHashingPair> red) {
        if ( map != null ) { red.add(map); }
        return red;
    }

    public void onTraversalDone(Set<LocHashingPair> normByLoc) {
        double[] normRanks = new double[(normByLoc.size())];
        int offset = 0;
        for ( LocHashingPair lhp : normByLoc ) {
            normRanks[offset]=(lhp.second);
            offset++;
        }
        Arrays.sort(normRanks);
        for ( LocHashingPair lhp : normByLoc ) {
            int rank = Arrays.binarySearch(normRanks,lhp.second);
            double prob = ((double) rank)/((double)normRanks.length);
            int qual = AssociationTestRunner.pToQ(prob);
            out.printf("%s\t%d\t%d\t%d\t%.2e\t%d",lhp.getFirst().getContig(),lhp.getFirst().getStart(),lhp.getFirst().getStop(),
                    qual,prob,rank);
        }
    }

    private static DataWhitener getWhitenerFromFile(File file) {
        // todo -- implement me
        return null;
    }

    private static double calculateNorm(double[] vec, double p) {
        double sumExp = 0.0;
        for ( double v : vec ) {
            sumExp += Math.pow(Math.abs(v),p);
        }

        return Math.pow(sumExp,1.0/p);
    }

    class DataWhitener {

        private final Algebra ALGEBRA = new Algebra();
        private final DoubleDoubleFunction SUBTRACT = new DoubleDoubleFunction() {
            @Override
            public double apply(double v, double v1) {
                return v - v1;
            }
        };

        private DenseDoubleMatrix1D means;
        private DenseDoubleMatrix1D maxAbs;
        private DoubleMatrix2D transform;

        public DataWhitener() {
        }

        public void setMeans(double[] inMeans) {
            means = new DenseDoubleMatrix1D(inMeans);
        }

        public void setMaxs(double[] maxs) {
            maxAbs = new DenseDoubleMatrix1D(maxs);
        }

        public void setCov(double[][] cov) {
            EigenvalueDecomposition decomp = new EigenvalueDecomposition(new DenseDoubleMatrix2D(cov));
            DoubleMatrix2D diag = decomp.getD();
            for ( int i = 0; i < diag.size(); i ++ ) {
                diag.set(i,i, Math.pow(diag.get(i,i),-0.5));
            }
            transform = ALGEBRA.mult(diag,ALGEBRA.transpose(decomp.getV()));
        }

        public double[] whiten(double[] value) {
            DenseDoubleMatrix1D vec = new DenseDoubleMatrix1D(value);
            return ALGEBRA.mult(transform, vec.assign(means, SUBTRACT) ).toArray();
        }
    }

    class LocHashingPair extends Pair<GenomeLoc,Double> {

        public LocHashingPair(GenomeLoc g, Double d) {
            super(g,d);
        }

        public int hashCode() { return first.hashCode(); }
    }
}
