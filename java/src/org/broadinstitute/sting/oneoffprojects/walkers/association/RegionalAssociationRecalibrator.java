package org.broadinstitute.sting.oneoffprojects.walkers.association;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.table.BedTableCodec;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/25/11
 * Time: 11:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class RegionalAssociationRecalibrator extends RodWalker<RegionalAssociationRecalibrator.LocHashingPair,Set<RegionalAssociationRecalibrator.LocHashingPair>> implements TreeReducible<Set<RegionalAssociationRecalibrator.LocHashingPair>> {

    @Argument(shortName="w",fullName="whiteningFile",doc="File containing the mean, max vectors and covariance matrix",required=false)
    File whitenFile = null;
    @Argument(shortName="p",fullName="normParameter",doc="Exponent for vector norm; p > 4 will induce a p=infinity approximation",required=false)
    double pNorm = 1.0;

    @Output
    public PrintStream out;

    private DataWhitener whitener;
    private double[] dataHolder;
    private int boundTables = 0;

    public void initialize() {

        for ( ReferenceOrderedDataSource source : getToolkit().getRodDataSources() ) {
            logger.debug(source.getType().getSimpleName());
            if ( source.getType().equals(BedTableCodec.class)) {
                ++boundTables;
                if ( ! source.getFile().getName().endsWith(".bedgraph") ) {
                    throw new UserException("Regional association requires bedgraph files. The file "+source.getFile().getAbsolutePath()+" does not have the proper extension.");
                }
            }
        }

        if ( whitenFile != null ) {
            whitener = getWhitenerFromFile(whitenFile);
        }
    }

    public Set<LocHashingPair> reduceInit() { return new TreeSet<LocHashingPair>(); }

    public LocHashingPair map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext ali ) {
        if ( tracker == null ) { return null; }

        ArrayList<TableFeature> features = new ArrayList<TableFeature>(boundTables);

        for ( GATKFeature feature : tracker.getAllRods()) {
            Object uo = feature.getUnderlyingObject();
            if ( uo instanceof TableFeature && feature.getLocation().getStart() == ref.getLocus().getStart()) {
                features.add((TableFeature)uo);
            }
        }

        if ( features.size() == 0 ) { return null; }
        if ( features.size() != boundTables ) {
            throw new UserException("Features do not line up at position "+ref.getLocus().toString());
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

    public Set<LocHashingPair> treeReduce(Set<LocHashingPair> left, Set<LocHashingPair> right) {
        if ( left == null ) {
            return right;
        } else if ( right == null ) {
            return left;
        } else {
            left.addAll(right);
            return left;
        }
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
            // note - equal values will always be assigned the same rank -- however
            // it is proper, no need for random assignment.
            double prob = (1.0 + normRanks.length - rank)/(1.0+normRanks.length);
            int qual = AssociationTestRunner.pToQ(prob);
            out.printf("%s\t%d\t%d\t%d\t%.2e\t%d\t%.2e%n",lhp.getFirst().getContig(),lhp.getFirst().getStart(),lhp.getFirst().getStop(),
                    qual,prob,rank,lhp.getSecond());
        }
    }

    private DataWhitener getWhitenerFromFile(File file){
        XReadLines lineReader;
        try {
            lineReader = new XReadLines(file);
        } catch (FileNotFoundException e) {
            throw new UserException("The provided whiten file could not be found",e);
        }

        // first line is the mean vector
        String[] meanLine = lineReader.next().split("\t");
        if ( meanLine.length != boundTables ) {
            String msg = String.format("The number of bound bedgraphs does not match the size of the mean value in the whitening file. bound: %d mean : %d",boundTables,meanLine.length);
            throw new UserException(msg);
        }
        double[] mean = new double[meanLine.length];
        for ( int i = 0 ; i < meanLine.length; i ++ ) {
            mean[i] = Double.parseDouble(meanLine[i]);
        }
        // skip
        lineReader.next();
        // now the max
        String[] maxLine = lineReader.next().split("\t");
        double[] max = new double[maxLine.length];
        for ( int i = 0; i < maxLine.length ; i++ ) {
            max[i] = Double.parseDouble(maxLine[i]);
        }
        // skip
        lineReader.next();
        // now the covariance matrix
        double[][] cov = new double[max.length][max.length];
        for ( int i = 0; i < max.length; i ++ ) {
            String[] row = lineReader.next().split("\t");
            for ( int j = 0; j < max.length; j ++ ) {
                cov[i][j] = Double.parseDouble(row[j]);
            }
        }

        DataWhitener dataWhitener = new DataWhitener();
        dataWhitener.setCov(cov);
        dataWhitener.setMaxs(max);
        dataWhitener.setMeans(mean);

        return dataWhitener;
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
            // note: want vectors that have NaN to have norm calculated within the
            //    n-k dimensional subspace. This is equivalent to setting the NaNs to
            //    a post-whitening value of zero.
            @Override
            public double apply(double v, double v1) {
                if ( Double.compare(v,Double.NaN) == 0 || Double.compare(v1,Double.NaN) == 0 ) {
                    return 0.0;
                }

		if ( Double.compare(v,Double.POSITIVE_INFINITY) == 0 ) { return 50.0; } // todo -- fixme
		if ( Double.compare(v,Double.NEGATIVE_INFINITY) == 0 ) { return -50.0; } // todo -- fixme, really should be an exception
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
            for ( int i = 0; i < diag.rows(); i ++ ) {
                // do not artificially inflate signal that is not there
                if ( diag.get(i,i) < 1.0 ) {
                    diag.set(i,i, 1.0);
                } else {
                    diag.set(i,i, Math.pow(diag.get(i,i),-0.5));
                }
            }
            transform = ALGEBRA.mult(diag,ALGEBRA.transpose(decomp.getV()));
            logger.debug("TRANSFORM:");
            logger.debug(transform);
        }

        public double[] whiten(double[] value) {
            DenseDoubleMatrix1D vec = new DenseDoubleMatrix1D(value);
            return ALGEBRA.mult(transform, vec.assign(means, SUBTRACT) ).toArray();
        }
    }

    class LocHashingPair extends Pair<GenomeLoc,Double> implements Comparable {

        public LocHashingPair(GenomeLoc g, Double d) {
            super(g,d);
        }

        public int hashCode() { return first.hashCode(); }

        public int compareTo(Object o) {
            if ( o instanceof LocHashingPair ) {
                return this.first.compareTo(((LocHashingPair) o).first);
            } else {
                return Integer.MIN_VALUE;
            }
        }
    }
}
