package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.MathUtils;

import java.util.LinkedList;
import java.util.List;

/**
 * Short one line description of the walker.
 *
 * <p> [Long description of the walker] </p>
 *
 *
 * <h2>Input</h2> <p> [Description of the Input] </p>
 *
 * <h2>Output</h2> <p> [Description of the Output] </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T [walker name]
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 4/17/12
 */
public class AccuracyDatum extends RecalDatum {
    private final List<Double> accuracy = new LinkedList<Double>();
    private final List<Byte> reportedQualities = new LinkedList<Byte>();

    public AccuracyDatum(final RecalDatum recalDatum, final byte originalQuality) {
        super(recalDatum);
        accuracy.add(calculateAccuracy(recalDatum, originalQuality));
        reportedQualities.add(originalQuality);
    }

    public void combine(final RecalDatum recalDatum, final byte originalQuality) {
        this.combine(recalDatum);
        accuracy.add(calculateAccuracy(recalDatum, originalQuality));
        reportedQualities.add(originalQuality);
    }

    @Override
    public String toString() {
        return String.format("%s,%.2f,%.2f", super.toString(), MathUtils.average(reportedQualities), MathUtils.average(accuracy));
    }

    private static double calculateAccuracy(final RecalDatum recalDatum, final byte originalQuality) {
        return recalDatum.getEmpiricalQuality() - originalQuality;
    }
}
