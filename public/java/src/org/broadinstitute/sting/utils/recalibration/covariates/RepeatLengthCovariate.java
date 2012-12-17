package org.broadinstitute.sting.utils.recalibration.covariates;

import org.broadinstitute.sting.gatk.walkers.bqsr.RecalibrationArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.TandemRepeat;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 11/3/12
 */

public class RepeatLengthCovariate implements ExperimentalCovariate {
    final int MAX_REPEAT_LENGTH = 20;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {}

    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {
        byte[] readBytes = read.getReadBases();
        for (int i = 0; i < readBytes.length; i++) {
            int maxRL = 0;
            for (int str = 1; str <= 8; str++) {
                if (i + str <= readBytes.length) {
                    maxRL = Math.max(maxRL, VariantContextUtils.findNumberofRepetitions(
                            Arrays.copyOfRange(readBytes,i,i + str),
                            Arrays.copyOfRange(readBytes,i,readBytes.length)
                    ));
                }
            }
            if(maxRL > MAX_REPEAT_LENGTH) { maxRL = MAX_REPEAT_LENGTH; }
            values.addCovariate(maxRL, maxRL, maxRL, i);
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Byte.parseByte(str);
    }

    @Override
    public String formatKey(final int key) {
        return String.format("%d", key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? Integer.parseInt((String) value) : (Integer) value;
    }

    @Override
    public int maximumKeyValue() {
        return MAX_REPEAT_LENGTH + 1;
    }

}
