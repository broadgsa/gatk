package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.util.ArrayList;
import java.util.List;
import org.broadinstitute.sting.oneoffprojects.walkers.association.interfaces.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 11:48:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class AssociationTestRunner {

    public static List<String> runTests(AssociationContext context) {
        List<String> results = new ArrayList<String>();
        if ( context.getClass().isInstance(TStatistic.class)) {
            results.add(runStudentT(context));
        }

        if ( context.getClass().isInstance(ZStatistic.class)) {
            results.add(runZ(context));
        }

        if ( context.getClass().isInstance(FisherExact.class) ) {
            results.add(runFisherExact(context));
        }

        return results;
    }

    public static String runStudentT(AssociationContext context) {
        return "Test not yet implemented";
    }

    public static String runZ(AssociationContext context) {
        return "Test not yet implemented";
    }

    public static String runFisherExact(AssociationContext context) {
        return "Test not yet implemented";
    }
}
