package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * User: depristo
 * Date: 8/3/12
 * Time: 12:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class AdaptiveContextUnitTest {
    // TODO
    // TODO actually need unit tests when we have validated the value of this approach
    // TODO particularly before we attempt to optimize the algorithm
    // TODO

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class AdaptiveContextTestProvider extends BaseTest.TestDataProvider {
        final RecalDatumNode<ContextDatum> pruned;
        final RecalDatumNode<ContextDatum> full;

        private AdaptiveContextTestProvider(Class c, RecalDatumNode<ContextDatum> pruned, RecalDatumNode<ContextDatum> full) {
            super(AdaptiveContextTestProvider.class);
            this.pruned = pruned;
            this.full = full;
        }
    }

    private RecalDatumNode<ContextDatum> makeTree(final String context, final int N, final int M,
                                                  final RecalDatumNode<ContextDatum> ... sub) {
        final ContextDatum contextDatum = new ContextDatum(context, N, M);
        final RecalDatumNode<ContextDatum> node = new RecalDatumNode<ContextDatum>(contextDatum);
        for ( final RecalDatumNode<ContextDatum> sub1 : sub ) {
            node.addSubnode(sub1);
        }
        return node;
    }

    @DataProvider(name = "AdaptiveContextTestProvider")
    public Object[][] makeRecalDatumTestProvider() {
//        final RecalDatumNode<ContextDatum> prune1 =
//                makeTree("A", 10, 1,
//                        makeTree("AA", 11, 2),
//                        makeTree("AC", 12, 3),
//                        makeTree("AG", 13, 4),
//                        makeTree("AT", 14, 5));
//
//        new AdaptiveContextTestProvider(pruned, full);

        return AdaptiveContextTestProvider.getTests(AdaptiveContextTestProvider.class);
    }

    @Test(dataProvider = "AdaptiveContextTestProvider")
    public void testAdaptiveContextFill(AdaptiveContextTestProvider cfg) {

    }
}
