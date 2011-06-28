package org.broadinstitute.sting;

import org.testng.reporters.TextReporter;

/**
 * HACK: Create a variant of the TestNG TextReporter that can be run with no
 *       arguments, and can therefore be added to the TestNG listener list.
 *
 * @author hanna
 * @version 0.1
 */
public class StingTextReporter extends TextReporter {
    public StingTextReporter() {
        super("Ant suite",2);
    }
}
