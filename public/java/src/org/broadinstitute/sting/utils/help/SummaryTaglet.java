package org.broadinstitute.sting.utils.help;

import com.sun.tools.doclets.Taglet;
import com.sun.javadoc.Tag;

import java.util.Map;

/**
 * Provide an alternate brief summary for this walker / package.
 * Acts as an alternative to the first sentence employed by default.
 * @author mhanna
 * @version 0.1
 */
public class SummaryTaglet extends HelpTaglet {
    /**
     * The key tag for this taglet.
     */
    public static final String NAME = "help.summary";

    /**
     * Return the name of this custom tag.
     */
    @Override
    public String getName() {
        return NAME;
    }

    /**
     * Will return false since overviews are always named
     * by the <code>@WalkerName</code> tag.
     * @return false always
     */
    @Override
    public boolean inOverview() {
        return true;
    }

    /**
     * Will return true to indicate that packages can be given useful summary.
     * @return true always
     */
    @Override
    public boolean inPackage() {
        return true;
    }

    /**
     * Register this Taglet.
     * @param tagletMap  the map to register this tag to.
     */
    public static void register(Map tagletMap) {
       SummaryTaglet tag = new SummaryTaglet();
       Taglet t = (Taglet)tagletMap.get(tag.getName());
       if (t != null) {
           tagletMap.remove(tag.getName());
       }
       tagletMap.put(tag.getName(), tag);
    }
}