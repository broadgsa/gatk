package org.broadinstitute.sting.utils.help;

import com.sun.tools.doclets.Taglet;
import com.sun.javadoc.Tag;

import java.util.Map;

/**
 * Provide a display name in the help for packages
 *
 * @author mhanna
 * @version 0.1
 */
public class DisplayNameTaglet extends HelpTaglet {
    /**
     * The display name for this taglet.
     */
    public static final String NAME = "help.display.name";

    /**
     * Return the name of this custom tag.
     */
    @Override
    public String getName() {
        return NAME;
    }

    /**
     * Will return true to indicate that packages can be given useful
     * display text.
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
       DisplayNameTaglet tag = new DisplayNameTaglet();
       Taglet t = (Taglet)tagletMap.get(tag.getName());
       if (t != null) {
           tagletMap.remove(tag.getName());
       }
       tagletMap.put(tag.getName(), tag);
    }
}
