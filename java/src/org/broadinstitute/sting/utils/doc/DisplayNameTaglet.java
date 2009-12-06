package org.broadinstitute.sting.utils.doc;

import com.sun.tools.doclets.Taglet;
import com.sun.javadoc.Tag;

import java.util.Map;

/**
 * Provide a display name in the help for packages
 *
 * @author mhanna
 * @version 0.1
 */
public class DisplayNameTaglet implements Taglet {
    /**
     * The display name for this taglet.
     */
    public static final String NAME = "display.name";

    /**
     * Return the name of this custom tag.
     */
    public String getName() {
        return NAME;
    }

    /**
     * Will return false since this tag cannot be applied
     * to a field.
     * @return false since this tag cannot be applied to a field.
     */
    public boolean inField() {
        return false;
    }

    /**
     * Will return false since this tag cannot be applied
     * to a constructor.
     * @return false since this tag cannot be applied to a constructor.
     */
    public boolean inConstructor() {
        return false;
    }

    /**
     * Will return false since this tag cannot be applied
     * to a method.
     * @return false since this tag cannot be applied to a method.
     */
    public boolean inMethod() {
        return false;
    }

    /**
     * Will return false since overviews are always named
     * by the <code>@WalkerName</code> tag.
     * @return false always
     */
    public boolean inOverview() {
        return false;
    }

    /**
     * Will return true to indicate that packages can be given useful
     * display text.
     * @return true always
     */
    public boolean inPackage() {
        return true;
    }

    /**
     * Will return false indicating that types cannot be given
     * alternate display names.
     * @return false always.
     */
    public boolean inType() {
        return false;
    }

    /**
     * Will return false since <code>@todo</code>
     * is not an inline tag.
     * @return false since <code>@todo</code>
     * is not an inline tag.
     */

    public boolean isInlineTag() {
        return false;
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

    /**
     * Create a string representation of this tag.  Since this tag is only
     * used by the help system, don't output any HTML.
     */
    public String toString(Tag tag) {
        return null;
    }

    /**
     * Create a string representation of this tag.  Since this tag is only
     * used by the help system, don't output any HTML.
     */
    public String toString(Tag[] tags) {
        return null;
    }
}
