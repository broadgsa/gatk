package org.broadinstitute.sting.utils.help;

import com.sun.tools.doclets.Taglet;
import com.sun.javadoc.Tag;

import java.util.Map;

/**
 * Basic functionality for the help taglet.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class HelpTaglet implements Taglet {
    /**
     * Return the name of this custom tag.
     */
    public abstract String getName();

    /**
     * Will return false since this tag cannot be applied
     * to a field.
     * @return false since this tag cannot be applied to a field.
     */
    public boolean inField() {
        return false;
    }

    /**
     * Will return false since by default, help tags cannot be applied to a constructor.
     * @return false since by default, help tags cannot be applied to a constructor.
     */
    public boolean inConstructor() {
        return false;
    }

    /**
     * Will return false since by default, help tags cannot be applied to a method.
     * @return false since by default, this tag cannot be applied to a method.
     */
    public boolean inMethod() {
        return false;
    }

    /**
     * Will return false since by default, help tags cannot be applied to an overview.
     * @return false since by default, help tags cannot be applied to an overview.
     */
    public boolean inOverview() {
        return false;
    }

    /**
     * Will return false since by default, help tags cannot be applied to a package.
     * description.
     * @return false since by default, help tags cannot be applied to a package.
     */
    public boolean inPackage() {
        return false;
    }

    /**
     * Will return false since help tags are by default not inline.
     * @return false since help tags are by default not inline.
     */
    public boolean inType() {
        return false;
    }

    /**
     * Will return false since help tags are by default not inline.
     * @return false since help tags are by default not inline.
     */
    public boolean isInlineTag() {
        return false;
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
