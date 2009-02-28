package edu.mit.broad.picard.metrics;

import edu.mit.broad.picard.util.StringUtil;

/**
 * Header that stores information about the version of some piece of software or
 * data used to create the metrics file.  Payload consists of a name or description
 * of the versioned item and a version string.
 *
 * @author Tim Fennell
 */
public class VersionHeader implements Header {
    private String versionedItem;
    private String versionString;

    public void parse(String in) {
        String[] fields = in.split("\t");
        this.versionedItem = fields[0];
        this.versionString = fields[1];
    }

    public String toString() {
        return this.versionedItem + "\t" + this.versionString;
    }

    public String getVersionedItem() { return versionedItem; }
    public void setVersionedItem(String versionedItem) {
        this.versionedItem = StringUtil.assertCharactersNotInString(versionedItem, '\t', '\n');
    }

    public String getVersionString() { return versionString; }
    public void setVersionString(String versionString) {
        this.versionString = StringUtil.assertCharactersNotInString(versionString, '\t', '\n');
    }

    /** Equals method that checks that both the item and version string are equal. */
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VersionHeader that = (VersionHeader) o;

        if (versionString != null ? !versionString.equals(that.versionString) : that.versionString != null)
            return false;
        if (versionedItem != null ? !versionedItem.equals(that.versionedItem) : that.versionedItem != null)
            return false;

        return true;
    }
}
