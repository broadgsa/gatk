package org.broadinstitute.sting.gatk.datasources.sample;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Aug 13, 2010
 * Time: 5:13:46 PM
 */
public class SampleAlias {

    String mainId;

    String[] otherIds;

    public String getMainId() {
        return mainId;
    }

    public void setMainId(String mainId) {
        this.mainId = mainId;
    }

    public String[] getOtherIds() {
        return otherIds;
    }

    public void setOtherIds(String[] otherIds) {
        this.otherIds = otherIds;
    }

}
