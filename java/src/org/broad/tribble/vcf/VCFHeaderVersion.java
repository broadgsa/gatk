package org.broad.tribble.vcf;

/**
 * information that identifies each header version
 */
public enum VCFHeaderVersion {
    VCF3_2("VCRv3.2","format"),
    VCF3_3("VCFv3.3","fileformat"),
    VCF4_0("VCFv4.0","fileformat");

    private final String versionString;
    private final String formatString;

    /**
     * create the enum, privately, using:
     * @param vString the version string
     * @param fString the format string
     */
    VCFHeaderVersion(String vString, String fString) {
        this.versionString = vString;
        this.formatString = fString;
    }

    /**
     * get the header version
     * @param version the version string
     * @param format the format string
     * @return a VCFHeaderVersion object
     */
    public static VCFHeaderVersion toHeaderVersion(String version, String format) {
        for (VCFHeaderVersion hv : VCFHeaderVersion.values())
            if (hv.versionString.equals(version) && hv.formatString.equals(format))
                    return hv;
        return null;
    }

    /**
     * get the header version
     * @param version the version string
     * @return a VCFHeaderVersion object
     */
    public static VCFHeaderVersion toHeaderVersion(String version) {
        for (VCFHeaderVersion hv : VCFHeaderVersion.values())
            if (hv.versionString.equals(version))
                    return hv;
        return null;
    }

    /**
     * are we a valid version string of some type
     * @param version the version string
     * @return true if we're valid of some type, false otherwise
     */
    public static boolean isVersionString(String version){
        return toHeaderVersion(version) != null;
    }

    /**
     * are we a valid format string for some type
     * @param format the format string
     * @return true if we're valid of some type, false otherwise
     */
    public static boolean isFormatString(String format){
        for (VCFHeaderVersion hv : VCFHeaderVersion.values())
            if (hv.formatString.equals(format))
                    return true;
        return false;
    }


    public String getVersionString() {
        return versionString;
    }

    public String getFormatString() {
        return formatString;
    }
}
