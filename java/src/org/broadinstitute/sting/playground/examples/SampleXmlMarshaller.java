package org.broadinstitute.sting.playground.examples;

import org.apache.log4j.BasicConfigurator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.simpleframework.xml.Element;
import org.simpleframework.xml.Root;
import org.simpleframework.xml.Serializer;
import org.simpleframework.xml.stream.Format;
import org.simpleframework.xml.stream.HyphenStyle;
import org.simpleframework.xml.core.Persister;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 7, 2009
 * Time: 2:16:43 PM
 * To change this template use File | Settings | File Templates.
 */

@Root // this is a root output object for simpleXML
public class SampleXmlMarshaller {
    public enum RMDType { HapMap };

    @Element // simpleXML tag to say make it an element, you could also specifiy Attribute
    private RMDType type;

    @Element
    private String fileVersion;

    @Element
    private String fileName;

    public SampleXmlMarshaller() {
    }

    public SampleXmlMarshaller( RMDType type, String fileName, String fileVersion ) {
        setType( type );
        setFileName( fileName );
        setFileVersion( fileVersion );
    }

    public RMDType getType() {
        return type;
    }

    public void setType( RMDType type ) {
        this.type = type;
    }

    public String getFileVersion() {
        return fileVersion;
    }

    public void setFileVersion( String fileVersion ) {
        this.fileVersion = fileVersion;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName( String fileName ) {
        this.fileName = fileName;
    }

    public String toString() {
        return String.format("Type = %s, Name = %s, Version = %s%n", type, fileName, fileVersion );
    }

    public boolean equals(SampleXmlMarshaller other) {
        if (other == null) { return false; }
        if (other.getType() == this.getType() &&
                other.getFileVersion().equals(this.getFileVersion()) &&
                other.getFileName().equals(this.getFileName())) {
            return true;
        }
        return false;
    }

    public static void main( String argv[] ) {
        if (argv.length != 1) {
            System.err.println("You must specify a filename for the XML output.");
            System.err.println("\nUsage: SampleXmlMarshaller outputfile\n");       
        }
        // our SampleXmlMarshallers
        SampleXmlMarshaller startWith = new SampleXmlMarshaller(RMDType.HapMap, "testFile.xml", "0.1");
        SampleXmlMarshaller endWith = null;

        // where to read and write to
        String writeTo = argv[0];

        BasicConfigurator.configure();

        marshal(startWith, writeTo);
        endWith = unmarshal(writeTo);

        if (startWith.equals(endWith)) {
            System.out.println("they're equal, check " + writeTo.toString() + " for the output");
        } else {
            System.out.println("they're NOT equal, check " + writeTo.toString() + " for the output.  Something must of gone wrong");
        }

    }

    public static void marshal(SampleXmlMarshaller sample, String filename) {
       Serializer serializer = new Persister(new Format(new HyphenStyle()));
        File out = new File(filename);
        try {
            serializer.write(sample, out);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public static SampleXmlMarshaller unmarshal(String filename) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        File source = new File(filename);
        try {
            SampleXmlMarshaller example = serializer.read(SampleXmlMarshaller.class, source);
            return example;
        } catch (Exception e) {
            throw new ReviewedStingException("Failed to marshal the data to file " + filename,e);
        }
    }
}
