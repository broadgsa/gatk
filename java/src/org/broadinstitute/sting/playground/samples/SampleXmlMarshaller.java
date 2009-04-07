package org.broadinstitute.sting.playground.samples;

import org.exolab.castor.xml.*;
import org.apache.log4j.BasicConfigurator;

import java.io.FileWriter;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 7, 2009
 * Time: 2:16:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class SampleXmlMarshaller {
    public enum RMDType { HapMap };

    private RMDType type;
    private String fileVersion;
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

    public static void main( String argv[] ) {
        if( argv.length == 0 || (!argv[0].equalsIgnoreCase("marshal") && !argv[0].equalsIgnoreCase("unmarshal")) ) {
            System.out.println( "USAGE: java org.broadinstitute.sting.playground.samples.SampleXmlMarshaller {marshal | unmarshal}");
            System.exit(-1);
        }

        BasicConfigurator.configure();

        if( argv[0].equalsIgnoreCase("marshal") )
            marshal();
        else if( argv[0].equalsIgnoreCase("unmarshal") )
            unmarshal();
    }

    public static void marshal() {
        SampleXmlMarshaller descriptor = new SampleXmlMarshaller( RMDType.HapMap, "foo.hapmap", "1.0" );
        System.out.printf( "Marshalling descriptor = %s%n", descriptor );
        try {
            FileWriter writer = new FileWriter( "foo.xml" );
            Marshaller.marshal( descriptor, writer );
        }
        catch( IOException ex ) {
            System.out.println("Unable to open writer:" + ex);
        }
        catch( MarshalException ex ) {
            System.out.println("Unable to marshal data:" + ex);
        }
        catch( ValidationException ex ) {
            System.out.println("Unable to validate data: " + ex);
        }
    }

    public static void unmarshal() {
        SampleXmlMarshaller descriptor = null;
        try {
            FileReader reader = new FileReader( "foo.xml" );
            descriptor = (SampleXmlMarshaller)Unmarshaller.unmarshal( SampleXmlMarshaller.class, reader );
        }
        catch( FileNotFoundException ex ) {
            System.out.println("Unable to open reader:" + ex);
        }
        catch( MarshalException ex ) {
            System.out.println("Unable to unmarshal data:" + ex);
        }
        catch( ValidationException ex ) {
            System.out.println("Unable to validate data: " + ex);
        }

        System.out.printf( "Unmarshalled descriptor = %s%n", descriptor );

    }
}
