package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;

import java.io.*;
import java.util.*;

public class VcfToGeliText extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Converts VCF files to simple Geli text files.";
    @Option(shortName="I", doc="Input file (vcf) to convert.", optional=false) public File IN = null;
    @Option(shortName="O",doc="Output file (gelitext). If not specified, output is printed to stdout.", optional=true) public File OUT = null;
    @Option(shortName="sample",doc="Sample number to extract from the VCF. If not specified, it takes the firt one.", optional=true) public Integer sample = 1;

    public static void main(final String[] argv) {
        System.exit(new VcfToGeliText().instanceMain(argv));
    }

    protected int doWork() {

        if ( IN == null )
            throw new RuntimeException("No input VCF file provided!");

        FileReader in;
        BufferedReader vcf;
        try {
            in = new FileReader(IN);
            vcf = new BufferedReader(in);
        } catch ( FileNotFoundException ie ) {
            System.out.println("Failed to open input file "+IN+": "+ie.getCause());
            return 1;
        }

        PrintStream out;
        if ( OUT == null ) out = System.out;
        else {
            try {
                out =  new PrintStream(OUT);
            } catch ( FileNotFoundException ie ) {
                System.out.println("Failed to open output file "+OUT+": "+ie.getCause());
                return 1;
            }
        }

        String currentline;
        try {
            while ( (currentline = vcf.readLine()) != null ) {
                if ( currentline.length() == 0 || currentline.charAt(0) == '#' )
                    continue;

                StringTokenizer st = new StringTokenizer(currentline);
                String chr = st.nextToken();
                String pos = st.nextToken();
                st.nextToken();
                String ref = st.nextToken();
                String altStr = st.nextToken();
                for (int i = 0; i < 4; i++)
                    st.nextToken();
                for (int i = 1; i < sample; i++)
                    st.nextToken();
                String sampleStr = st.nextToken();

                HashMap<String,String> alleles = new HashMap<String,String>();
                alleles.put("0", ref);
                StringTokenizer stAlt = new StringTokenizer(altStr, ",");
                int index = 1;
                while ( stAlt.hasMoreTokens() )
                    alleles.put(String.valueOf(index++), stAlt.nextToken());

                StringTokenizer stSample = new StringTokenizer(sampleStr, "/:");
                String genotype1 = stSample.nextToken();
                String genotype2 = stSample.nextToken();
                if ( genotype1.equals("0") && genotype2.equals("0") )
                    continue;

                out.println(chr + "\t" + pos + "\t" + ref + "\t0\t0\t"
                        + alleles.get(genotype1) + alleles.get(genotype2)
                        + "\t30\t30\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");
            }
            vcf.close();
            in.close();
            out.close();
        } catch (IOException e) {
            System.out.println("Jaa I/O exception"+e.getCause());
            return 1;
        }

        return 0;
    }


}
