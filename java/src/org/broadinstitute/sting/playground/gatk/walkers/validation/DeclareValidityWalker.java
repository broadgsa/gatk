/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.validation;

import com.sun.tools.corba.se.idl.constExpr.Not;
import com.sun.tools.internal.xjc.reader.xmlschema.bindinfo.BIConversion;
import org.broad.tribble.Feature;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.bed.BedParser;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import javax.activation.*;
import java.beans.VetoableChangeSupport;
import java.io.PrintStream;
import java.io.Writer;
import java.util.*;


/**
 * Declares the validity of variants in a vcf as either true or false. For use with the IGV crowd-sourcing bed generation
 */

@Requires(value={},referenceMetaData=@RMD(name="validated", type=VariantContext.class))
public class DeclareValidityWalker   extends RodWalker<Integer, Integer>{


    @Output(doc = "Write bed to this file instead of standard out")
    public PrintStream out;

    @Argument(fullName = "validity", shortName = "V",
            doc = "Rank of variant validity on a 0-4 scale where 0 is definitely false positive; 4 is definitely true positive.")
    int validity;

    @Argument(fullName = "Note", shortName = "N", doc = "Annotation to be included in FP/TP field", required = false)
    String Note =".";

    @Argument(fullName = "Source", shortName = "s", doc = "Institutional source of annotation", required = false)
    String Source = ".";

    @Argument(fullName = "Build", shortName = "bld", doc = "Genome build", required = false)
    String build = ".";


    @Override
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }



    public class ValidityDeclaration{
        protected String contig;
        protected long start;
        protected long stop;
        protected int valid;
        protected Allele refBase;
        protected Allele altBase;
        protected String Note;
        protected String Source;
        protected String Build;
        protected String user;

        public String getBuild(){
        String refPath = getToolkit().getArguments().referenceFile.getPath();
        if (refPath.contains("19")) {return "hg19";}
        else if (refPath.contains("18")) {return "hg18";}
        else if (refPath.contains("36")) {return "b36";}
        else if (refPath.contains("37")) {return "b37";}
        else {return "unknown";}
    }

        public ValidityDeclaration(GenomeLoc Loc, VariantContext Con, int validity, String Note, String Source, String Build){   //Constructor expects  1 based
            this.contig=Loc.getContig();
            this.start=Loc.getStart()-1;
            this.stop=Loc.getStop();
            this.valid =validity;
            this.altBase = Con.getAlternateAllele(0);
            if (Con.getAlternateAlleles().toArray().length >1)
            {
                logger.warn("***NOTE: Only the first alternate allele in a VCF will be declared as " +valid+"***");
            }
            this.refBase = Con.getReference();
            this.Note = Note;
            this.Source = Source;
            this.user = System.getenv("USER");
            if (Build == "."){
                this.Build=getBuild();
            }
            else{this.Build = Build; }

        }
        public String toString() {
            return String.format("%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s", contig, start, stop, refBase,altBase, validity, user, Build, Note, Source);
        }
    }
 /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
     public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context){
         if ( tracker == null )
            return 0;


         VariantContext current = tracker.getVariantContext(ref, "validated", context.getLocation());
         if (current == null) {
             return 0;}


         ValidityDeclaration bedLine = new ValidityDeclaration(ref.getLocus(), current, validity, Note, Source, build);
         out.println(bedLine);
         return 1;
     }

    public Integer reduce(Integer counter, Integer sum) {
        return counter+sum;
    }

}
