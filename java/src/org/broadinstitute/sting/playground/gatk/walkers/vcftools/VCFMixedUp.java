package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;

import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.ArrayList;

public class VCFMixedUp extends RefWalker<Integer, Integer> {
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for (ReferenceOrderedDatum rod : tracker.getAllRods()) {
            if (rod instanceof RodVCF) {
                for (String rodfile : this.getToolkit().getArguments().RODBindings) {
                    String[] rodpieces = rodfile.split(",");

                    if (rodpieces[2].contains(".vcf")) {
                        try {
                            FileReader reader = new FileReader(rodpieces[2]);
                            BufferedReader breader = new BufferedReader(reader);

                            ArrayList<String> header = new ArrayList<String>();
                            ArrayList<String> variant = new ArrayList<String>();

                            String line;
                            while ((line = breader.readLine()) != null) {
                                String[] pieces = line.split("\\s+");

                                if (line.contains("##")) {
                                } else if (line.contains("#CHROM")) {
                                    for (int i = 0; i < pieces.length; i++) {
                                        header.add(i, pieces[i]);
                                    }
                                } else {
                                    for (int i = 0; i < pieces.length; i++) {
                                        variant.add(i, pieces[i]);
                                    }

                                    for (int i = 0; i < pieces.length; i++) {
                                        if (i < 9 || header.get(i).contains("NA12892") || header.get(i).contains("NA10851")) {
                                            out.printf("%s => %s\n", header.get(i), variant.get(i));
                                        }
                                    }
                                }
                            }
                        } catch (IOException e) {}
                    }
                }

                RodVCF vcfrod = (RodVCF) rod;
                VCFRecord record = vcfrod.mCurrentRecord;

                out.println();

                for (VCFGenotypeRecord gr : record.getVCFGenotypeRecords()) {
                    if (gr.getSampleName().equalsIgnoreCase("NA12892") || gr.getSampleName().equalsIgnoreCase("NA10851")) {
                        out.printf("%s: %s:%s:%s:%s\n", gr.getSampleName(),
                                                        gr.toGenotypeString(record.getAlternateAlleles()),
                                                        gr.getFields().get("GQ"),
                                                        gr.getFields().get("DP"),
                                                        gr.getFields().get("HQ"));
                    }
                }
            }
        }

        return 0;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }
}
