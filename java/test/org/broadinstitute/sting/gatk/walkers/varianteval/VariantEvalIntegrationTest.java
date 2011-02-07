package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = validationDataLocation + "/VariantEval";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf";

    private static String rootGZ = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf.gz" +
            " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf.gz";

    private static String[] testsEnumerations = {root, rootGZ};

    private String cmdLineBuilder(String ... arguments) {
        String cmdline = "";

        for ( int argIndex = 0; argIndex < arguments.length; argIndex++ ) {
            cmdline += arguments[argIndex];

            if (argIndex < arguments.length - 1) {
                cmdline += " ";
            }
        }

        return cmdline;
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() {
//        nProcessedLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | awk '{ print length($4) }' | ~kiran/bin/SimpleStats = 38
//        nCalledLoci    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -c PASS = 9
//        nRefLoci       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 == ".") print $0 }' | wc -l = 4
//        nVariantLoci   = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != ".") print $0 }' | wc -l = 5
//        variantRate    = nVariantLoci / nProcessedLoci = 0.131578947
//        nSNPs          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($4) == 1 && length($5) == 1) print $0 }' | wc -l = 3
//        nInsertions    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($5) > 1) print $0 }' | wc -l = 1
//        nDeletions     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($4) > 1) print $0 }' | wc -l = 1
//        nNoCalls       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 4
//        nHets          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/1" || $i ~ "1/0") print $0 }' | wc -l = 8
//        nHomRef        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/0") print $0 }' | wc -l = 10
//        nHomVar        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "1/1") print $0 }' | wc -l = 5

        WalkerTestSpec spec = new WalkerTestSpec(
                                cmdLineBuilder(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-D " + b37dbSNP129,
                                        "-B:eval,VCF " + fundamentalTestVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-BTI eval",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("476b495de54e1a377c6895c02a6fdf6a")
                              );
        executeTest("testFundamentalsCountVariantsSNPsandIndels", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() {
//        nProcessedLociKnown = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | awk '{ print length($4) }' | ~kiran/bin/SimpleStats = 38
//        nCalledLociKnown    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }'= 3
//        nVariantLociKnown   = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }' | wc -l = 3
//        variantRateKnown    = nVariantLoci / nProcessedLoci = 0.0789473684
//        nSNPsKnown          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != "." && length($5) == 1) print $0 }' | wc -l = 3
//        nInsertionsKnown    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != "." && length($5) > 1) print $0 }' | wc -l = 0
//        nDeletionsKnown     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != "." && length($4) > 1) print $0 }' | wc -l = 0
//        nNoCallsKnown       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }' | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 0
//        nHetsKnown          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }' | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/1" || $i ~ "1/0") print $0 }' | wc -l = 3
//        nHomRefKnown        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }' | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/0") print $0 }' | wc -l = 1
//        nHomVarKnown        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | grep DBSNP129 | awk '{ if ($5 != ".") print $0 }' | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "1/1") print $0 }' | wc -l = 5

        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("9f4e4fff339e725f42d65063e43e7d1c")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNovelty", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() {
//        nProcessedLociFail  = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | awk '{ print length($4) }' | ~kiran/bin/SimpleStats = 38
//        nCalledLociFail     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -vc PASS = 3
//        nRefLociFail        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ if ($5 == ".") print $0 }' | wc -l = 1
//        nVariantLociFail    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ if ($5 != ".") print $0 }' | wc -l = 2
//        nSNPsFail           = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ if ($5 != "." && length($4) == 1 && length($5) == 1) print $0 }' | wc -l = 1
//        nInsertionsFail     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ if ($5 != "." && length($5) > 1) print $0 }' | wc -l = 0
//        nDeletionsFail      = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ if ($5 != "." && length($4) > 1) print $0 }' | wc -l = 1
//        nNoCallsFail        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 3
//        nHetsFail           = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/1" || $i ~ "1/0") print $0 }' | wc -l = 1
//        nHomRefFail         = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/0") print $0 }' | wc -l = 2
//        nHomVarFail         = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -v PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "1/1") print $0 }' | wc -l = 3

        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-ST Filter",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("369fa4f37bcc03b8a0bc1e58bf22bf0a")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNoveltyAndFilter", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() {
//        nProcessedLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | awk '{ print length($4) }' | ~kiran/bin/SimpleStats = 38
//        nCalledLoci    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep -c PASS = 8
//        nRefLoci       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ if ($5 == ".") print $0 }' | wc -l = 3
//        nVariantLoci   = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ if ($5 != ".") print $0 }' | wc -l = 5
//        nSNPs          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ if ($5 != "." && length($4) == 1 && length($5) == 1) print $0 }' | wc -l = 3
//        nInsertions    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ if ($5 != "." && length($5) > 1) print $0 }' | wc -l = 1
//        nDeletions     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ if ($5 != "." && length($4) > 1) print $0 }' | wc -l = 1
//        nNoCalls       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 4
//        nHets          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/1" || $i ~ "1/0") print $0 }' | wc -l = 8
//        nHomRef        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/0") print $0 }' | wc -l = 10
//        nHomVar        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep NONCPG | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "1/1") print $0 }' | wc -l = 5

        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST CpG",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("891ad0d38f1a1b08b31fe1cb6a3afc04")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithCpG", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() {
        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST FunctionalClass",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("d588179e2d9ed6e92a6ae1a80ac04270")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithFunctionalClass", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() {
        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Degeneracy",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("ceb0f5d9e0ea99eb8d00bce2f7bc1b73")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithDegeneracy", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() {
//        HG00513
//          nSNPs = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if (length($4) == 1 && length($5) == 1) print $0 }' | awk '{ print $10 }' | grep -v '0/0' | grep -v '\.\/\.' | wc -l = 3
//          nInsertions = $ grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if (length($4) == 1 && length($5) > 1) print $0 }' | awk '{ print $10 }' | grep -v '0/0' | grep -v '\.\/\.' | wc -l = 1
//          nDeletions = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if (length($4) > 1 && length($5) == 1) print $0 }' | awk '{ print $10 }' | grep -v '0/0' | grep -v '\.\/\.' | wc -l = 0
//          nNoCalls = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($10 ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 2
//          nHets = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($10 ~ "0/1") print $0 }' | wc -l = 2
//          nHomRef = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($10 ~ "0/0") print $0 }' | wc -l = 3
//          nHomVar = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($10 ~ "1/1") print $0 }' | wc -l = 2

        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Sample",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("94ce29b34b9e2e4304fc1bbf3f971a7d")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithSample", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() {
//        nSNPs = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $0 }' | wc -l = 7
//        nRefLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && $6 == ".") print $0 }' | wc -l = 4
//        nVariantLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && $6 != ".") print $0 }' | wc -l = 3
//        nSNPs = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && length($5) == 1 && length($6) == 1 && $6 != ".") print $0 }' | wc -l = 3
//        nHomRef = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '0/0' = 9
//        nHets = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '0/1' = 3
//        nHomVar = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '1/1' = 5
//        nNoCalls = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '\.\/\.' = 4
        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthSelect",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("96de32970b204816ecd9a120b9d8782b")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithJexlExpression", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() {
//        nSNPs = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $0 }' | wc -l = 7
//        nRefLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && $6 == ".") print $0 }' | wc -l = 4
//        nVariantLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && $6 != ".") print $0 }' | wc -l = 3
//        nSNPs = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20 && length($5) == 1 && length($6) == 1 && $6 != ".") print $0 }' | wc -l = 3
//        nHomRef = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '0/0' = 9
//        nHets = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '0/1' = 3
//        nHomVar = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '1/1' = 5
//        nNoCalls = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk -F"[\t;]" '{ for (i = 1; i < NF; i++) if ($i ~ "DP=") print $i, $0 }' | sed 's/^DP=//' | awk '{ if ($1 < 20) print $11 "\n" $12 "\n" $13 }' | grep -c '\.\/\.' = 4
        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-D " + b37dbSNP129,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthLt20",
                        "-select 'DP > 20'",
                        "-selectName DepthGt20",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("aea882132eb6afdc93fbc70e8d6c50e2")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithMultipleJexlExpressions", spec);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() {
//        nProcessedLoci = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | awk '{ print length($4) }' | ~kiran/bin/SimpleStats = 38
//        nCalledLoci    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep -c PASS = 9
//        nRefLoci       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 == ".") print $0 }' | wc -l = 4
//        nVariantLoci   = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != ".") print $0 }' | wc -l = 5
//        variantRate    = nVariantLoci / nProcessedLoci = 0.131578947
//        nSNPs          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($4) == 1 && length($5) == 1) print $0 }' | wc -l = 3
//        nInsertions    = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($5) > 1) print $0 }' | wc -l = 1
//        nDeletions     = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ if ($5 != "." && length($4) > 1) print $0 }' | wc -l = 1
//        nNoCalls       = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "[[:punct:]]/[[:punct:]]") print $0 }' | wc -l = 4
//        nHets          = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/1" || $i ~ "1/0") print $0 }' | wc -l = 8
//        nHomRef        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "0/0") print $0 }' | wc -l = 10
//        nHomVar        = grep -v '#' FundamentalsTest.annotated.db.subset.snps_and_indels.vcf | grep PASS | awk '{ for (i = 10; i <= 12; i++) if ($i ~ "1/1") print $0 }' | wc -l = 5

        WalkerTestSpec spec = new WalkerTestSpec(
                cmdLineBuilder(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("2a6215c86918a2f5c6748dd74057e79e")
        );
        executeTest("testFundamentalsCountVariantsNoCompRod", spec);
    }

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -o %s",
                    1, Arrays.asList("4184a8d44f8c559c904e41edf464a467"));
            executeTestParallel("testSelect1", spec);
            //executeTest("testSelect1", spec);
        }
    }

//    @Test
//    public void testSelect2() {
//        String extraArgs = "-L 1:1-10,000,000";
//        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
//                1, Arrays.asList(""));
//        //executeTestParallel("testSelect2", spec);
//        executeTest("testSelect2", spec);
//    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFiles[] = {"GenotypeConcordanceEval.vcf", "GenotypeConcordanceEval.vcf.gz"};
        for (String vcfFile : vcfFiles) {
            WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -B:eval,VCF " + validationDataLocation + vcfFile + " -B:comp,VCF " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                    1,
                    Arrays.asList("1387fcf8d5c53ff2c820fe79cc999bcf"));
            executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
            //executeTest("testVEGenotypeConcordance" + vcfFile, spec);
        }

    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "b28516b4d3627d2eb017a5449284a4e4");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -mvq 0 -EV MendelianViolationEvaluator", "a49d897095586ddb72bfe9faf0291312");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();
            for (String tests : testsEnumerations) {
                WalkerTestSpec spec = new WalkerTestSpec( tests + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
                executeTestParallel("testVESimple", spec);
                //executeTest("testVESimple", spec);
            }
        }
    }

    @Test
    public void testVEComplex() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -mvq 30 -EV MendelianViolationEvaluator" +
                " -B:dbsnp_130,dbSNP " + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B:comp_hapmap,VCF " + validationDataLocation + "CEU_hapmap_nogt_23.vcf";


        expectations.put("", "65d74eb7eea2355989c389e8fa886c06");
        expectations.put(" -knownName comp_hapmap -knownName dbsnp", "776827cd5a6fcaf8b8508813e8dc023c");
        expectations.put(" -knownName comp_hapmap", "e99a89cbca2027c983edc00d31ea4ec9");
        for (String tests : testsEnumerations) {
            for (Map.Entry<String, String> entry : expectations.entrySet()) {
                String extraArgs2 = entry.getKey();
                String md5 = entry.getValue();

                WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs1 + extraArgs2 + " -o %s",
                        1, // just one output file
                        Arrays.asList(md5));
                executeTestParallel("testVEComplex", spec);
                //executeTest("testVEComplex", spec);
            }
        }
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -EV GenotypeConcordance -B:evalYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.ug.very.few.lines.vcf -B:compYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("49cb4a6126c5383abd9a49a6c22b8d93"));
        executeTestParallel("testCompVsEvalAC",spec);
        //executeTest("testCompVsEvalAC",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" -B:eval,vcf " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -tf " + testDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("b018ca01b374def71fcd1bf164326485"));
        executeTestParallel("testTranches",spec);
        //executeTest("testTranches",spec);
    }

    @Test
    public void testCompOverlap() {
        String extraArgs = "-T VariantEval -R " + b37KGReference + " -L " + validationDataLocation + "VariantEval/pacbio.hg19.intervals -B:comphapmap,vcf " + comparisonDataLocation + "Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf -B:eval,vcf " + validationDataLocation + "VariantEval/pacbio.ts.recalibrated.vcf -noEV -EV CompOverlap -sn NA12878 -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("e8b5561eb60ea98a9be4a45abee00e07"));
        executeTestParallel("testCompOverlap",spec);
        //executeTest("testCompOverlap",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() {
        String dbsnp = GATKDataLocation + "dbsnp_129_b37.rod";

        String extraArgs = "-T VariantEval -R " +
                           b37KGReference +
                           " -L 20" +
                           " -D " + dbsnp +
                           " -B:evalBI,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                           " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("9323a6ad62dedbdb08752411960db60f"));
        executeTestParallel("testEvalTrackWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() {
        String dbsnp = GATKDataLocation + "dbsnp_129_b37.rod";

        String extraArgs = "-T VariantEval -R " + b37KGReference +
                " -L 20" +
                " -D " + dbsnp +
                " -B:evalBI,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                " -B:evalBC,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bc.sites.vcf" +
                " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("ef23195331affd332af1de9d261fdd0a"));
        executeTestParallel("testMultipleEvalTracksWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleCompTracks() {
        String dbsnp = GATKDataLocation + "dbsnp_132_b37.vcf";

        String extraArgs =  "-T VariantEval" +
                           " -R " + b37KGReference +
                           " -B:comp,VCF " + validationDataLocation + "/VariantEval/ALL.phase1.chr20.broad.snps.genotypes.subset.vcf" +
                           " -B:eval,VCF " + validationDataLocation + "/VariantEval/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf" +
                           " -B:dbsnp,VCF " + dbsnp +
                           " -L 20:10000000-10100000" +
                           " -noST -noEV -ST Novelty -EV CompOverlap" +
                           " -o %s";

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("3fced8e5fa7a1c952d08fead0accd3fb"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

//    @Test
//    public void testVEGenomicallyAnnotated() {
//        String vecmd = "-T VariantEval" +
//                       " -R " + b36KGReference +
//                       " -L 21" +
//                       " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
//                       " -EV CountFunctionalClasses -noEV" +
//                       " -B:eval,VCF " + validationDataLocation + "test.filtered.maf_annotated.vcf" +
//                       " -o %s";
//        String md5 = "";
//
//        WalkerTestSpec spec = new WalkerTestSpec(vecmd, 1, Arrays.asList(md5));
//        executeTestParallel("testVEGenomicallyAnnotated", spec);
//        //executeTest("testVEGenomicallyAnnotated", spec);
//    }
//
//    @Test
//    public void testVEWriteVCF() {
//        String extraArgs = "-L 1:1-10,000,000 -NO_HEADER -family NA19238+NA19239=NA19240 -mvq 30 -EV MendelianViolationEvaluator";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs + " -o %s -outputVCF %s -NO_HEADER",
//                    2,
//                    Arrays.asList("50321436a65ef7d574286cb0a1c55f7e", "d4bdd06ed5cb1aff1dfee8b69d5d17b8"));
//            executeTestParallel("testVEWriteVCF", spec);
//            //executeTest("testVEWriteVCF", spec);
//        }
//    }
//
//    @Test
//    public void testVEValidatePass() {
//        String extraArgs = "-L 1:1-10,000,000";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(withValidateTiTv(withSelect(tests, "DP < 50", "DP50"), 1.0, 4.0) + " " + extraArgs + " -o %s",
//                    1, Arrays.asList("8a0203f0533b628ad7f1f230a43f105f"));
//            executeTestParallel("testVEValidatePass", spec);
//        }
//    }
//
//    @Test
//    public void testVEValidateFail() {
//        String extraArgs = "-L 1:1-10,000,000";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(withValidateTiTv(withSelect(tests, "DP < 50", "DP50"), 1.0, 1.2) + " " + extraArgs + " -o %s",
//                    1, UserException.class);
//            executeTestParallel("testVEValidateFail", spec);
//        }
//    }
//
//    private static String withValidateTiTv(String cmd, double min, double max) {
//        return String.format("%s -validate 'eval.comp_genotypes.all.called.all.titv.tiTvRatio >= %2$s' -validate 'eval.comp_genotypes.all.called.all.titv.tiTvRatio <= %3$s'", cmd, min, max);
//    }

}
