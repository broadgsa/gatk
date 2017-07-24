## Collected FAQs about VCF files

http://gatkforums.broadinstitute.org/gatk/discussion/1318/collected-faqs-about-vcf-files

<h3>1. What file formats do you support for variant callsets?</h3>
<p>We support the <a href="http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0">Variant Call Format (VCF)</a>  for variant callsets.  No other file formats are supported.</p>
<h3>2. How can I know if my VCF file is valid?</h3>
<p><a href="http://vcftools.sourceforge.net/">VCFTools</a> contains a <a href="http://vcftools.sourceforge.net/docs.html#validator">validation tool</a> that will allow you to verify it.</p>
<h3>3. Are you planning to include any converters from different formats or allow different input formats than VCF?</h3>
<p>No, we like VCF and we think it's important to have a good standard format. Multiplying formats just makes life hard for everyone, both developers and analysts. </p>