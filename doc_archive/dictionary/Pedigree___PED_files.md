## Pedigree / PED files

http://gatkforums.broadinstitute.org/gatk/discussion/7696/pedigree-ped-files

<p>A pedigree is a structured description of the familial relationships between samples. </p>
<p>Some GATK tools are capable of incorporating pedigree information in the analysis they perform if provided in the form of a PED file through the <code>--pedigree</code> (or <code>-ped</code>) argument. </p>
<hr />
<h3>PED file format</h3>
<p>PED files are tabular text files describing meta-data about the samples. See <a href="http://www.broadinstitute.org/mpg/tagger/faq.html"><a href="http://www.broadinstitute.org/mpg/tagger/faq.html">http://www.broadinstitute.org/mpg/tagger/faq.html</a></a> and <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped"><a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped">http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped</a></a> for more information.</p>
<p>The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:</p>
<ul>
<li>Family ID</li>
<li>Individual ID</li>
<li>Paternal ID</li>
<li>Maternal ID</li>
<li>Sex (1=male; 2=female; other=unknown)</li>
<li>Phenotype</li>
</ul>
<p>The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person. If an individual's sex is unknown, then any character other than 1 or 2 can be used in the fifth column.</p>
<p>A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a quantitative trait or an &quot;affected status&quot; column: GATK will automatically detect which type (i.e. based on whether a value other than 0, 1, 2 or the missing genotype code is observed). </p>
<p>Affected status should be coded as follows:</p>
<ul>
<li>-9 missing</li>
<li>0 missing</li>
<li>1 unaffected</li>
<li>2 affected</li>
</ul>
<p>If any value outside of -9,0,1,2 is detected, then the samples are assumed to have phenotype values, interpreted as string phenotype values.</p>
<p>Note that genotypes (column 7 onwards) cannot be specified to the GATK.</p>
<p>You can add a comment to a PED or MAP file by starting the line with a # character. The rest of that line will be ignored, so make sure none of the IDs start with this character.</p>
<p>Each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to tell the GATK PED parser that the corresponding fields are missing from the ped file.</p>
<h4>Example</h4>
<p>Here are two individuals (one row = one person):</p>
<pre>
FAM001  1  0 0  1  2
FAM001  2  0 0  1  2
</pre>