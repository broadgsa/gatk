## OxoG oxidative artifacts

http://gatkforums.broadinstitute.org/gatk/discussion/6328/oxog-oxidative-artifacts

<p>Oxidation of guanine to 8-oxoguanine is one of the most common <strong>pre-adapter artifacts</strong> associated with genomic library preparation, arising from a combination of heat, shearing, and metal contaminates in a sample (doi: 10.1093/nar/gks1443). The 8-oxoguanine base can pair with either cytosine or adenine, ultimately leading to G→T transversion mutations during PCR amplification. </p>
<p>This occurs when a G on the template strand is oxidized, giving it an affinity for binding to A rather than the usual C. Thus, PCR will introduce apparent G&gt;T substitutions in read 1 and C&gt;A in read 2. In the resulting alignments, a given G&gt;T or C&gt;A observation could either be: </p>
<ol>
<li>a true mutation </li>
<li>an 8-oxoguanine artifact</li>
<li>some other kind of artifact.  </li>
</ol>
<p>The variants (C→A)/(G→T) tend to occur in specific sequence contexts e.g. CCG→CAG (doi:10.1093/nar/gks1443).  Although occurring at relatively low frequencies, these artifacts can have profound impacts on variant calling fidelity (doi:10.1093/nar/gks1443).</p>