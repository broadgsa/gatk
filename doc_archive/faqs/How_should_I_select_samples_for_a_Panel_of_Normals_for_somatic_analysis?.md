## How should I select samples for a Panel of Normals for somatic analysis?

http://gatkforums.broadinstitute.org/gatk/discussion/7366/how-should-i-select-samples-for-a-panel-of-normals-for-somatic-analysis

<p>The Panel of Normals (PoN) plays two important roles in somatic variant analysis: </p>
<ol>
<li>Exclude germline variant sites that are found in the normals to avoid calling them as potential somatic variants in the tumor;</li>
<li>Exclude technical artifacts that arise from particular techniques (eg sample preservation) and technologies (eg library capture, sequencing chemistry).</li>
</ol>
<p>Given these roles, the most important selection criteria are the technical properties of how the normal data was generated. It's very important to use normals that are as technically similar as possible to the tumor. Also, the samples should come from subjects that were young and healthy (to minimize the chance of using as normal a sample from someone who has an undiagnosed tumor).</p>
<p>If possible it is better to use normals generated from the same type of tissue because if the tissues were preserved differently, the artifact patterns may be different. </p>