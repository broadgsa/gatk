## Phred-scaled Quality Scores

http://gatkforums.broadinstitute.org/gatk/discussion/4260/phred-scaled-quality-scores

<p>You may have noticed that a lot of the scores that are output by the GATK are in Phred scale. The Phred scale was originally used to represent base quality scores emitted by the Phred program in the early days of the Human Genome Project (see <a href="http://en.wikipedia.org/wiki/Phred_quality_score">this Wikipedia article</a> for more historical background). Now they are widely used to represent probabilities and confidence scores in other contexts of genome science.</p>
<h3>Phred scale in context</h3>
<p>In the context of sequencing, Phred-scaled quality scores are used to represent how confident we are in the assignment of each base call by the sequencer. </p>
<p>In the context of variant calling, Phred-scaled quality scores can be used to represent many types of probabilities. The most commonly used in GATK is the QUAL score, or variant quality score. It is used in much the same way as the base quality score: the variant quality score is a Phred-scaled estimate of how confident we are that the variant caller correctly identified that a given genome position displays variation in at least one sample. </p>
<h3>Phred scale in practice</h3>
<p>In today’s sequencing output, by convention, most useable Phred-scaled base quality scores range from 2 to 40, with some variations in the range depending on the origin of the sequence data (see the <a href="https://en.wikipedia.org/wiki/FASTQ_format#Encoding">FASTQ format</a> documentation for details). However, Phred-scaled quality scores in general can range anywhere from 0 to infinity. A <strong>higher score</strong> indicates a higher probability that a particular decision is <strong>correct</strong>, while conversely, a <strong>lower score</strong> indicates a higher probability that the decision is <strong>incorrect</strong>. </p>
<p>The Phred quality score (Q) is logarithmically related to the error probability (E).</p>
<p>$$ Q = -10 \log E $$</p>
<p>So we can interpret this score as an estimate of <strong>error</strong>, where the error is <em>e.g.</em> the probability that the base is called <strong>incorrectly</strong> by the sequencer, but we can also interpret it as an estimate of <strong>accuracy</strong>, where the accuracy is <em>e.g.</em> the probability that the base was identified <strong>correctly</strong> by the sequencer. Depending on how we decide to express it, we can make the following calculations:</p>
<p>If we want the probability of error (E), we take:</p>
<p>$$ E = 10 ^{-\left(\frac{Q}{10}\right)} $$ </p>
<p>And conversely, if we want to express this as the estimate of accuracy (A), we simply take </p>
<p>$$
\begin{eqnarray}
A &amp;=&amp; 1 - E  \nonumber \
&amp;=&amp; 1 - 10 ^{-\left(\frac{Q}{10}\right)}  \nonumber \
\end{eqnarray}
$$</p>
<p>Here is a table of how to interpret a range of Phred Quality Scores. It is largely adapted from the Wikipedia page for Phred Quality Score.</p>
<p>For many purposes, a Phred Score of 20 or above is acceptable, because this means that whatever it qualifies is 99% accurate, with a 1% chance of error. </p>
<table class="table table-striped">
<thead>
<tr>
<th>Phred Quality Score</th>
<th>Error</th>
<th>Accuracy (1 - Error)</th>
</tr>
</thead>
<tbody>
<tr>
<td>10</td>
<td>1/10 = 10%</td>
<td>90%</td>
</tr>
<tr>
<td>20</td>
<td>1/100 = 1%</td>
<td>99%</td>
</tr>
<tr>
<td>30</td>
<td>1/1000 = 0.1%</td>
<td>99.9%</td>
</tr>
<tr>
<td>40</td>
<td>1/10000 = 0.01%</td>
<td>99.99%</td>
</tr>
<tr>
<td>50</td>
<td>1/100000 = 0.001%</td>
<td>99.999%</td>
</tr>
<tr>
<td>60</td>
<td>1/1000000 = 0.0001%</td>
<td>99.9999%</td>
</tr>
</tbody>
</table>
<p>And finally, here is a graphical representation of the Phred scores showing their relationship to accuracy and error probabilities. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/78/663145e9df43db3efe5df4d0b88cf4.png" />
<p>The red line shows the error, and the blue line shows the accuracy. Of course, as error decreases, accuracy increases symmetrically.  </p>
<p>Note: You can see that below Q20 (which is how we usually refer to a Phred score of 20), the curve is really steep, meaning that as the Phred score decreases, you lose confidence very rapidly. In contrast, above Q20, both of the graphs level out. This is why Q20 is a good cutoff score for many basic purposes.</p>