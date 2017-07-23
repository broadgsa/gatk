## What is Map/Reduce and why are GATK tools called "walkers"?

http://gatkforums.broadinstitute.org/gatk/discussion/1754/what-is-map-reduce-and-why-are-gatk-tools-called-walkers

<h3>Overview</h3>
<p>One of the key challenges of working with next-gen sequence data is that input files are usually very large. We can’t just make the program open the files, load all the data into memory and perform whatever analysis is needed on all of it in one go. It’s just too much work, even for supercomputers.</p>
<p>Instead, we make the program cut the job into smaller tasks that the computer can easily process separately. Then we have it combine the results of each step into the final result.</p>
<h3>Map/Reduce</h3>
<p><strong>Map/Reduce</strong> is the technique we use to achieve this. It consists of three steps formally called <code>filter</code>, <code>map</code> and <code>reduce</code>. Let’s apply it to an example case where we want to find out what is the average depth of coverage in our dataset for a certain region of the genome.</p>
<ul>
<li>
<p><code>filter</code> determines what subset of the data needs to be processed in each task. In our example, the program lists all the reference positions in our region of interest.</p>
</li>
<li>
<p><code>map</code> applies the function, <em>i.e.</em> performs the analysis on each subset of data. In our example, for each position in the list, the program looks into the BAM file, pulls out the pileup of bases and outputs the depth of coverage at that position.</p>
</li>
<li><code>reduce</code> combines the elements in the list of results output by the <code>map</code> function. In our example, the program takes the coverage numbers that were calculated separately for all the reference positions and calculates their average, which is the final result we want.</li>
</ul>
<p>This may seem trivial for such a simple example, but it is a very powerful method with many advantages. Among other things, it makes it relatively easy to parallelize operations, which makes the tools run much faster on large datasets.</p>
<h3>Walkers, filters and traversal types</h3>
<p>All the tools in the GATK are built from the ground up to take advantage of this method. That’s why we call them <strong>walkers</strong>: because they “walk” across the genome, getting things done.</p>
<p>Note that even though it’s not included in the Map/Reduce technique’s name, the <code>filter</code> step is very important. It determines what data get presented to the tool for analysis, selecting only the appropriate data for each task and discarding anything that’s not relevant. This is a key part of the Map/Reduce technique, because that’s what makes each task “bite-sized” enough for the computer to handle easily.</p>
<p>Each tool has filters that are tailored specifically for the type of analysis it performs. The filters rely on <strong>traversal engines</strong>, which are little programs that are designed to “traverse” the data (<em>i.e.</em> walk through the data) in specific ways.</p>
<p>There are three major types of traversal: <strong>Locus Traversal</strong>, <strong>Read Traversal</strong> and <strong>Active Region Traversal</strong>.  In our interval coverage example, the tool’s filter uses the <strong>Locus Traversal</strong> engine, which walks through the data by locus, <em>i.e.</em> by position along the reference genome. Because of that, the tool is classified as a <strong>Locus Walker</strong>.  Similarly, the <strong>Read Traversal</strong> engine is used, you’ve guessed it, by <strong>Read Walkers</strong>. </p>
<p>The GATK engine comes packed with many other ways to walk through the genome and get the job done seamlessly, but those are the ones you’ll encounter most often. </p>
<h3>Further reading</h3>
<p><a href="http://www.broadinstitute.org/gatk/guide/article?id=1988">A primer on parallelism with the GATK</a>
<a href="http://www.broadinstitute.org/gatk/guide/article?id=1975">How can I use parallelism to make GATK tools run faster?</a></p>