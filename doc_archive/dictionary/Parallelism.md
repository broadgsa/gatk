## Parallelism

http://gatkforums.broadinstitute.org/gatk/discussion/1988/parallelism

<p><em>This document explains the concepts involved and how they are applied within the GATK (and Crom+WDL or Queue where applicable). For specific configuration recommendations, see the companion document on <a href="http://www.broadinstitute.org/gatk/guide/article?id=1975">parallelizing GATK tools</a>.</em></p>
<hr />
<h2>1. The concept of parallelism</h2>
<p>Parallelism is a way to make a program finish faster by performing several operations in parallel, rather than sequentially (<em>i.e.</em> waiting for each operation to finish before starting the next one).</p>
<p>Imagine you need to cook rice for sixty-four people, but your rice cooker can only make enough rice for four people at a time. If you have to cook all the batches of rice sequentially, it's going to take all night. But if you have eight rice cookers that you can use in parallel, you can finish up to eight times faster.</p>
<p>This is a very simple idea but it has a key requirement: you have to be able to break down the job into smaller tasks that can be done independently. It's easy enough to divide portions of rice because rice itself is a collection of discrete units. In contrast, let's look at a case where you can't make that kind of division: it takes one pregnant woman nine months to grow a baby, but you can't do it in one month by having nine women share the work. </p>
<p>The good news is that most GATK runs are more like rice than like babies. Because GATK tools are built to use the Map/Reduce method (see <a href="http://www.broadinstitute.org/gatk/guide/article?id=1754">doc</a> for details), most GATK runs essentially consist of a series of many small independent operations that can be parallelized.</p>
<h3>A quick warning about tradeoffs</h3>
<p>Parallelism is a great way to speed up processing on large amounts of data, but it has &quot;overhead&quot; costs. Without getting too technical at this point, let's just say that parallelized jobs need to be managed, you have to set aside memory for them, regulate file access, collect results and so on. So it's important to balance the costs against the benefits, and avoid dividing the overall work into too many small jobs.</p>
<p>Going back to the introductory example, you wouldn't want to use a million tiny rice cookers that each boil a single grain of rice. They would take way too much space on your countertop, and the time it would take to distribute each grain then collect it when it's cooked would negate any benefits from parallelizing in the first place.</p>
<h3>Parallel computing in practice (sort of)</h3>
<p>OK, parallelism sounds great (despite the tradeoffs caveat), but how do we get from cooking rice to executing programs? What actually happens in the computer?</p>
<p>Consider that when you run a program like the GATK, you're just telling the computer to execute a set of instructions.</p>
<p>Let's say we have a text file and we want to count the number of lines in it. The set of instructions to do this can be as simple as:</p>
<ul>
<li><code>open the file, count the number of lines in the file, tell us the number, close the file</code></li>
</ul>
<p><em>Note that <code>tell us the number</code> can mean writing it to the console, or storing it somewhere for use later on.</em></p>
<p>Now let's say we want to know the number of words on each line. The set of instructions would be:</p>
<ul>
<li><code>open the file, read the first line, count the number of words, tell us the number, read the second line, count the number of words, tell us the number, read the third line, count the number of words, tell us the number</code></li>
</ul>
<p>And so on until we've read all the lines, and finally we can close the file. It's pretty straightforward, but if our file has a lot of lines, it will take a long time, and it will probably not use all the computing power we have available.</p>
<p>So to parallelize this program and save time, we just cut up this set of instructions into separate subsets like this:</p>
<ul>
<li>
<p><code>open the file, index the lines</code>  </p>
</li>
<li><code>read the first line, count the number of words, tell us the number</code></li>
<li><code>read the second line, count the number of words, tell us the number</code></li>
<li><code>read the third line, count the number of words, tell us the number</code></li>
<li>
<p><code>[repeat for all lines]</code></p>
</li>
<li><code>collect final results and close the file</code></li>
</ul>
<p>Here, the <code>read the Nth line</code> steps can be performed in parallel, because they are all independent operations.</p>
<p>You'll notice that we added a step, <code>index the lines</code>. That's a little bit of peliminary work that allows us to perform the <code>read the Nth line</code> steps in parallel (or in any order we want) because it tells us how many lines there are and where to find each one within the file. It makes the whole process much more efficient. As you may know, the GATK requires index files for the main data files (reference, BAMs and VCFs); the reason is essentially to have that indexing step already done.</p>
<p>Anyway, that's the general principle: you transform your linear set of instructions into several subsets of instructions. There's usually one subset that has to be run first and one that has to be run last, but all the subsets in the middle can be run at the same time (in parallel) or in whatever order you want.</p>
<hr />
<h2>2. Parallelizing the GATK</h2>
<p>There are three different modes of parallelism offered by the GATK, and to really understand the difference you first need to understand what are the different <em>levels of computing</em> that are involved.</p>
<h3>A quick word about levels of computing</h3>
<p>By <em>levels of computing</em>, we mean the computing units in terms of hardware: the core, the machine (or CPU) and the cluster or cloud.</p>
<ul>
<li>
<p><strong>Core:</strong> the level below the machine. On your laptop or desktop, the CPU (central processing unit, or processor) contains one or more cores. If you have a recent machine, your CPU probably has at least two cores, and is therefore called dual-core. If it has four, it's a quad-core, and so on. High-end consumer machines like the latest Mac Pro have up to twelve-core CPUs (which should be called dodeca-core if we follow the Latin terminology) but the CPUs on some professional-grade machines can have tens or hundreds of cores.</p>
</li>
<li>
<p><strong>Machine:</strong> the middle of the scale. For most of us, the machine is the laptop or desktop computer.  Really we should refer to the CPU specifically, since that's the relevant part that does the processing, but the most common usage is to say <strong>machine</strong>. Except if the machine is part of a cluster, in which case it's called a <strong>node</strong>.</p>
</li>
<li><strong>Cluster or cloud:</strong> the level above the machine. This is a high-performance computing structure made of a bunch of machines (usually called <strong>nodes</strong>) networked together. If you have access to a cluster, chances are it either belongs to your institution, or your company is renting time on it. A cluster can also be called a <strong>server farm</strong> or a <strong>load-sharing facility</strong>.</li>
</ul>
<p>Parallelism can be applied at all three of these levels, but in different ways of course, and under different names. Parallelism takes the name of <strong>multi-threading</strong> at the core and machine levels, and <strong>scatter-gather</strong> at the cluster level.</p>
<h3>Multi-threading</h3>
<p>In computing, a <strong>thread of execution</strong> is a set of instructions that the program issues to the processor to get work done. In <strong>single-threading mode</strong>, a program only sends a single thread at a time to the processor and waits for it to be finished before sending another one. In <strong>multi-threading mode</strong>, the program may send several threads to the processor at the same time.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/2e/0f426b616b548a3f11b6928a20a324.png" />
<p>Not making sense? Let's go back to our earlier example, in which we wanted to count the number of words in each line of our text document. Hopefully it is clear that the first version of our little program (one long set of sequential instructions) is what you would run in single-threaded mode. And the second version (several subsets of instructions) is what you would run in multi-threaded mode, with each subset forming a separate thread. You would send out the first thread, which performs the preliminary work; then once it's done you would send the &quot;middle&quot; threads, which can be run in parallel; then finally once they're all done you would send out the final thread to clean up and collect final results.  </p>
<p>If you're still having a hard time visualizing what the different threads are like, just imagine that you're doing cross-stitching. If you're a regular human, you're working with just one hand. You're pulling a needle and thread (a single thread!) through the canvas, making one stitch after another, one row after another. Now try to imagine an octopus doing cross-stitching. He can make several rows of stitches at the same time using a different needle and thread for each. Multi-threading in computers is surprisingly similar to that.</p>
<p><em>Hey, if you have a better example, let us know in the forum and we'll use that instead.</em></p>
<p>Alright, now that you understand the idea of multithreading, let's get practical: how do we do get the GATK to use multi-threading?</p>
<p>There are two options for multi-threading with the GATK, controlled by the arguments <code>-nt</code> and <code>-nct</code>, respectively. They can be combined, since they act at different levels of computing:</p>
<ul>
<li>
<p><code>-nt</code> / <code>--num_threads</code> controls the number of <strong>data threads</strong> sent to the processor (acting at the <strong>machine</strong> level)</p>
</li>
<li><code>-nct</code> / <code>--num_cpu_threads_per_data_thread</code> controls the number of <strong>CPU threads</strong> allocated to each data thread (acting at the <strong>core</strong> level).</li>
</ul>
<p>Not all GATK tools can use these options due to the nature of the analyses that they perform and how they traverse the data. Even in the case of tools that are used sequentially to perform a multi-step process, the individual tools may not support the same options. For example, at time of writing (Dec. 2012), of the tools involved in local realignment around indels, RealignerTargetCreator supports <code>-nt</code> but not <code>-nct</code>, while IndelRealigner does not support either of these options. </p>
<p>In addition, there are some important technical details that affect how these options can be used with optimal results. Those are explained along with specific recommendations for the main GATK tools in a <a href="http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools">companion document</a> on parallelizing the GATK.</p>
<h3>Scatter-gather</h3>
<p>If you Google it, you'll find that the term <strong>scatter-gather</strong> can refer to a lot of different things, including strategies to get the best price quotes from online vendors, methods to control memory allocation and… an indie-rock band. What all of those things have in common (except possibly the band) is that they involve breaking up a task into smaller, parallelized tasks (scattering) then collecting and integrating the results (gathering). That should sound really familiar to you by now, since it's the general principle of parallel computing.</p>
<p>So yes, &quot;scatter-gather&quot; is really just another way to say we're parallelizing things. OK, but how is it different from multithreading, and why do we need yet another name?</p>
<p>As you know by now, multithreading specifically refers to what happens internally when the program (in our case, the GATK) sends several sets of instructions to the processor to achieve the instructions that you originally gave it in a single command-line. In contrast, the scatter-gather strategy as used by the GATK involves separate programs. There are two pipelining solutions that we support for scatter-gathering GATK jobs, Crom+WDL and Queue. They are quite different, but both are able to generate separate GATK jobs (each with its own command-line) to achieve the instructions given in a script.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/14/96d5bba42167a599a60ed7ac58602f.png" />
<p>At the simplest level, the script can involve a single GATK tool*. In that case, the execution engine (Cromwell or Queue) will create separate GATK commands that will each run that tool on a portion of the input data (= the scatter step). The results of each run will be stored in temporary files. Then once all the runs are done, the engine will collate all the results into the final output files, as if the tool had been run as a single command (= the gather step).</p>
<p><em>Note that Queue and Cromwell have additional capabilities, such as managing the use of multiple GATK tools in a dependency-aware manner to run complex pipelines, but that is outside the scope of this article. To learn more about pipelining the GATK with Queue, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1306">Queue documentation</a>. To learn more about Crom+WDL, see the <a href="https://software.broadinstitute.org/wdl/">WDL website</a>.</em></p>
<h3>Compare and combine</h3>
<p>So you see, scatter-gather is a very different process from multi-threading because the parallelization happens <strong>outside</strong> of the program itself. The big advantage is that this opens up the upper level of computing: the cluster level. Remember, the GATK program is limited to dispatching threads to the processor of the machine on which it is run – it cannot by itself send threads to a different machine. But an execution engine like Queue or Cromwell can dispatch scattered GATK jobs to different machines in a computing cluster or on a cloud platform by interfacing with the appropriate job management software.</p>
<p>That being said, multithreading has the great advantage that cores and machines all have access to shared machine memory with very high bandwidth capacity. In contrast, the multiple machines on a network used for scatter-gather are fundamentally limited by network costs.  </p>
<p>The good news is that you can combine scatter-gather and multithreading: use Queue or Cromwell to scatter GATK jobs to different nodes on your cluster or cloud platform, then use the GATK's internal multithreading capabilities to parallelize the jobs running on each node.</p>
<p>Going back to the rice-cooking example, it's as if instead of cooking the rice yourself, you hired a catering company to do it for you. The company assigns the work to several people, who each have their own cooking station with multiple rice cookers. Now you can feed a lot more people in the same amount of time! And you don't even have to clean the dishes. </p>