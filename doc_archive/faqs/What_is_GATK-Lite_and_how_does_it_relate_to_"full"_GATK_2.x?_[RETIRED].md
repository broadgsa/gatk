## What is GATK-Lite and how does it relate to "full" GATK 2.x? [RETIRED]

http://gatkforums.broadinstitute.org/gatk/discussion/1720/what-is-gatk-lite-and-how-does-it-relate-to-full-gatk-2-x-retired

<p><strong>Please note that GATK-Lite was retired in February 2013 when version 2.4 was released. See the announcement <a href="http://www.broadinstitute.org/gatk/guide/article?id=2091">here</a>.</strong></p>
<hr />
<p>You probably know by now that GATK-Lite is a free-for-everyone and completely open-source version of the GATK (licensed under the original [MIT license]( <a href="http://en.wikipedia.org/wiki/MIT_License">http://en.wikipedia.org/wiki/MIT_License</a>)). </p>
<p>But what's in the box? What can GATK-Lite do -- or rather, what can it <strong>not</strong> do that the full version (let's call it GATK-Full) can? And what does that mean exactly, in terms of functionality, reliability and power?  </p>
<p>To really understand the differences between GATK-Lite and GATK-Full, you need some more information on how the GATK works, and how we work to develop and improve it.</p>
<h3>First you need to understand what are the two core components of the GATK: the engine and tools (see picture below).</h3>
<p>As explained <a href="http://www.broadinstitute.org/gatk/about/#what-is-the-gatk">here</a>, the <strong>engine</strong> handles all the common work that's related to data access, conversion and traversal, as well as high-performance computing features. The engine is supported by an infrastructure of software libraries. If the GATK was a car, that would be the engine and chassis. What we call the *<em>tools</em> are attached on top of that, and they provide the various analytical and processing functionalities like variant calling and base or variant recalibration. On your car, that would be headlights, airbags and so on.</p>
<p><img src="http://www.broadinstitute.org/gatk/img/core_gatk2.png" alt="Core GATK components" /></p>
<h3>Second is how we work on developing the GATK, and what it means for how improvements are shared (or not) between Lite and Full.</h3>
<p>We do all our development work on a single codebase. This means that everything --the engine and all tools-- is on one common workbench. There are <strong>not</strong> different versions that we work on in parallel -- that would be crazy to manage! That's why the version numbers of GATK-Lite and GATK-Full always match: if the latest GATK-Full version is numbered 2.1-13, then the latest GATK-Lite is also numbered 2.1-13.</p>
<p>The most important consequence of this setup is that when we make improvements to the infrastructure and engine, the same improvements will end up in GATK Lite and in GATK Full. So for the purposes of power, speed and robustness of the GATK that is determined by the engine, there is no difference between them. </p>
<p>For the tools, it's a little more complicated -- but not much. When we &quot;build&quot; the GATK binaries (the <code>.jar</code> files), we put everything from the workbench into the Full build, but we only put a subset into the Lite build. Note that this Lite subset is pretty big -- it contains all the tools that were previously available in GATK 1.x versions, and always will. We also  reserve the right to add previews or not-fully-featured versions of the new tools that are in Full, at our discretion, to the Lite build.</p>
<h3>So there are two basic types of differences between the tools available in the Lite and Full builds (see picture below).</h3>
<ol>
<li>
<p>We have a new tool that performs a brand new function (which wasn't available in GATK 1.x), and we only include it in the Full build.</p>
</li>
<li>We have a tool that has some new add-on capabilities (which weren't possible in GATK 1.x); we put the tool in both the Lite and the Full build, but the add-ons are only available in the Full build.</li>
</ol>
<p><img src="http://www.broadinstitute.org/gatk/img/lite_vs_2x.png" alt="Tools in Lite vs. Full" /></p>
<p>Reprising the car analogy, GATK-Lite and GATK-Full are like two versions of the same car -- the basic version and the fully-equipped one. They both have the exact same engine, and most of the equipment (tools) is the same -- for example, they both have the same airbag system, and they both have headlights. But there are a few important differences: </p>
<ol>
<li>
<p>The GATK-Full car comes with a GPS (sat-nav for our UK friends), for which the Lite car has no equivalent. You could buy a portable GPS unit from a third-party store for your Lite car, but it might not be as good, and certainly not as convenient, as the Full car's built-in one.</p>
</li>
<li>Both cars have windows of course, but the Full car has power windows, while the Lite car doesn't. The Lite windows can open and close, but you have to operate them by hand, which is much slower. </li>
</ol>
<h3>So, to summarize:</h3>
<p>The underlying engine is exactly the same in both GATK-Lite and GATK-Full. Most functionalities are available in both builds, performed by the same tools. Some functionalities are available in both builds, but they are performed by different tools, and the tool in the Full build is better. New, cutting-edge functionalities are only available in the Full build, and there is no equivalent in the Lite build. </p>
<p>We hope this clears up some of the confusion surrounding GATK-Lite. If not, please leave a comment and we'll do our best to clarify further! </p>