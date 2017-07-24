## GATK development process and coding standards

http://gatkforums.broadinstitute.org/gatk/discussion/2129/gatk-development-process-and-coding-standards

<h2>Introduction</h2>
<p>This document describes the current GATK coding standards for documentation and unit testing.  The overall goal is that all functions be well documented, have unit tests, and conform to the coding conventions described in this guideline. It is primarily meant as an internal reference for team members, but we are making it public to provide an example of how we work. There are a few mentions of specific team member responsibilities and who to contact with questions; please just disregard those as they will not be applicable to you.</p>
<h2>Coding conventions</h2>
<h3>General conventions</h3>
<p>The Genome Analysis Toolkit generally follows Java coding standards and good practices, which can be viewed <a href="http://www.oracle.com/technetwork/java/codeconvtoc-136057.html">at Sun's site</a>. </p>
<p>The original coding standard document for the GATK was developed in early 2009.  It remains a reasonable starting point but may be superseded by statements on this page (<a href="https://us.v-cdn.net/5019796/uploads/FileUpload/18/a199e46fbc5c5e08866e8136db7192.pdf">available as a PDF</a>).</p>
<h3>Size of functions and functional programming style</h3>
<p>Code in the GATK should be structured into clear, simple, and testable functions.  Clear means that the function takes a limited number of arguments, most of which are values not modified, and in general should return newly allocated results, as opposed to directly modifying the input arguments (functional style).  The max. size of functions should be approximately one screen's worth of real estate (no more than 80 lines), including inline comments.  If you are writing functions that are much larger than this, you must refactor your code into modular components.</p>
<h3>Code duplication</h3>
<p>Do not duplicate code.  If you are finding yourself wanting to make a copy of functionality, refactor the code you want to duplicate and enhance it.  Duplicating code introduces bugs, makes the system harder to maintain, and will require more work since you will have a new function that must be tested, as opposed to expanding the tests on the existing functionality.</p>
<h3>Documentation</h3>
<p>Functions must be documented following the javadoc conventions.  That means that the first line of the comment should be a simple statement of the purpose of the function.  Following that is an expanded description of the function, such as edge case conditions, requirements on the argument, state changes, etc.  Finally comes the @param and @return fields, that should describe the meaning of each function argument, restrictions on the values allowed or returned.  In general, the return field should be about types and ranges of those values, not the meaning of the result, as this should be in the body of the documentation.</p>
<h3>Testing for valid inputs and contracts</h3>
<p>The GATK uses Contracts for Java to help us enforce code quality during testing.  See <a href="http://code.google.com/p/cofoja/">CoFoJa</a> for more information.  If you've never programmed with contracts, read their excellent description <a href="http://code.google.com/p/cofoja/wiki/AddContracts">Adding contracts to a stack</a>.  Contracts are only enabled when we are testing the code (unittests and integration tests) and not during normal execution, so contracts can be reasonably expensive to compute.  They are best used to enforce assumptions about the status of class variables and return results.  </p>
<p>Contracts are tricky when it comes to input arguments.  The best practice is simple:</p>
<ul>
<li>Public functions with arguments should explicitly test those input arguments for good values with live java code (such as in the example below).  Because the function is public, you don't know what the caller will be passing in, so you have to check and ensure quality.</li>
<li>Private functions with arguments should use contracts instead.  Because the function is private, the author of the code controls use of the function, and the contracts enforce good use.  But in principal the quality of the inputs should be assumed at runtime since only the author controlled calls to the function and input QC should have happened elsewhere</li>
</ul>
<p>Below is an example private function that makes good use of input argument contracts:</p>
<pre><code class="pre_md">/**
 * Helper function to write out a IGV formatted line to out, at loc, with values
 *
 * http://www.broadinstitute.org/software/igv/IGV
 *
 * @param out a non-null PrintStream where we'll write our line
 * @param loc the location of values
 * @param featureName string name of this feature (see IGV format)
 * @param values the floating point values to associate with loc and feature name in out
 */
@Requires({
        "out != null",
        "loc != null",
        "values.length &gt; 0"
})
private void printIGVFormatRow(final PrintStream out, final GenomeLoc loc, final String featureName, final double ... values) {
    // note that start and stop are 0 based, but the stop is exclusive so we don't subtract 1
    out.printf("%s\t%d\t%d\t%s", loc.getContig(), loc.getStart() - 1, loc.getStop(), featureName);
    for ( final double value : values )
        out.print(String.format("\t%.3f", value));
    out.println();
} </code class="pre_md"></pre>
<h3>Final variables</h3>
<p>Final java fields cannot be reassigned once set.  Nearly all variables you write should be final, unless they are obviously accumulator results or other things you actually want to modify.  Nearly all of your function arguments should be final.  Being final stops incorrect reassigns (a major bug source) as well as more clearly captures the flow of information through the code. </p>
<h3>An example high-quality GATK function</h3>
<pre><code class="pre_md">/**
 * Get the reference bases from referenceReader spanned by the extended location of this active region,
 * including additional padding bp on either side.  If this expanded region would exceed the boundaries
 * of the active region's contig, the returned result will be truncated to only include on-genome reference
 * bases
 * @param referenceReader the source of the reference genome bases
 * @param padding the padding, in BP, we want to add to either side of this active region extended region
 * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
 * @return a non-null array of bytes holding the reference bases in referenceReader
 */
@Ensures("result != null")
public byte[] getReference( final IndexedFastaSequenceFile referenceReader, final int padding, final GenomeLoc genomeLoc ) {
    if ( referenceReader == null ) throw new IllegalArgumentException("referenceReader cannot be null");
    if ( padding &lt; 0 ) throw new IllegalArgumentException("padding must be a positive integer but got " + padding);
    if ( genomeLoc == null ) throw new IllegalArgumentException("genomeLoc cannot be null");
    if ( genomeLoc.size() == 0 ) throw new IllegalArgumentException("GenomeLoc must have size &gt; 0 but got " + genomeLoc);

    final byte[] reference =  referenceReader.getSubsequenceAt( genomeLoc.getContig(),
            Math.max(1, genomeLoc.getStart() - padding),
            Math.min(referenceReader.getSequenceDictionary().getSequence(genomeLoc.getContig()).getSequenceLength(), genomeLoc.getStop() + padding) ).getBases();

    return reference;
}</code class="pre_md"></pre>
<h2>Unit testing</h2>
<p>All classes and methods in the GATK should have unit tests to ensure that they work properly, and to protect yourself and others who may want to extend, modify, enhance, or optimize you code.  That GATK development team assumes that anything that isn't unit tested is broken.  Perhaps right now they aren't broken, but with a team of 10 people they will become broken soon if you don't ensure they are correct going forward with unit tests.</p>
<p>Walkers are a particularly complex issue.  UnitTesting the map and reduce results is very hard, and in my view largely unnecessary.  That said, you should write your walkers and supporting classes in such a way that all of the complex data processing functions are separated from the map and reduce functions, and those should be unit tested properly.  </p>
<p>Code coverage tells you how much of your class, at the statement or function level, has unit testing coverage.  The GATK development standard is to reach something &gt;80% method coverage (and ideally &gt;80% statement coverage).  The target is flexible as some methods are trivial (they just call into another method) so perhaps don't need coverage.  At the statement level, you get deducted from 100% for branches that check for things that perhaps you don't care about, such as illegal arguments, so reaching 100% statement level coverage is unrealistic for most clases.</p>
<p>You can find out more information about generating code coverage results at <a href="http://gatkforums.broadinstitute.org/discussion/2002/clover-coverage-analysis-with-ant#latest">Analyzing coverage with clover</a> </p>
<p>We've created a unit testing example template in the GATK codebase that provides examples of creating core GATK data structures from scratch for unit testing.  The code is in class ExampleToCopyUnitTest and can be viewed here in github directly <a href="https://github.com/broadinstitute/gsa-unstable/blob/master/public/java/test/org/broadinstitute/sting/ExampleToCopyUnitTest.java">ExampleToCopyUnitTest</a>.</p>
<h2>The GSA-Workflow</h2>
<p>As of GATK 2.5, we are moving to a full code review process, which has the following benefits:</p>
<ul>
<li>Reducing obvious coding bugs seen by other eyes</li>
<li>Reducing code duplication, as reviewers will be able to see duplicated code within the commit and potentially across the codebase</li>
<li>Ensure that coding quality standards are met (style and unit testing)</li>
<li>Setting a higher code quality standard for the master GATK unstable branch</li>
<li>Providing detailed coding feedback to newer developers, so they can improve their skills over time</li>
</ul>
<h3>The GSA workflow in words :</h3>
<ul>
<li>Create a new branch to start any work. Never work on master.
<ul>
<li>branch names have to follow the convention of [author prefix]<em>[feature name]</em>[JIRA ticket] (e.g. rp_pairhmm_GSA-232)</li>
</ul></li>
<li>Make frequent commits.</li>
<li>Push frequently your branch to origin (branch -&gt; branch)</li>
<li>When you're done -- rewrite your commit history to tell a compelling story <a href="http://git-scm.com/book/en/Git-Tools-Rewriting-History">Git Tools Rewriting History</a></li>
<li>Push your rewritten history, and request a code review. 
<ul>
<li>The entire GSA team will review your code</li>
<li>Mark DePristo assigns the reviewer responsible for making the judgment based on all reviews and merge your code into master. </li>
</ul></li>
<li>If your pull-request gets rejected, follow the comments from the team to fix it and repeat the workflow until you're ready to submit a new pull request.</li>
<li>If your pull-request is accepted, the reviewer will merge and remove your remote branch.</li>
</ul>
<h3>Example GSA workflow in the command line:</h3>
<pre><code class="pre_md"># starting a new feature
git checkout -b rp_pairhmm_GSA-332
git commit -av 
git push -u origin rp_pairhmm_GSA-332

# doing work on existing feature
git commit -av
git push

# ready to submit pull-request
git fetch origin
git rebase -i origin/master
git push -f

# after being accepted, delete your branch
git checkout master 
git pull
git branch -d rp_pairhmm_GSA-332
(the reviewer will remove your github branch)</code class="pre_md"></pre>
<h3>Commit histories and rebasing</h3>
<p>You must commit your code in small commit blocks with commit messages that follow the git best practices, which require the first line of the commit to summarize the purpose of the commit, followed by -- lines that describe the changes in more detail.  For example, here's a recent commit that meets this criteria that added unit tests to the GenomeLocParser:</p>
<pre><code class="pre_md">Refactoring and unit testing GenomeLocParser

-- Moved previously inner class to MRUCachingSAMSequenceDictionary, and unit test to 100% coverage
-- Fully document all functions in GenomeLocParser
-- Unit tests for things like parsePosition (shocking it wasn't tested!)
-- Removed function to specifically create GenomeLocs for VariantContexts.  The fact that you must incorporate END attributes in the context means that createGenomeLoc(Feature) works correctly
-- Depreciated (and moved functionality) of setStart, setStop, and incPos to GenomeLoc
-- Unit test coverage at like 80%, moving to 100% with next commit</code class="pre_md"></pre>
<p>Now, git encourages you to commit code often, and develop your code in whatever order or what is best for you.  So it's common to end up with 20 commits, all with strange, brief commit messages, that you want to push into the master branch.  It is not acceptable to push such changes.  You need to use the git command rebase to reorganize your commit history so satisfy the small number of clear commits with clear messages.  </p>
<p>Here is a recommended git workflow using rebase:</p>
<ol>
<li>
<p>Start every project by creating a new branch for it. From your master branch, type the following command (replacing &quot;myBranch&quot; with an appropriate name for the new branch):</p>
<pre><code class="pre_md">git checkout -b myBranch</code class="pre_md"></pre>
<p>Note that you only include the <em>-b</em> when you're first creating the branch. After a branch is already created, you can switch to it by typing the checkout command without the <em>-b</em>: &quot;git checkout myBranch&quot;</p>
<p>Also note that since you're always starting a new branch from master, you should keep your master branch up-to-date by occasionally doing a &quot;git pull&quot; while your master branch is checked out. You shouldn't do any actual work on your master branch, however.</p>
</li>
<li>
<p>When you want to update your branch with the latest commits from the central repo, type this while your branch is checked out:</p>
<pre><code class="pre_md">git fetch &amp;&amp; git rebase origin/master</code class="pre_md"></pre>
<p>If there are conflicts while updating your branch, git will tell you what additional commands to use.</p>
<p>If you need to combine or reorder your commits, add &quot;-i&quot; to the above command, like so:</p>
<pre><code class="pre_md">git fetch &amp;&amp; git rebase -i origin/master</code class="pre_md"></pre>
<p>If you want to edit your commits without also retrieving any new commits, omit the &quot;git fetch&quot; from the above command.</p>
</li>
</ol>
<p>If you find the above commands cumbersome or hard to remember, create aliases for them using the following commands:</p>
<pre><code class="pre_md">    git config --global alias.up '!git fetch &amp;&amp; git rebase origin/master'
    git config --global alias.edit '!git fetch &amp;&amp; git rebase -i origin/master'
    git config --global alias.done '!git push origin HEAD:master'</code class="pre_md"></pre>
<p>Then you can type &quot;git up&quot; to update your branch, &quot;git edit&quot; to combine/reorder commits, and &quot;git done&quot; to push your branch.</p>
<p>Here are more useful tutorials on how to use rebase:</p>
<ul>
<li><a href="http://git-scm.com/book/en/Git-Tools-Rewriting-History">Git Tools Rewriting History</a></li>
<li><a href="http://www.reviewboard.org/docs/codebase/dev/git/clean-commits/">Keeping commit histories clean</a></li>
<li><a href="http://darwinweb.net/articles/the-case-for-git-rebase">The case for git rebase</a></li>
<li><a href="http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html">Squashing commits with rebase</a></li>
</ul>
<p>If you need help with rebasing, talk to Mauricio or David and they will help you out.</p>