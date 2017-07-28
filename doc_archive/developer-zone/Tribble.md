## Tribble

http://gatkforums.broadinstitute.org/gatk/discussion/1349/tribble

<h2>1. Overview</h2>
<p>The Tribble project was started as an effort to overhaul our reference-ordered data system; we had many different formats that were shoehorned into a common framework that didn't really work as intended.  What we wanted was a common framework that allowed for searching of reference ordered data, regardless of the underlying type.  Jim Robinson had developed indexing schemes for text-based files, which was incorporated into the Tribble library.</p>
<h2>2. Architecture Overview</h2>
<p>Tribble provides a lightweight interface and API for querying features and creating indexes from feature files, while allowing iteration over know feature files that we're unable to create indexes for.   The main entry point for external users is the BasicFeatureReader class. It takes in a codec, an index file, and a file containing the features to be processed.  With an instance of a <code>BasicFeatureReader</code>, you can query for features that span a specific location, or get an iterator over all the records in the file. </p>
<h2>3. Developer Overview</h2>
<p>For developers, there are two important classes to implement: the FeatureCodec, which decodes lines of text and produces features, and the feature class, which is your underlying record type.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/cc/f41b5df64878ee361ba5e4b78047ce.png" />
<p>For developers there are two classes that are important:</p>
<ul>
<li>
<p><strong>Feature</strong></p>
<p>This is the genomicly oriented feature that represents the underlying data in the input file. For instance in the VCF format, this is the variant call including quality information, the reference base, and the alternate base.  The required information to implement a feature is the chromosome name, the start position (one based), and the stop position.  The start and stop position represent a closed, one-based interval.  I.e. the first base in chromosome one would be chr1:1-1. </p>
</li>
<li>
<p><strong>FeatureCodec</strong>  </p>
<p>This class takes in a line of text (from an input source, whether it's a file, compressed file, or a http link), and produces the above feature. </p>
</li>
</ul>
<p>To implement your new format into Tribble, you need to implement the two above classes (in an appropriately named subfolder in the Tribble check-out).  The Feature object should know nothing about the file representation; it should represent the data as an in-memory object.  The interface for a feature looks like:</p>
<pre><code class="pre_md">public interface Feature {

    /**
     * Return the features reference sequence name, e.g chromosome or contig
     */
    public String getChr();

    /**
     * Return the start position in 1-based coordinates (first base is 1)
     */
    public int getStart();

    /**
     * Return the end position following 1-based fully closed conventions.  The length of a feature is
     * end - start + 1;
     */
    public int getEnd();
}</code class="pre_md"></pre>
<p>And the interface for FeatureCodec:</p>
<pre><code class="pre_md">/**
 * the base interface for classes that read in features.
 * @param &lt;T&gt; The feature type this codec reads
 */
public interface FeatureCodec&lt;T extends Feature&gt; {
    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(String line);

    /**
     * Decode a line as a Feature.
     *
     * @param line the input line to decode
     * @return  Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public T decode(String line);

    /**
     * This function returns the object the codec generates.  This is allowed to be Feature in the case where
     * conditionally different types are generated.  Be as specific as you can though.
     *
     * This function is used by reflections based tools, so we can know the underlying type
     *
     * @return the feature type this codec generates.
     */
    public Class&lt;T&gt; getFeatureType();

    /**  Read and return the header, or null if there is no header.
     *
     * @return header object
     */
    public Object readHeader(LineReader reader);
}</code class="pre_md"></pre>
<h2>4. Supported Formats</h2>
<p>The following formats are supported in Tribble:</p>
<ul>
<li>VCF Format</li>
<li>DbSNP Format</li>
<li>BED Format</li>
<li>GATK Interval Format</li>
</ul>
<h2>5. Updating the Tribble, htsjdk, and/or Picard library</h2>
<p>Updating the revision of Tribble on the system is a relatively straightforward task if the following steps are taken.</p>
<p><em>NOTE:</em> Any directory starting with <code>~</code> may be different on your machine, depending on where you cloned the various repositories for gsa-unstable, picard, and htsjdk.</p>
<p>A Maven script to install picard into the local repository is located under <code>gsa-unstable/private/picard-maven</code>. To operate, it requires a symbolic link named <code>picard</code> pointing to a working checkout of the <a href="http://github.com/broadinstitute/picard">picard github repository</a>. <em>NOTE:</em> compiling picard <a href="http://broadinstitute.github.io/picard">requires</a> an <a href="http://github.com/samtools/htsjdk">htsjdk github repository</a> checkout available at <code>picard/htsjdk</code>, either as a subdirectory or another symbolic link. The final full path should be <code>gsa-unstable/private/picard-maven/picard/htsjdk</code>.</p>
<pre><code class="pre_md">cd ~/src/gsa-unstable
cd private/picard-maven
ln -s ~/src/picard picard</code class="pre_md"></pre>
<p>Create a git branch of Picard and/or htsjdk and make your changes. To install your changes into the GATK you must run <code>mvn install</code> in the <code>private/picard-maven</code> directory. This will compile and copy the jars into <code>gsa-unstable/public/repo</code>, and update <code>gsa-unstable/gatk-root/pom.xml</code> with the corresponding version. While making changes your revision of picard and htslib will be labeled with <code>-SNAPSHOT</code>.</p>
<pre><code class="pre_md">cd ~/src/gsa-unstable
cd private/picard-maven
mvn install</code class="pre_md"></pre>
<p>Continue testing in the GATK. Once your changes and updated tests for picard/htsjdk are complete, push your branch and submit your pull request to the Picard and/or htsjdk github. After your Picard/htsjdk patches are accepted, switch your Picard/htsjdk branches back to the master branch. <em>NOTE:</em> Leave your gsa-unstable branch on your development branch!</p>
<pre><code class="pre_md">cd ~/src/picard
ant clean
git checkout master
git fetch
git rebase
cd htsjdk
git checkout master
git fetch
git rebase</code class="pre_md"></pre>
<p><em>NOTE:</em> The version number of old and new Picard/htsjdk will vary, and during active development will end with <code>-SNAPSHOT</code>. While, if needed, you may push <code>-SNAPSHOT</code> version for testing on Bamboo, you should NOT submit a pull request with a <code>-SNAPSHOT</code> version. <code>-SNAPSHOT</code> indicates your local changes are not reproducible from source control.</p>
<p>When ready, run <code>mvn install</code> once more to create the non <code>-SNAPSHOT</code> versions under <code>gsa-unstable/public/repo</code>. In that directory, <code>git add</code> the new version, and <code>git rm</code> the old versions.</p>
<pre><code class="pre_md">cd ~/src/gsa-unstable
cd public/repo
git add picard/picard/1.115.1499/
git add samtools/htsjdk/1.115.1509/
git rm -r picard/picard/1.112.1452/
git rm -r samtools/htsjdk/1.112.1452/</code class="pre_md"></pre>
<p>Commit and then push your gsa-unstable branch, then issue a pull request for review.</p>