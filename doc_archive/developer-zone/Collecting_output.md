## Collecting output

http://gatkforums.broadinstitute.org/gatk/discussion/1341/collecting-output

<h2>1. Analysis output overview</h2>
<p>In theory, any class implementing the <code>OutputStream</code> interface.  In practice, three types of classes are commonly used: <code>PrintStreams</code> for plain text files, <code>SAMFileWriters</code> for BAM files, and <code>VCFWriters</code> for VCF files.</p>
<h2>2. PrintStream</h2>
<p>To declare a basic <code>PrintStream</code> for output, use the following declaration syntax:</p>
<pre><code class="pre_md">@Output
public PrintStream out;</code class="pre_md"></pre>
<p>And use it just as you would any other PrintStream:</p>
<pre><code class="pre_md">out.println("Hello, world!");</code class="pre_md"></pre>
<p>By default, <code>@Output</code> streams prepopulate <code>fullName</code>, <code>shortName</code>, <code>required</code>, and <code>doc</code>.  <code>required</code> in this context means that the GATK will always fill in the contents of the <code>out</code> field for you.  If the user specifies no <code>--out</code> command-line argument, the 'out' field will be prepopulated with a stream pointing to <code>System.out</code>.</p>
<p>If your walker outputs a custom format that requires more than simple concatenation by [Queue]() you should also implement a custom <code>Gatherer</code>.</p>
<h2>3. SAMFileWriter</h2>
<p>For some applications, you might need to manage their own SAM readers and writers directly from inside your walker.  Current best practice for creating these Readers / Writers is to declare arguments of type <code>SAMFileReader</code> or <code>SAMFileWriter</code> as in the following example:</p>
<pre><code class="pre_md">@Output
SAMFileWriter outputBamFile = null;</code class="pre_md"></pre>
<p>If you do not specify the full name and short name, the writer will provide system default names for these arguments.  Creating a <code>SAMFileWriter</code> in this way will create the type of writer most commonly used by members of the GSA group at the Broad Institute -- it will use the same header as the input BAM and require presorted data.  To change either of these attributes, use the <code>StingSAMIterator</code> interface instead:</p>
<pre><code class="pre_md">@Output
StingSAMFileWriter outputBamFile = null;</code class="pre_md"></pre>
<p>and later, in <code>initialize()</code>, run one or both of the following methods:</p>
<p>outputBAMFile.writeHeader(customHeader);
outputBAMFile.setPresorted(false);</p>
<p>You can change the header or presorted state until the first alignment is written to the file.</p>
<h2>4. VCFWriter</h2>
<p><code>VCFWriter</code> outputs behave similarly to <code>PrintStreams</code> and <code>SAMFileWriters</code>.  Declare a <code>VCFWriter</code> as follows:</p>
<p>@Output(doc=&quot;File to which variants should be written&quot;,required=true)
protected VCFWriter writer = null;</p>
<h2>5. Debugging Output</h2>
<p>The walkers provide a protected logger instance. Users can adjust the debug level of the walkers using the <code>-l</code> command line option.</p>
<p>Turning on verbose logging can produce more output than is really necessary. To selectively turn on logging for a class or package, specify a <code>log4j.properties</code> property file from the command line as follows:</p>
<pre><code class="pre_md">-Dlog4j.configuration=file:///&lt;your development root&gt;/Sting/java/config/log4j.properties</code class="pre_md"></pre>
<p>An example <code>log4j.properties</code> file is available in the <code>java/config</code> directory of the Git repository.</p>