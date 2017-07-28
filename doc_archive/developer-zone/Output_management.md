## Output management

http://gatkforums.broadinstitute.org/gatk/discussion/1327/output-management

<h3>1. Introduction</h3>
<p>When running either single-threaded or in shared-memory parallelism mode, the GATK guarantees that output written to an output stream created via the <code>@Argument</code> mechanism will ultimately be assembled in genomic order.  In order to assemble the final output file, the GATK will write the output generated from each thread into a temporary output file, ultimately assembling the data via a central coordinating thread.  There are three major elements in the GATK that facilitate this functionality:</p>
<ul>
<li>
<p>Stub</p>
<p>The front-end interface to the output management system.  Stubs will be injected into the walker by the command-line argument system and relay information from the walker to the output management system.  There will be one stub per invocation of the GATK.</p>
</li>
<li>
<p>Storage</p>
<p>The back end interface, responsible for creating, writing and deleting temporary output files as well as merging their contents back into the primary output file.  One Storage object will exist per shard processed in the GATK.</p>
</li>
<li>
<p>OutputTracker</p>
<p>The dispatcher; ultimately connects the stub object's output creation request back to the most appropriate storage object to satisfy that request.  One OutputTracker will exist per GATK invocation.</p>
</li>
</ul>
<h3>2. Basic Mechanism</h3>
<p>Stubs are directly injected into the walker through the GATK's command-line argument parser as a go-between from walker to output management system.  When a walker calls into the stub it's first responsibility is to call into the output tracker to retrieve an appropriate storage object.  The behavior of the OutputTracker from this point forward depends mainly on the parallelization mode of this traversal of the GATK.</p>
<h4>If the traversal is single-threaded:</h4>
<ul>
<li>
<p>the OutputTracker (implemented as DirectOutputTracker) will create the storage object if necessary and return it to the stub.</p>
</li>
<li>
<p>The stub will forward the request to the provided storage object.  </p>
</li>
<li>At the end of the traversal, the microscheduler will request that the OutputTracker finalize and close the file.</li>
</ul>
<h4>If the traversal is multi-threaded using shared-memory parallelism:</h4>
<ul>
<li>
<p>The OutputTracker (implemented as ThreadLocalOutputTracker) will look for a storage object associated with this thread via a ThreadLocal.  </p>
</li>
<li>
<p>If no such storage object exists, it will be created pointing to a temporary file.  </p>
</li>
<li>
<p>At the end of <strong>each shard processed</strong>, that file will be closed and an OutputMergeTask will be created so that the shared-memory parallelism code can merge the output at its leisure.</p>
</li>
<li>The shared-memory parallelism code will merge when a fixed number of temporary files appear in the input queue.  The constant used to determine this frequency is fixed at compile time (see <code>HierarchicalMicroScheduler.MAX_OUTSTANDING_OUTPUT_MERGES</code>).</li>
</ul>
<h3>3. Using output management</h3>
<p>To use the output management system, declare a field in your walker of one of the existing core output types, coupled with either an <code>@Argument</code> or <code>@Output</code> annotation.</p>
<pre><code class="pre_md">@Output(doc="Write output to this BAM filename instead of STDOUT")
SAMFileWriter out;</code class="pre_md"></pre>
<p>Currently supported output types are SAM/BAM (declare SAMFileWriter), VCF (declare VCFWriter), and any non-buffering stream extending from OutputStream.</p>
<h3>4. Implementing a new output type</h3>
<p>To create a new output type, three types must be implemented: Stub, Storage, and ArgumentTypeDescriptor.</p>
<h4>To implement Stub</h4>
<p>Create a new Stub class, extending/inheriting the core output type's interface and implementing the Stub interface.</p>
<pre><code class="pre_md">OutputStreamStub extends OutputStream implements Stub&lt;OutputStream&gt; {</code class="pre_md"></pre>
<p>Implement a register function so that the engine can provide the stub with the session's OutputTracker.</p>
<pre><code class="pre_md">public void register( OutputTracker outputTracker ) {
    this.outputTracker = outputTracker;
}</code class="pre_md"></pre>
<p>Add as fields any parameters necessary for the storage object to create temporary storage.</p>
<pre><code class="pre_md">private final File targetFile;
public File getOutputFile() { return targetFile; }</code class="pre_md"></pre>
<p>Implement/override every method in the core output type's interface to pass along calls to the appropriate storage object via the OutputTracker.</p>
<pre><code class="pre_md">public void write( byte[] b, int off, int len ) throws IOException {
    outputTracker.getStorage(this).write(b, off, len);
}</code class="pre_md"></pre>
<h4>To implement Storage</h4>
<p>Create a Storage class, again extending inheriting the core output type's interface and implementing the Storage interface.</p>
<pre><code class="pre_md">public class OutputStreamStorage extends OutputStream implements Storage&lt;OutputStream&gt; {</code class="pre_md"></pre>
<p>Implement constructors that will accept just the Stub or Stub + alternate file path and create a repository for data, and a close function that will close that repository.</p>
<pre><code class="pre_md">public OutputStreamStorage( OutputStreamStub stub ) { ... }
public OutputStreamStorage( OutputStreamStub stub, File file ) { ... }
public void close() { ... }</code class="pre_md"></pre>
<p>Implement a <code>mergeInto</code> function capable of reconstituting the file created by the constructor, dumping it back into the core output type's interface, and removing the source file.</p>
<pre><code class="pre_md">public void mergeInto( OutputStream targetStream ) { ... }</code class="pre_md"></pre>
<p>Add a block to <code>StorageFactory.createStorage()</code> capable of creating the new storage object.  <strong>TODO: use reflection to generate the storage classes.</strong></p>
<pre><code class="pre_md">    if(stub instanceof OutputStreamStub) {
        if( file != null )
            storage = new OutputStreamStorage((OutputStreamStub)stub,file);
        else
            storage = new OutputStreamStorage((OutputStreamStub)stub);
    }</code class="pre_md"></pre>
<h4>To implement ArgumentTypeDescriptor</h4>
<p>Create a new object inheriting from type <code>ArgumentTypeDescriptor</code>.  Note that the <code>ArgumentTypeDescriptor</code> does NOT need to support the core output type's interface.</p>
<pre><code class="pre_md">public class OutputStreamArgumentTypeDescriptor extends ArgumentTypeDescriptor {</code class="pre_md"></pre>
<p>Implement a truth function indicating which types this <code>ArgumentTypeDescriptor</code> can service.</p>
<pre><code class="pre_md"> @Override
 public boolean supports( Class type ) {
     return SAMFileWriter.class.equals(type) || StingSAMFileWriter.class.equals(type);
 }</code class="pre_md"></pre>
<p>Implement a parse function that constructs the new Stub object.  The function should register this type as an output by caling <code>engine.addOutput(stub)</code>.</p>
<pre><code class="pre_md"> public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches )  {
     ...
     OutputStreamStub stub = new OutputStreamStub(new File(fileName));
     ...
     engine.addOutput(stub);
     ....
     return stub;
}</code class="pre_md"></pre>
<p>Add a creator for this new ArgumentTypeDescriptor in <code>CommandLineExecutable.getArgumentTypeDescriptors()</code>.</p>
<pre><code class="pre_md"> protected Collection&lt;ArgumentTypeDescriptor&gt; getArgumentTypeDescriptors() {
     return Arrays.asList( new VCFWriterArgumentTypeDescriptor(engine,System.out,argumentSources),
                           new SAMFileWriterArgumentTypeDescriptor(engine,System.out),
                           new OutputStreamArgumentTypeDescriptor(engine,System.out) );
 }</code class="pre_md"></pre>
<p>After creating these three objects, the new output type should be ready for usage as described above.</p>
<h3>5. Outstanding issues</h3>
<ul>
<li>
<p>Only non-buffering iterators are currently supported by the GATK.  Of particular note, <code>PrintWriter</code> will appear to drop records if created by the command-line argument system; use <code>PrintStream</code> instead.</p>
</li>
<li>For efficiency, the GATK does not reduce output files together following the tree pattern used by shared-memory parallelism; output merges happen via an independent queue.  Because of this, output merges happening during a <code>treeReduce</code> may not behave correctly.</li>
</ul>