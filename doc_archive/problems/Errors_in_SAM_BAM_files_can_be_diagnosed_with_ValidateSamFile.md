## Errors in SAM/BAM files can be diagnosed with ValidateSamFile

http://gatkforums.broadinstitute.org/gatk/discussion/7571/errors-in-sam-bam-files-can-be-diagnosed-with-validatesamfile

<h3>The problem</h3>
<p>You're trying to run a GATK or Picard tool that operates on a SAM or BAM file, and getting some cryptic error that doesn't clearly tell you what's wrong. Bits of the stack trace (the pile of lines in the output log that the program outputs when there is a problem) may contain the following: <code>java.lang.String</code>, <code>Error Type Count</code>, <code>NullPointerException</code> -- or maybe something else that doesn't mean anything to you. </p>
<h3>Why this happens</h3>
<p>The most frequent cause of these unexplained problems is not a bug in the program -- it's an invalid or malformed SAM/BAM file. This means that there is something wrong either with the content of the file (something important is missing) or with its format (something is written the wrong way). Invalid SAM/BAM files generally have one or more errors in the following sections: the header tags, the alignment fields, or the optional alignment tags. In addition, the SAM/BAM index file can be a source of errors as well. </p>
<p>The source of these errors is usually introduced by upstream processing tools, such as the genome mapper/aligner or any other data processing tools you may have applied before feeding the data to <a href="http://broadinstitute.github.io/picard/">Picard</a> or <a href="https://www.broadinstitute.org/gatk/">GATK</a>.  </p>
<h3>The solution</h3>
<p>To fix these problems, you first have to know what's wrong. Fortunately there's a handy <a href="http://broadinstitute.github.io/picard/">Picard</a> tool that can test for (almost) all possible SAM/BAM format errors, called <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a>. </p>
<p>We recommend the workflow included below for diagnosing problems with <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a>. This workflow will help you tackle the problem efficiently and set priorities for dealing with multiple errors (which often happens). We also outline typical solutions for common errors, but note that this is not meant to be an exhaustive list -- there are too many possible problems to tackle all of them in this document. To be clear, here we focus on diagnostics, not treatment. </p>
<p><em>In some cases, it may not be possible to fix some problems that are too severe, and you may need to redo the genome alignment/mapping from scratch! Consider running <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> proactively at all key steps of your analysis pipeline to catch errors early!</em></p>
<hr />
<h2>Workflow for diagnosing SAM/BAM file errors with <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a></h2>
<div>
<center>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ee/b00ad090cf37733ce86ac1e5999b87.png" />
</center>
</div>
<h3>1. Generate summary of errors</h3>
<p>First, run <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> in <code>SUMMARY</code> mode in order to get a summary of everything that is missing or improperly formatted in your input file. We set <code>MODE=SUMMARY</code> explicitly because by default the tool would just emit details about the 100 first problems it finds then quit. If you have some minor formatting issues that don't really matter but affect every read record, you won't get to see more important problems that occur later in the file.</p>
<pre><code>$ java -jar picard.jar ValidateSamFile \ 
        I=input.bam \ 
        MODE=SUMMARY </code></pre>
<p>If this outputs <code>No errors found</code>, then your SAM/BAM file is completely valid. If you were running this purely as a preventative measure, then you're good to go and proceed to the next step in your pipeline. If you were doing this to diagnose a problem, then you're back to square one -- but at least now you know it's not likely to be a SAM/BAM file format issue. One exception: some analysis tools require Read Group tags like <code>SM</code> that not required by the format specification itself, so the input files will pass validation but the analysis tools will still error out. If that happens to you, check whether your files have <code>SM</code> tags in the <code>@RG</code> lines in their BAM header. That is the most common culprit. </p>
<p>However, if the command above outputs one or more of the 8 possible <code>WARNING</code> or 48 possible <code>ERROR</code> messages (see tables at the end of this document), you must proceed to the next step in the diagnostic workflow.</p>
<p>When run in <code>SUMMARY</code> mode, <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> outputs a table that differentiates between two levels of error: <code>ERROR</code> proper and <code>WARNING</code>, based on the severity of problems that they would cause in downstream analysis. All problems that fall in the <code>ERROR</code> category <strong>must</strong> be addressed to in order to proceed with other <a href="http://broadinstitute.github.io/picard/">Picard</a> or <a href="https://www.broadinstitute.org/gatk/">GATK</a> tools, while those that fall in the <code>WARNING</code> category may often be ignored for some, if not all subsequent analyses. </p>
<h4>Example of error summary</h4>
<table class="table table-striped">

<tr><th><b>ValidateSamFile (SUMMARY)    </b></th><th><b>    Count   </b></th></tr>
<tr><td>ERROR:MISSING_READ_GROUP                                       </td><td>1           </td></tr>
<tr><td>ERROR:MISMATCH_MATE_ALIGNMENT_START            </td><td>4            </td></tr>
<tr><td>ERROR:MATES_ARE_SAME_END                                      </td><td>894289  </td></tr>
<tr><td>ERROR:CIGAR_MAPS_OFF_REFERENCE                     </td><td>354        </td></tr>
<tr><td>ERROR:MATE_NOT_FOUND                                               </td><td>1            </td></tr>
<tr><td>ERROR:MISMATCH_FLAG_MATE_UNMAPPED              </td><td>46672    </td></tr>
<tr><td>ERROR:MISMATCH_READ_LENGTH_AND_E2_LENGTH  </td><td> 1           </td></tr>
<tr><td>WARNING:RECORD_MISSING_READ_GROUP              </td><td>54          </td></tr>
<tr><td>WARNING:MISSING_TAG_NM                                             </td><td>33          </td></tr>
</table>
<p>This table, generated by <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> from a real BAM file, indicates that this file has a total of 1 <code>MISSING_READ_GROUP</code> error, 4 <code>MISMATCH_MATE_ALIGNMENT_START</code> errors, 894,289 <code>MATES_ARE_SAME_END</code> errors, and so on.  Moreover, this output also indicates that there are 54 <code>RECORD_MISSING_READ_GROUP</code> warnings and 33 <code>MISSING_TAG_NM</code> warnings.  </p>
<h3>2. Generate detailed list of ERROR records</h3>
<p>Since <code>ERRORs</code> are more severe than <code>WARNINGs</code>, we focus on diagnosing and fixing them first.  From the first step we only had a summary of errors, so now we generate a more detailed report with this command:</p>
<pre><code>$ java -jar picard.jar ValidateSamFile \ 
        I=input.bam \ 
        IGNORE_WARNINGS=true \
        MODE=VERBOSE </code></pre>
<p>Note that we invoked the <code>MODE=VERBOSE</code> and the <code>IGNORE_WARNINGS=true</code> arguments. </p>
<p>The former is technically not necessary as <code>VERBOSE</code> is the tool's default mode, but we specify it here to make it clear that that's the behavior we want. This produces a complete list of every problematic record, as well as a more descriptive explanation for each type of <code>ERROR</code> than is given in the <code>SUMMARY</code> output. </p>
<p>The <code>IGNORE_WARNINGS</code> option enables us to specifically examine only the records with <code>ERRORs</code>. When working with large files, this feature can be quite helpful, because there may be many records with <code>WARNINGs</code> that are not immediately important, and we don't want them flooding the log output.  </p>
<h4>Example of VERBOSE report for ERRORs only</h4>
<table class="table table-striped">
<tr><th><b>ValidateSamFile (VERBOSE)    </b></th><th><b>    Error Description   </b></th></tr>
<tr><td>      ERROR: Read groups is empty          </td><td> Empty read group field for multiple records </td></tr>
<tr><td>      ERROR: Record 1, Read name 20FUKAAXX100202:6:27:4968:125377</td><td>Mate alignment does not match alignment start of mate</td></tr>
<tr><td>      ERROR: Record 3, Read name 20FUKAAXX100202:6:27:4986:125375 </td><td>Both mates are marked as second of pair</td></tr>
<tr><td>ERROR: Record 6, Read name 20GAVAAXX100126:4:47:18102:194445 </td><td> Read CIGAR M operator maps off end of reference </td></tr>
<tr><td>ERROR: Read name 30PPJAAXX090125:1:60:1109:517#0 </td><td>Mate not found for paired read </td></tr>
<tr><td>ERROR: Record 402, Read name 20GAVAAXX100126:3:44:17022:23968</td><td> Mate unmapped flag does not match read unmapped flag of mate</td></tr>
<tr><td>ERROR: Record 12, Read name HWI-ST1041:151:C7BJEACXX:1:1101:1128:82805</td><td> Read length does not match quals length</td></tr>
</table>
<p>These <code>ERRORs</code> are all problems that we must address before using this BAM file as input for further analysis. Most <code>ERRORs</code> can typically be fixed using <a href="http://broadinstitute.github.io/picard/">Picard</a> tools to either correct the formatting or fill in missing information, although sometimes you may want to simply filter out malformed reads using Samtools.  </p>
<p>For example, <code>MISSING_READ_GROUP</code> errors can be solved by adding the read group information to your data using the <a href="http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a> tool. Most mate pair information errors can be fixed with <a href="http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation">FixMateInformation</a>.</p>
<p>Once you have attempted to fix the errors in your file, you should put your new SAM/BAM file through the first validation step in the workflow, running <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> in <code>SUMMARY</code> mode again. We do this to evaluate whether our attempted fix has solved the original <code>ERRORs</code>, and/or any of the original <code>WARNINGs</code>, and/or introduced any <em>new</em> <code>ERRORs</code> or <code>WARNINGs</code> (sadly, this does happen). </p>
<p>If you still have <code>ERRORs</code>, you'll have to loop through this part of the workflow until no more <code>ERRORs</code> are detected. </p>
<p>If you have no more <code>ERRORs</code>, congratulations! It's time to look at the <code>WARNINGs</code> (assuming there are still some -- if not, you're off to the races).</p>
<h3>3. Generate detailed list of WARNING records</h3>
<p>To obtain more detailed information about the warnings, we invoke the following command: </p>
<pre><code>$ java -jar picard.jar ValidateSamFile \ 
        I=input.bam \ 
        IGNORE=type \
        MODE=VERBOSE </code></pre>
<p>At this time we often use the <code>IGNORE</code> option to tell the program to ignore a specific type of <code>WARNING</code> that we consider less important, in order to focus on the rest. In some cases we may even decide to not try to address some <code>WARNINGs</code> at all because we know they are harmless (for example, <code>MATE_NOT_FOUND</code> warnings are expected when working with a small snippet of data). But in general we do strongly recommend that you address all of them to avoid any downstream complications, unless you're sure you know what you're doing. </p>
<h4>Example of VERBOSE report for WARNINGs only</h4>
<table class="table table-striped">
<tr><th><b>ValidateSamFile (VERBOSE)    </b></th><th><b>    Warning Description </b></th></tr>
<tr><td>WARNING: Read name H0164ALXX140820:2:1204:13829:66057</td><td> A record is missing a read group</td></tr>
<tr><td>WARNING: Record 1, Read name HARMONIA-H16:1253:0:7:1208:15900:108776    </td><td>   NM tag (nucleotide differences) is missing  </td></tr>
</table>
<p>Here we see a read group-related <code>WARNING</code> which would probably be fixed when we fix the  <code>MISSING_READ_GROUP</code> error we encountered earlier, hence the prioritization strategy of tackling <code>ERRORs</code> first and <code>WARNINGs</code> second.</p>
<p>We also see a <code>WARNING</code> about missing <code>NM</code> tags. This is an alignment tag that is added by some but not all genome aligners, and is not used by the downstream tools that we care about, so you may decide to ignore this warning by adding <code>IGNORE=MISSING_TAG_NM</code> from now on when you run <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> on this file. </p>
<p>Once you have attempted to fix all the <code>WARNINGs</code> that you care about in your file, you put your new SAM/BAM file through the first validation step in the workflow again, running <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> in <code>SUMMARY</code> mode. Again, we check that no new <code>ERRORs</code> have been introduced and that the only <code>WARNINGs</code> that remain are the ones we feel comfortable ignoring. If that's not the case we run through the workflow again. If it's all good, we can proceed with our analysis. </p>
<hr />
<h3>Appendix: List of all WARNINGs and ERRORs emitted by ValidateSamFile</h3>
<p>We are currently in the process of updating the Picard website to include the following two tables, describing <code>WARNING</code> (Table I) and <code>ERROR</code> (Table II) cases. Until that's done, you can find them here.</p>
<div>
<table class="table table-striped">
<tr><th><b>Table I</b></th><th></th></tr>
<tr><th><b>WARNING</b></th><th><b>Description</b></th></tr>
<tr><td><em>Header Issues</em></td><td> </td></tr>  
<tr><td>INVALID_DATE_STRING</td><td>Date string is not ISO-8601</td></tr>
<tr><td>INVALID_QUALITY_FORMAT  </td><td>   Quality encodings out of range; appear to be Solexa or Illumina when Phred expected.  Avoid exception being thrown as a result of no qualities being read.  </td></tr>
<tr><td><em>General Alignment Record Issues</em></td><td> </td></tr>
<tr><td>ADJACENT_INDEL_IN_CIGAR </td><td>   CIGAR string contains an insertion (I) followed by deletion (D), or vice versa  </td></tr>
<tr><td>RECORD_MISSING_READ_GROUP   </td><td>   A SAMRecord is found with no read group id  </td></tr>
<tr><td><em>Mate Pair Issues</em></td><td> </td></tr>
<tr><td>PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND   </td><td>   Pair flag set but not marked as first or second of pair </td></tr>
<tr><td><em>Optional Alignment Tag Issues</em></td><td> </td></tr>
<tr><td>MISSING_TAG_NM  </td><td>   The NM tag (nucleotide differences) is missing  </td></tr>
<tr><td>E2_BASE_EQUALS_PRIMARY_BASE </td><td>   Secondary base calls should not be the same as primary, unless one or the other is N    </td></tr>
<tr><td><em>General File, Index or Sequence Dictionary Issues</em></td><td> </td></tr>  
<tr><td>BAM_FILE_MISSING_TERMINATOR_BLOCK   </td><td>   BAM appears to be healthy, but is an older file so doesn't have terminator block    </td></tr>
</table>
</div>
<div>
<table class="table table-striped">
<tr><th><b>Table II</b></th><th></th></tr>
<tr><th><b> ERROR</b></th><th><b>   Description </b></th></tr>
<tr><td><em>Header Issues</em></td><td> </td></tr>  
<tr><td>    DUPLICATE_PROGRAM_GROUP_ID  </td><td>   Same program group id appears more than once    </td></tr>
<tr><td>    DUPLICATE_READ_GROUP_ID </td><td>   Same read group id appears more than once   </td></tr>
<tr><td>    HEADER_RECORD_MISSING_REQUIRED_TAG  </td><td>   Header tag missing in header line   </td></tr>
<tr><td>    HEADER_TAG_MULTIPLY_DEFINED </td><td>   Header tag appears more than once in header line with different value   </td></tr>
<tr><td>    INVALID_PLATFORM_VALUE  </td><td>   The read group has an invalid value set for its PL field    </td></tr>
<tr><td>    INVALID_VERSION_NUMBER  </td><td>   Does not match any of the acceptable versions   </td></tr>
<tr><td>    MISSING_HEADER  </td><td>   The SAM/BAM file is missing the header  </td></tr>
<tr><td>    MISSING_PLATFORM_VALUE  </td><td>   The read group is missing its PL (platform unit) field  </td></tr>
<tr><td>    MISSING_READ_GROUP  </td><td>   The header is missing read group information    </td></tr>
<tr><td>    MISSING_SEQUENCE_DICTIONARY </td><td>   There is no sequence dictionary in the header   </td></tr>
<tr><td>    MISSING_VERSION_NUMBER  </td><td>   Header has no version number    </td></tr>
<tr><td>    POORLY_FORMATTED_HEADER_TAG </td><td>   Header tag does not have colon  </td></tr>
<tr><td>    READ_GROUP_NOT_FOUND    </td><td>   A read group ID on a SAMRecord is not found in the header   </td></tr>
<tr><td>    UNRECOGNIZED_HEADER_TYPE    </td><td>   Header record is not one of the standard types  </td></tr>
<tr><td><em>General Alignment Record Issues</em></td><td> </td></tr>
<tr><td>    CIGAR_MAPS_OFF_REFERENCE    </td><td>   Bases corresponding to M operator in CIGAR extend beyond reference  </td></tr>
<tr><td>    INVALID_ALIGNMENT_START </td><td>   Alignment start position is incorrect   </td></tr>
<tr><td>    INVALID_CIGAR   </td><td>   CIGAR string error for either read or mate  </td></tr>
<tr><td>    INVALID_FLAG_FIRST_OF_PAIR  </td><td>   First of pair flag set for unpaired read    </td></tr>
<tr><td>    INVALID_FLAG_SECOND_OF_PAIR </td><td>   Second of pair flag set for unpaired read   </td></tr>
<tr><td>    INVALID_FLAG_PROPER_PAIR    </td><td>   Proper pair flag set for unpaired read  </td></tr>
<tr><td>    INVALID_FLAG_MATE_NEG_STRAND    </td><td>   Mate negative strand flag set for unpaired read </td></tr>
<tr><td>    INVALID_FLAG_NOT_PRIM_ALIGNMENT     </td><td>   Not primary alignment flag set for unmapped read    </td></tr>
<tr><td>    INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT    </td><td>   Supplementary alignment flag set for unmapped read  </td></tr>
<tr><td>    INVALID_FLAG_READ_UNMAPPED  </td><td>   Mapped read flat not set for mapped read    </td></tr>
<tr><td>    INVALID_INSERT_SIZE </td><td>   Inferred insert size is out of range    </td></tr>
<tr><td>    INVALID_MAPPING_QUALITY </td><td>   Mapping quality set for unmapped read or is >= 256  </td></tr>
<tr><td>    INVALID_PREDICTED_MEDIAN_INSERT_SIZE    </td><td>   PI tag value is not numeric </td></tr>
<tr><td>    MISMATCH_READ_LENGTH_AND_QUALS_LENGTH   </td><td>   Length of sequence string and length of base quality string do not match </td></tr>
<tr><td>    TAG_VALUE_TOO_LARGE </td><td>   Unsigned integer tag value is deprecated in BAM.  Template length   </td></tr>
<tr><td><em>Mate Pair Issues</em></td><td> </td></tr>
<tr><td>    INVALID_FLAG_MATE_UNMAPPED  </td><td>   Mate unmapped flag is incorrectly set   </td></tr>
<tr><td>    MATE_NOT_FOUND  </td><td>   Read is marked as paired, but its pair was not found    </td></tr>
<tr><td>    MATE_CIGAR_STRING_INVALID_PRESENCE  </td><td>   A cigar string for a read whose mate is NOT mapped  </td></tr>
<tr><td>    MATE_FIELD_MISMATCH </td><td>   Read alignment fields do not match its mate </td></tr>
<tr><td>    MATES_ARE_SAME_END  </td><td>   Both mates of a pair are marked either as first or second mates </td></tr>  
<tr><td>    MISMATCH_FLAG_MATE_UNMAPPED </td><td>   Mate unmapped flag does not match read unmapped flag of mate    </td></tr>
<tr><td>    MISMATCH_FLAG_MATE_NEG_STRAND   </td><td>   Mate negative strand flag does not match read strand flag   </td></tr>
<tr><td>    MISMATCH_MATE_ALIGNMENT_START   </td><td>   Mate alignment does not match alignment start of mate   </td></tr>
<tr><td>    MISMATCH_MATE_CIGAR_STRING  </td><td>   The mate cigar tag does not match its mate's cigar string   </td></tr>
<tr><td>    MISMATCH_MATE_REF_INDEX </td><td>   Mate reference index (MRNM) does not match reference index of mate  </td></tr>
<tr><td><em>Optional Alignment Tag Issues</em></td><td> </td></tr>
<tr><td>    INVALID_MATE_REF_INDEX  </td><td>   Mate reference index (MRNM) set for unpaired read   </td></tr>  
<tr><td>    INVALID_TAG_NM  </td><td>   The NM tag (nucleotide differences) is incorrect    </td></tr>
<tr><td>    MISMATCH_READ_LENGTH_AND_E2_LENGTH  </td><td>   Lengths of secondary base calls tag values and read should match    </td></tr>
<tr><td>    MISMATCH_READ_LENGTH_AND_U2_LENGTH  </td><td>   Secondary base quals tag values should match read length    </td></tr>
<tr><td>    EMPTY_READ  </td><td>   Indicates that a read corresponding to the first strand has a length of zero and/or lacks flow signal intensities (FZ)      </td></tr>
<tr><td>    INVALID_INDEXING_BIN    </td><td>   Indexing bin set on SAMRecord does not agree with computed value    </td></tr>
<tr><td><em>General File, Index or Sequence Dictionary Issues</em></td><td> </td></tr>  
<tr><td>    INVALID_INDEX_FILE_POINTER  </td><td>   Invalid virtualFilePointer in index     </td></tr>
<tr><td>    INVALID_REFERENCE_INDEX </td><td>   Reference index not found in sequence dictionary    </td></tr>
<tr><td>    RECORD_OUT_OF_ORDER     </td><td>   The record is out of order  </td></tr>
<tr><td>    TRUNCATED_FILE  </td><td>   BAM file does not have terminator block </td></tr>
</table>
</div>