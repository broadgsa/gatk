/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.qscripts.lib

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.library.ipf.vcf.VCFExtractIntervals
import scala.collection.JavaConversions._
import org.broadinstitute.gatk.utils.text.XReadLines
import java.io.PrintStream
import org.broadinstitute.gatk.queue.extensions.gatk.SelectVariants

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/2/12
 * Time: 12:13 PM
 * To change this template use File | Settings | File Templates.
 */

class ChunkVCF extends QScript {

  @Input(shortName="V",fullName="VCF",doc="The VCF you want to chunk",required=true)
  var inVCF : File = _

  @Input(shortName="N",fullName="numEntriesInChunk",doc="The number of variants per chunk",required=true)
  var numEntries : Int = _

  @Input(shortName="I",fullName="Intervals",doc="The SNP interval list to chunk. If not provided, one will be created for you to provide in a second run.",required=false)
  var intervals : File = _

  @Input(fullName="preserveChromosomes",doc="Restrict chunks to one chromosome (smaller chunk at end of chromosome)",required=false)
  var preserve : Boolean = false

  @Input(fullName="reference",doc="The reference file",required=false)
  var ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Input(fullName="samples",doc="A file of sample IDs to condense VCF file to",required=false)
  var extractSamples : File = _

  val tmpdir : File = System.getProperty("java.io.tmpdir")

  def script = {
    if ( intervals == null ) {
      // create an interval list from the VCF
      val ivals : File = swapExt(inVCF,".vcf",".intervals.list")
      val extract : VCFExtractIntervals = new VCFExtractIntervals(inVCF,ivals,false)
      add(extract)
    } else {
      var chunkNum = 1
      var numLinesInChunk = 0
      var chromosome : String = asScalaIterator(new XReadLines(intervals)).next().split(":")(0)
      var chunkFile : File = new File(tmpdir,"ChunkVCF.chunk%d.intervals.list".format(chunkNum))
      var chunkWriter = new PrintStream(chunkFile)
      asScalaIterator(new XReadLines(intervals)).foreach( int => {
        // check new chromosome or full chunk
        if ( ( preserve && ! int.split(":")(0).equals(chromosome) ) || numLinesInChunk > numEntries ) {
          chunkWriter.close()
          val chunkSelect : SelectVariants = new SelectVariants
          chunkSelect.variant = inVCF
          chunkSelect.reference_sequence = ref
          chunkSelect.memoryLimit = 2
          chunkSelect.intervals :+= chunkFile
          if ( extractSamples != null )
            chunkSelect.sample_file :+= extractSamples
          chunkSelect.out = swapExt(inVCF,".vcf",".chunk%d.vcf".format(chunkNum))
          add(chunkSelect)
          chunkNum += 1
          numLinesInChunk = 0
          chromosome = int.split(":")(0)
          chunkFile = new File(tmpdir,"ChunkVCF.chunk%d.intervals.list".format(chunkNum))
          chunkWriter = new PrintStream(chunkFile)
        }
        chunkWriter.printf("%s%n",int)
        numLinesInChunk += 1
      })
      // last chunk
      if ( numLinesInChunk > 0 ) {
        // some work to do
        val chunkSelect : SelectVariants = new SelectVariants
        chunkSelect.variant = inVCF
        chunkSelect.reference_sequence = ref
        chunkSelect.memoryLimit = 2
        chunkSelect.intervals :+= chunkFile
        chunkWriter.close()
        if ( extractSamples != null )
          chunkSelect.sample_file :+= extractSamples
        chunkSelect.out = swapExt(inVCF,".vcf",".chunk%d.vcf".format(chunkNum))
        add(chunkSelect)
      }
    }
  }
}