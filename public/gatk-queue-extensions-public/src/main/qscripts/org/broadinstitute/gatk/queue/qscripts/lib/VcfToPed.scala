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
import org.broadinstitute.gatk.utils.commandline.Input
import org.broadinstitute.gatk.queue.library.ipf.vcf.VCFExtractIntervals
import org.broadinstitute.gatk.utils.text.XReadLines
import collection.JavaConversions._
import java.io._
import org.broadinstitute.gatk.queue.extensions.gatk.{SelectVariants, VariantsToPed}

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/31/12
 * Time: 10:46 PM
 * To change this template use File | Settings | File Templates.
 */

class VcfToPed extends QScript {

  @Input(shortName = "V", fullName="Variants", required=true,doc="VCF to convert to ped")
  var variants : File = _

  @Output(shortName = "B", fullName="Bed",required=true,doc="Name of the ped output file (fam and bim will use the root of this file)")
  var bed : File = _

  @Input(shortName = "M", fullName="Meta",required=true,doc="The sample metadata file, can be a .fam or [NAME]\\tkey1=val1;key2=val2")
  var meta : File = _

  @Input(shortName = "Int", fullName="Intervals",required=false,doc="Intervals. If not specified script will produce them and exit.")
  var intervals : File = _

  @Argument(shortName="R",fullName="Ref",required=false,doc="Reference file")
  var ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName="D",fullName="dbsnp",required=false,doc="dbsnp file")
  var dbsnp : File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf")

  @Argument(shortName="sf",fullName="sampleFile",required=false,doc="sample file")
  var samFile : File = _

  val tmpdir : File = System.getProperty("java.io.tmpdir")

  def script = {
    if ( intervals == null ) {
      val ivals : File = swapExt(variants,".vcf",".intervals.list")
      val extract : VCFExtractIntervals = new VCFExtractIntervals(variants,ivals,false)
      add(extract)
    } else {
      val IS_GZ : Boolean = variants.getName.endsWith(".vcf.gz")
      var iXRL = new XReadLines(intervals)
      var chunk = 1;
      var subListFile : File = null
      if ( IS_GZ )
        subListFile = swapExt(tmpdir,variants,".vcf.gz",".chunk%d.list".format(chunk))
      else
        subListFile = swapExt(tmpdir,variants,".vcf",".chunk%d.list".format(chunk))
      var subList = new PrintStream(subListFile)
      var nL = 0;
      var bedOuts : List[File] = Nil;
      var bimOuts : List[File] = Nil
      var lastFam : File = null;
      while ( iXRL.hasNext ) {
        subList.printf("%s%n",iXRL.next())
        nL = nL + 1
        if ( nL > 10000 ) {
          val toPed : VariantsToPed = new VariantsToPed
          toPed.memoryLimit = 2
          toPed.reference_sequence = ref
          toPed.intervals :+= subListFile
          toPed.dbsnp = dbsnp
          if ( samFile != null ) {
            val base : String = bed.getName.stripSuffix(".bed")+"_%d".format(chunk)
            val extract : SelectVariants = new SelectVariants
            extract.reference_sequence = ref
            extract.memoryLimit = 2
            extract.intervals :+= subListFile
            extract.variant = variants
            extract.out = new File(tmpdir,base+"_extract%d.vcf".format(chunk))
            extract.sample_file :+= samFile
            add(extract)
            toPed.variant = extract.out
          } else {
            toPed.variant = variants
          }
          toPed.metaData = meta
          val base : String = bed.getName.stripSuffix(".bed")+"_%d".format(chunk)
          val tBed = new File(tmpdir,base+".bed")
          val bim = new File(tmpdir,base+".bim")
          val fam = new File(tmpdir,base+".fam")
          toPed.bed = tBed
          toPed.bim = bim
          toPed.fam = fam
          add(toPed)
          subList.close()
          chunk = chunk + 1
          if ( IS_GZ  )
            subListFile = swapExt(tmpdir,variants,".vcf.gz",".chunk%d.list".format(chunk))
          else
            subListFile = swapExt(tmpdir,variants,".vcf",".chunk%d.list".format(chunk))
          subList = new PrintStream(subListFile)
          bedOuts :+= tBed
          bimOuts :+= bim
          lastFam = fam
          nL = 0;
        }
      }

      if ( nL > 0 ) {
        val toPed : VariantsToPed = new VariantsToPed
        toPed.reference_sequence = ref
        toPed.intervals :+= new File(subListFile)
        toPed.dbsnp = dbsnp
        if ( samFile != null ) {
          val base : String = bed.getName.stripSuffix(".bed")+"_%d".format(chunk)
          val extract : SelectVariants = new SelectVariants
          extract.reference_sequence = ref
          extract.memoryLimit = 2
          extract.intervals :+= subListFile
          extract.variant = variants
          extract.out = new File(tmpdir,base+"_extract%d.vcf".format(chunk))
          extract.sample_file :+= samFile
          add(extract)
          toPed.variant = extract.out
        } else {
          toPed.variant = variants
        }
        toPed.metaData = meta
        toPed.memoryLimit = 2
        val base : String = bed.getName.stripSuffix(".bed")+"_%d".format(chunk)
        val tBed = new File(tmpdir,base+".bed")
        val bim = new File(tmpdir,base+".bim")
        val fam = new File(tmpdir,base+".fam")
        toPed.bed = tBed
        toPed.bim = bim
        toPed.fam = fam
        lastFam = fam
        add(toPed)
        subList.close()
        bedOuts :+= tBed
        bimOuts :+= bim
      }

      var gatherUP = new MyPedGather
      gatherUP.binPed = bedOuts
      gatherUP.bim = bimOuts
      gatherUP.outPed = bed
      gatherUP.outBim = swapExt(bed,".bed",".bim")

      add(gatherUP)

      class copyFam extends InProcessFunction {
        @Input(doc="fam") var inFam = lastFam
        @Output(doc="fam") var outFam = swapExt(bed,".bed",".fam")

        def run = {
          var stream = new PrintStream(outFam)
          asScalaIterator(new XReadLines(inFam)).foreach( u => {
            stream.printf("%s%n",u)
          })
          stream.close()
        }
      }

      add(new copyFam)
    }

  }

   class MyPedGather extends InProcessFunction {
    @Input(doc="Peds to be merged") var binPed: List[File] = Nil
    @Input(doc="Bims to be merged") var bim : List[File] = Nil
    @Output(doc="The final Ped to write to") var outPed : File = _
    @Output(doc="The final bim to write to") var outBim : File = _

    def run : Unit = {
      var stream : PrintStream = new PrintStream(outPed)
      stream.write((List[Byte](0x6c.toByte,0x1b.toByte,0x1.toByte)).toArray)
      binPed.map(u => new FileInputStream(u) ).foreach( u => {
        u.skip(3)
        var b = -1
        do {
          b = u.read()
          stream.write(b.toByte)
        } while ( b != -1 )
      })
      stream.close()

     stream = new PrintStream(outBim)
     bim.map(u => new XReadLines(u)).foreach( u => {
       asScalaIterator(u).foreach( x => {
         stream.printf("%s%n",x)
       })
     })

      stream.close()
    }
  }

}