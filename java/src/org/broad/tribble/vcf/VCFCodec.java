package org.broad.tribble.vcf;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.util.LineReader;
import org.broadinstitute.sting.gatk.refdata.features.vcf4.VCF4Codec;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class VCFCodec
 *
 * The codec for VCF, which relies on VCFReaderUtils to do most of the processing
 */
public class VCFCodec extends VCF4Codec {}
