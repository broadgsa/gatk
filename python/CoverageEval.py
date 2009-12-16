#!/usr/bin/env python

import sys, itertools, FlatFileTable
from enum import Enum

def subset_list_by_indices(indices, list):
    subset = []
    for index in indices:
        subset.append(list[index])
    return subset

def chunk_generator(record_gen, key_fields):
    """Input:
  record_gen: generator that produces dictionaries (records in database speak)
  key_fields: keys in each dictionary used to determine chunk membership
Output:
  locus_chunk: list of consecutive lines that have the same key_fields"""
  
    locus_chunk = []
    last_key = ""
    first_record = True
    for record in record_gen:
        key = [record[f] for f in key_fields]
        if key == last_key or first_record:
            locus_chunk.append(record)
            first_record = False
        else:
            if locus_chunk != []:
                yield locus_chunk
                locus_chunk = [record]
        last_key = key
    yield locus_chunk


class call_stats:
    def __init__(self, call_type, coverage): #, acc_conf_calls, conf_call_rate, cum_corr_calls, cum_calls, coverage):
        self.call_type = call_type
        self.coverage = coverage
        self.calls = 0
        self.conf_ref_calls = 0
        self.conf_het_calls = 0
        self.conf_hom_calls = 0
        self.conf_var_calls = 0
        self.conf_genotype_calls = 0
        self.conf_refvar_calls = 0
        self.correct_genotype = 0
        self.correct_refvar = 0


    def add_stat(self, calls, conf_ref_calls, conf_het_calls, conf_hom_calls, conf_var_calls, conf_genotype_calls, conf_refvar_calls, correct_genotype, correct_refvar):
        #print self, calls, conf_ref_calls, conf_het_calls, conf_hom_calls, conf_var_calls, conf_genotype_calls, conf_refvar_calls, correct_genotype, correct_refvar
        self.calls += calls
        self.conf_ref_calls += conf_ref_calls
        self.conf_het_calls += conf_het_calls
        self.conf_hom_calls += conf_hom_calls
        self.conf_var_calls += conf_var_calls
        self.conf_genotype_calls += conf_genotype_calls
        self.conf_refvar_calls += conf_refvar_calls
        self.correct_genotype += correct_genotype
        self.correct_refvar += correct_refvar

    @staticmethod
    def calc_discovery_stats(chunk):
        conf_calls = 0
        correct_genotype = 0
        calls = 0
        for record in chunk:
            calls += 1
            if float(record["BtrLod"]) >= 5 or call_type.from_record_2_state(record["ReferenceBase"][0],record["BestGenotype"]) == call_type.call_types_2_state.Reference and float(record["BtnbLod"]) >= 5:
                conf_calls += 1
                if call_type.discovery_call_correct(record):
                    correct_genotype += 1

        return correct_genotype, conf_calls, calls

    @staticmethod
    def calc_genotyping_stats(chunk):
        conf_calls = 0
        correct_genotype = 0
        calls = 0
        for record in chunk:
            calls += 1
            if float(record["BtnbLod"]) >= 5:
                conf_calls += 1
                if call_type.genotyping_call_correct(record):
                    correct_genotype += 1
                    
        return correct_genotype, conf_calls, calls

    @staticmethod
    def stats_header():
        return "TrueGenotype,Coverage,AccuracyConfidentGenotypingCalls,ConfidentGenotypingCallRate,AccuracyConfidentDiscoveryCalls,ConfidentDiscoveryCallRate,Calls,ConfRefCalls,ConfHetCalls,ConfHomCalls,ConfGenotypeCalls,CorrectGenotypes,ConfVarCalls,ConfDiscoveryCalls,CorrectDiscovery"
        
    def __str__(self):
        return ",".join(map(str, (self.calls, self.conf_ref_calls, self.conf_het_calls, self.conf_hom_calls, self.conf_genotype_calls, self.correct_genotype, self.conf_var_calls, self.conf_refvar_calls, self.correct_refvar)))

    def stats_str(self):
        return "%s,%d,%.5f,%.5f,%.5f,%.5f,%s" % (self.call_type, self.coverage, float(self.correct_genotype) / max(self.conf_genotype_calls,1), float(self.conf_genotype_calls) / max(self.calls,1), float(self.correct_refvar) / max(self.conf_refvar_calls,1), float(self.conf_refvar_calls) / max(self.calls,1), self.__str__())

class call_type:
    """Class that returns an Enum with the type of call provided by a record"""
    call_types_3_state = Enum("HomozygousSNP","HeterozygousSNP","HomozygousReference")
    call_types_3_state_short = Enum("Hom","Het","Ref")    
    call_types_2_state = Enum("Variant","Reference")
    
    @staticmethod
    def from_record_3_state(ref, genotype): # record):
        """Given reference base as string, determine whether called genotype is homref, het, homvar"""
        #ref = record["ReferenceBase"][0]
        #genotype = record["HapmapChipGenotype"]
        return call_type.call_types_3_state[genotype.count(ref)]

    @staticmethod
    def from_record_2_state(ref, genotype):
        """Given reference base as string, determine whether called genotype is ref or var"""
        #ref = record["ReferenceBase"][0]
        #genotype = record["HapmapChipGenotype"]
        return  call_type.call_types_2_state[0] if genotype.count(ref) < 2 else call_type.call_types_2_state[1]

    @staticmethod
    def genotyping_call_correct(record):
        return record["HapmapChipGenotype"] == record["BestGenotype"]
 

    @staticmethod
    def discovery_call_correct(record):
        return call_type.from_record_2_state(record["ReferenceBase"][0], record["HapmapChipGenotype"]) == call_type.from_record_2_state(record["ReferenceBase"][0], record["BestGenotype"])

def print_record_debug(rec):
    print "TEST Genotyping:", call_type.genotyping_call_correct(rec)
    print "TEST Discovery:", call_type.discovery_call_correct(rec)
    print "DIFF:", call_type.genotyping_call_correct(rec) != call_type.discovery_call_correct(rec)
    print "\n".join(["   %20s => '%s'" % (k,v) for k,v in sorted(rec.items())])
    print call_type.from_record_2_state(rec["ReferenceBase"][0],rec["HapmapChipGenotype"])
    print call_type.from_record_2_state(rec["ReferenceBase"][0],rec["BestGenotype"])
    print

def aggregate_stats(filename, max_loci, table_filename, debug):
    aggregate = dict()
    fout = open(table_filename,"w")
    fout.write(call_stats.stats_header()+"\n")

    locus_gen = chunk_generator(FlatFileTable.record_generator(filename, None), ("Sequence","Position"))
    for index, locus_chunk in enumerate(locus_gen):
        if index >= max_loci:
            break
        if (index % 1000) == 0:
            sys.stderr.write( str(index)+" loci processed, at: "+locus_chunk[0]["Sequence"]+":"+locus_chunk[0]["Position"]+"\n")

        covs = dict()
        coverage_chunk_gen = chunk_generator(locus_chunk, ("DownsampledCoverage", "Sequence", "Position"))
        for cov_chunk in coverage_chunk_gen:
            first_record = cov_chunk[0]
            #stat = call_stats.calc_stats(cov_chunk)
            
            for record in cov_chunk:
                hapmap_genotyping_call_type = call_type.from_record_3_state(record["ReferenceBase"][0],record["HapmapChipGenotype"])
                key = hapmap_genotyping_call_type, int(record["DownsampledCoverage"])
                #key = call_type.from_record_3_state(record)#, int(record["DownsampledCoverage"])
                record["DownsampledCoverage"] = int(record["DownsampledCoverage"])
                record["HapmapChipCallType"] = hapmap_genotyping_call_type
                value = record

                correct_genotype, conf_genotype_calls, genotype_calls = call_stats.calc_genotyping_stats([record])
                correct_refvar, conf_refvar_calls, refvar_calls = call_stats.calc_discovery_stats([record])
                assert(genotype_calls == refvar_calls)

                conf_ref_calls = 0
                conf_het_calls = 0
                conf_hom_calls = 0
                best_genotyping_call_type = call_type.from_record_3_state(record["ReferenceBase"][0],record["BestGenotype"])
                if conf_genotype_calls: 
                    if best_genotyping_call_type.index == 0: conf_hom_calls = 1
                    if best_genotyping_call_type.index == 1: conf_het_calls = 1
                    if best_genotyping_call_type.index == 2: conf_ref_calls = 1

                conf_var_calls = 0
                if conf_refvar_calls:
                    this_variant_call_type = call_type.from_record_2_state(record["ReferenceBase"][0],record["BestGenotype"])
                    conf_var_calls = 1 if this_variant_call_type.index == 0 else 0
                
                aggregate.setdefault(key, call_stats(*key))
                #print ",".join(map(str,(genotype_calls, conf_ref_calls, conf_het_calls, conf_hom_calls, conf_var_calls, conf_genotype_calls, conf_refvar_calls, correct_genotype, correct_refvar)))
                aggregate[key].add_stat(genotype_calls, conf_ref_calls, conf_het_calls, conf_hom_calls, conf_var_calls, conf_genotype_calls, conf_refvar_calls, correct_genotype, correct_refvar)

                if debug:# and conf_refvar_calls:
                    print "KEYS:",key
                    print_record_debug(record)


    for key, records in sorted(aggregate.items()):
        fout.write(records.stats_str()+"\n")

    fout.close()
    #return aggregate

class weighted_avg:

    def __init__(self):
        self.sum = 0.0
        self.count = 0

    def add(self, value, counts):
        self.sum += value*counts
        self.count += counts

    def return_avg(self):
        return float(self.sum) / max(self.count,1)

def stats_from_hist(options, depth_hist_filename, stats_filename, variant_eval_dir, depth_multiplier=1.0):

    #hist_zero = {"CallType" : ,"Coverage","AccuracyConfidentCalls","ConfidentCallRate","CumCorrectCalls","CumCalls","CumgCorrectFraction"}
    #prob_genotype = [1e-5, 1e-3, 1-1e-3]
    #prob_genotype = [0.37, 0.62, .0]
    #prob_genotype = [0.203, 0.304, .491]
    #prob_genotype = [0.216, 0.302, .481]
    prob_genotype = [0.205, 0.306, 0.491] # Based on CEU NA12878 actual hapmap chip calls
    #prob_genotype = [0.213, 0.313, 0.474] # Based on YRB NA19240 actual hapmap chip calls
    theta = 1.0/1850 # expected heterozygosity
    prob_genotype = [0.5 * theta, 1.0 * theta, 1 - 1.5*theta] # Based on CEU NA12878 actual hapmap chip calls
    
    
    hist = []
    hist_gen = FlatFileTable.record_generator(depth_hist_filename, sep=" ", skip_until_regex_line = "^depth count freq")
    for index, record in enumerate(hist_gen):
        assert int(record["depth"]) == index
        hist.append(int(record["count"]))
    
    # If upsampling is not done in the CoverageEval GATK module, the number of observations of reads 
    # with high depth of coverage can be very low giving 
        
    stats_dict = dict()
    stats_gen = FlatFileTable.record_generator(stats_filename, sep=",")
    for record in stats_gen:
        key1 = int(record["Coverage"])
        key2 = record["TrueGenotype"]
        #print key2
        stats_dict.setdefault(key1, dict()) # create a nested dict if it doesn't exist
        stats_dict[key1][key2] = record # create an entry for these keys

        #highest_homref_calls = 0
        #if record["TrueGenotype"] == "HomozygousReference":
        #    calls = int(record["Calls"])
        #    if calls > highest_homref_calls:
        #        highest_homref_calls = calls
        
    acc = dict()
    call_rate = dict()
    conf_calls = dict()
    
    start = 1
    end = 1000
    max_usable_depth = 40 # Depth of coverage beyond which stats are not sampled enough and we take the stat at this depth instead
    for depth, depth_count in enumerate(hist[start:end+1],start): # For Cd = depth count
        #print "DEPTH: "+str(depth)
        try:
            depth = max(int(float(depth*depth_multiplier)),1)
            if depth > max_usable_depth: # Ensure that high entries with bad stats use a good stat from a depth that we now is well sampled
                depth = max_usable_depth
            depth_entries = stats_dict[depth]
        except KeyError:
            print "Stopped on depth",depth
            break
        if True:
            for true_genotype, stat in depth_entries.items(): # For t (SNP type) = true_genotype
                #print "TRUE_GENOTYPE: "+str(true_genotype)
                for genotype in call_type.call_types_3_state:
                    conf_calls.setdefault(genotype, 0.0)
                    prob_conf_x_call = float(stat["Conf"+str(call_type.call_types_3_state_short[genotype.index])+"Calls"])/float(stat["Calls"])
                    conf_calls[genotype] += depth_count * prob_conf_x_call * prob_genotype[genotype.index]
                    #if genotype.index == 1:
                    #    print "%.5f " % prob_conf_x_call, depth, depth_count, conf_calls[genotype], int(stat["Conf"+str(call_type.call_types_3_state_short[genotype.index])+"Calls"]), int(stat["Calls"])

                acc.setdefault(true_genotype,weighted_avg())
                call_rate.setdefault(true_genotype,weighted_avg())
                acc[true_genotype].add(float(stat["AccuracyConfidentGenotypingCalls"]), depth_count)
                call_rate[true_genotype].add(float(stat["ConfidentGenotypingCallRate"]), depth_count)

    import numpy
    for genotype in call_type.call_types_3_state:
        print "%19s accuracy : %.3f" % (str(genotype), acc[str(genotype)].return_avg())
        print "%19s call rate: %.3f" % (str(genotype), call_rate[str(genotype)].return_avg())

    print "\nExpected performance given perfect accuracy and call rate:"
    print "%19s        %7s %7s %7s" % ("", "Actual", "Perfect", "Diff.")
    total_hist_sites = numpy.sum(hist)
    total_predicted = 0
    for genotype in call_type.call_types_3_state:
        predicted = conf_calls[genotype]
        total_predicted += predicted
        perfect = prob_genotype[genotype.index]*total_hist_sites
        diff = perfect - predicted
        percent_of_possible = predicted / perfect * 100
        print "%19s calls: %8.0f %8.0f %8.0f %8.1f" % (genotype, predicted, perfect, diff, percent_of_possible)
        #repl_string += "%s %.0f\n" % (genotype, predicted)
    print "              Total calls: %7d" % total_predicted
    
    #stats_gen = FlatFileTable.record_generator(stats_filename, sep=",")
    #for chunk in chunk_generator(stats_gen, key_fields=("True_Genotype")):
    
    print "\nCoverage histogram mean: %.2f" % numpy.average(range(len(hist)), weights=hist)
    
    # STEM AND LEAF PLOT
    
    # If VariantEval directory given, compare with those results
    if options.variant_eval_file != None:
        vareval_file = open(options.variant_eval_file)
        num_hets = None; num_homs = None
        for line in vareval_file:
            if "UNKNOWN_CALLED_VAR_HET_NO_SITES" in line: num_hets = line.rstrip().split()[2]
            if "UNKNOWN_CALLED_VAR_HOM_NO_SITES" in line: num_homs = line.rstrip().split()[2]
            
    if options.oneline_stats != None:
        oneline_stats = open(options.oneline_stats, "w")
        pred_hom = conf_calls[call_type.call_types_3_state[0]]
        pred_het = conf_calls[call_type.call_types_3_state[1]]
        print >>oneline_stats, "%s %.0f %s %.0f %s" % (options.oneline_stats, pred_hom, num_homs, pred_het, num_hets)  

def usage(parser):
    parser.print_usage()
    sys.exit()

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c", "--genotype_call_file", help="GELI file to use in generating coverage stat table", dest="genotype_filename", )
    parser.add_option("-s", "--stats_file", help="file to containing empirical genotyping stats", dest="stats_filename", default=None)
    parser.add_option("-g", "--histogram_file", help="file containing counts of each depth of coverage", dest="hist_filename")
    parser.add_option("-m", "--max_loci", help="maximum number of loci to parse (for debugging)", default=sys.maxint, dest="max_loci", type="int")
    parser.add_option("-e", "--evaluate", help="evaluate genotypes; requires a stats file and a histogram file", default=False, dest="evaluate_genotypes", action="store_true")
    parser.add_option("-d", "--debug", help="provide debugging output", default=False, dest="debug", action="store_true")
    parser.add_option("-p", "--depth_multiplier", help="multiply all depths in histogram by this value; for \"downsampling\" depth", default=1.0, dest="depth_multiplier", type=float)
    parser.add_option("-v", "--variant_eval_file", help="file with output of VariantEval genotype concordance to compare this prediction to", default=None, dest="variant_eval_file")
    parser.add_option("-o", "--oneline_stats_file", help="output single, tabular line of stats to this file", default=None, dest="oneline_stats")
    
    (options, args) = parser.parse_args()

    if not options.evaluate_genotypes:
        print "Creating performance tables from genotypes file"
        #if options.stats_filename == None:
        #    sys.exit("Must provide -s stats fliname option")
        if options.genotype_filename == None:
            sys.exit("Must provide -c genotype call filename option")
        stats_filename = options.stats_filename if options.stats_filename != None else options.genotype_filename+".stats"
        aggregate_stats(options.genotype_filename, options.max_loci, options.stats_filename, options.debug)
    else:
        print "Evaluating genotypes"
        print "Depth multiplier:",options.depth_multiplier
        if options.hist_filename == None:
            sys.exit("Must provide -g histogram filename option")
        if options.stats_filename == None:
            sys.exit("Must provide -s stats fliname option")
        stats_from_hist(options, options.hist_filename, options.stats_filename, options.variant_eval_file, options.depth_multiplier)
        


