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
  line_gen: generator that produces dictionaries
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
    def __init__(self, acc_conf_calls, conf_call_rate, cum_corr_calls, cum_calls, coverage):
        self.AccuracyConfidentCalls = acc_conf_calls
        self.ConfidentCallRate = conf_call_rate
        self.CumulativeConfidentCorrectCalls = cum_corr_calls
        self.CumulativeCalls = cum_calls
        self.Coverage = coverage
#    def stat_generator(chunk):

    @staticmethod
    def calc_discovery_stats(chunk):
        calls = 0
        conf_calls = 0
        correct_genotype = 0
        for record in chunk:
            if abs(float(record["BtrLod"])) >= 5:
                conf_calls += 1
                if call_type.discovery_call_correct(record):
                #if call_type.genotyping_call_correct(record):
                    correct_genotype += 1
                    
            calls += 1

        return correct_genotype, conf_calls, calls

    @staticmethod
    def calc_genotyping_stats(chunk):
        calls = 0
        conf_calls = 0
        correct_genotype = 0
        for record in chunk:
            if abs(float(record["BtnbLod"])) >= 5:
                conf_calls += 1
                if call_type.genotyping_call_correct(record):
                    correct_genotype += 1
                    
            calls += 1

        return correct_genotype, conf_calls, calls
        #return call_stats(float(correct_genotype) / max(conf_calls,1), float(conf_calls) / max(calls,1))

    def __str__(self):
        return "%d,%.5f,%.5f,%d,%d,%.5f" % (self.Coverage, self.AccuracyConfidentCalls, self.ConfidentCallRate, self.CumulativeConfidentCorrectCalls, self.CumulativeCalls, float(self.CumulativeConfidentCorrectCalls)/self.CumulativeCalls )

class call_type:
    """Class that returns an Enum with the type of call provided by a record"""
    call_types_3_state = Enum("HomozygousSNP","HeterozygousSNP","HomozygousReference")
    call_types_2_state = Enum("Variant","Reference")
    
    @staticmethod
    def from_record_3_state(record):
        """Given reference base as string, determine whether called genotype is homref, het, homvar"""
        ref = record["ReferenceBase"][0]
        genotype = record["HapmapChipGenotype"]
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


def aggregate_stats(filename, max_loci):
    aggregate = dict()

    locus_gen = chunk_generator(FlatFileTable.record_generator(filename, None), ("Sequence","Position"))
    #print "Fraction correct genotype\tCoverage sampled\tLocus\tReference base\tHapmap chip genotype (Max. coverage genotype call for reference calls)"
    for index, locus_chunk in enumerate(locus_gen):
        if index >= max_loci:
            break
        if (index % 1000) == 0:
            sys.stderr.write( str(index)+" loci processed, at: "+locus_chunk[0]["Sequence"]+":"+locus_chunk[0]["Position"]+"\n")

        covs = dict()
        coverage_chunk_gen = chunk_generator(locus_chunk, ("DownsampledCoverage", "Sequence", "Position"))
        for cov_chunk in coverage_chunk_gen:
            #record = cov_chunk[0]
            #stat = call_stats.calc_stats(cov_chunk)
            
            for record in cov_chunk:
                key = call_type.from_record_3_state(record), int(record["DownsampledCoverage"])
                #key = call_type.from_record_3_state(record)#, int(record["DownsampledCoverage"])
                record["DownsampledCoverage"] = int(record["DownsampledCoverage"])
                record["HapmapChipCallType"] = call_type.from_record_3_state(record)
                value = record
                if aggregate.has_key(key):
                    aggregate[key].append(value)
                else:
                    aggregate[key] = [value]

    #print "\t".join(map(str,("%.4f\t%.4f" % (stat.AccuracyConfidordentCalls, stat.ConfidentCallRate), record["DownsampledCoverage"], record["Sequence"]+":"+record["Position"],record["ReferenceBase"],record["HapmapChipGenotype"])))
    #print "\n".join(map(str,sorted(aggregate.items())))

    return aggregate

def create_coverage_stats_table(aggregate, table_filename, debug):
    fout = open(table_filename,"w")

    print >>fout, "CallType,Coverage,AccuracyConfidentCalls,ConfidentCallRate,CumCorrectCalls,CumCalls,CumCorrectFraction"

    cum_correct_calls = [0,0,0]
    cum_calls = [0,0,0]

    for key, records in sorted(aggregate.items()):
        if debug:
            print "KEYS:",key
            for rec in records:
                if True: #abs(float(rec["BtrLod"])) > 5:
                    print "TEST Genotyping:", call_type.genotyping_call_correct(rec)
                    print "TEST Discovery:", call_type.discovery_call_correct(rec)
                    print "DIFF:", call_type.genotyping_call_correct(rec) != call_type.discovery_call_correct(rec)
                    print "\n".join(["   %20s => '%s'" % (k,v) for k,v in sorted(rec.items())])
                    print call_type.from_record_2_state(rec["ReferenceBase"][0],rec["HapmapChipGenotype"])
                    print call_type.from_record_2_state(rec["ReferenceBase"][0],rec["BestGenotype"])
                    print
            print

        #print "\n".join(["%s => %s" % record.items() for record in records])
        if options.do_discovery:
            correct_genotype, conf_calls, calls = call_stats.calc_discovery_stats(records)
        else:
            correct_genotype, conf_calls, calls = call_stats.calc_genotyping_stats(records)
        this_call_type = call_type.from_record_3_state(records[0])
        cum_correct_calls[this_call_type.index] += correct_genotype
        cum_calls[this_call_type.index] += calls
        #yield call_stats(float(correct_genotype) / max(conf_calls,1), float(conf_calls) / max(calls,1))
        record = records[0]

        print >>fout, str(record["HapmapChipCallType"])+","+str(call_stats(float(correct_genotype) / max(conf_calls,1), float(conf_calls) / max(calls,1), cum_correct_calls[this_call_type.index], cum_calls[this_call_type.index], record["DownsampledCoverage"]))
        # record["HapmapChipCallType"])

class weighted_avg:

    def __init__(self):
        self.sum = 0.0
        self.count = 0

    def add(self, value, counts):
        self.sum += value*counts
        self.count += counts
        #print value, counts, self.sum, self.count

    def return_avg(self):
        return float(self.sum) / max(self.count,1)

def stats_from_hist(depth_hist_filename, stats_filename):

    #hist_zero = {"CallType" : ,"Coverage","AccuracyConfidentCalls","ConfidentCallRate","CumCorrectCalls","CumCalls","CumCorrectFraction"}

    hist = []
    hist_gen = FlatFileTable.record_generator(depth_hist_filename, sep=" ", skip_n_lines=2)
    for index, record in enumerate(hist_gen):
        assert int(record["depth"]) == index
        hist.append(int(record["count"]))
    
    stats_dict = dict()
    stats_gen = FlatFileTable.record_generator(stats_filename, sep=",")
    for record in stats_gen:
        key1 = int(record["Coverage"])
        key2 = record["CallType"]
        stats_dict.setdefault(key1, dict()) # create a nested dict if it doesn't exist
        stats_dict[key1][key2] = record # create an entry for these keys
        
    #print stats_dict

    acc = dict() #[weighted_avg()] * 3
    call_rate = dict() #[weighted_avg()] * 3
    
    start = 10
    end = 10
    for depth, depth_count in enumerate(hist[start:end+1],start):
        #print "DEPTH: "+str(depth)
        try:
            depth_entries = stats_dict[depth]
            for calltype, stat in depth_entries.items():
                acc.setdefault(calltype,weighted_avg())
                call_rate.setdefault(calltype,weighted_avg())
                acc[calltype].add(float(stat["AccuracyConfidentCalls"]), depth_count)
                call_rate[calltype].add(float(stat["ConfidentCallRate"]), depth_count)

        #    acc[calltype] = stat
            #ref = depth_entries["HomozygousReference"]
            #het = depth_entries["HeterozygousSNP"]
            #hom = depth_entries["HomozygousSNP"]
        except KeyError:
            break
        
        #acc.add(float(ref["AccuracyConfidentCalls"]), depth_count)
        #call_rate.add(float(ref["ConfidentCallRate"]), depth_count)

        #print(float(het["AccuracyConfidentCalls"]), depth_count)
        #print(float(het["ConfidentCallRate"]), depth_count)

    for calltype in ("HomozygousSNP","HeterozygousSNP","HomozygousReference"):
        print "%25s accuracy : %.3f" % (calltype, acc[calltype].return_avg())
        print "%25s call rate: %.3f" % (calltype, call_rate[calltype].return_avg())

        
def usage(parser):
    #print "Usage: CoverageEval.py geli_file OPTIONS"
    parser.print_usage()
    sys.exit()

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-c", "--genotype_call_file", help="GELI file to use in generating coverage stat table", dest="genotype_filename", )
    parser.add_option("-s", "--stats_file", help="file to containing empirical genotyping stats", dest="stats_filename")
    parser.add_option("-g", "--histogram_file", help="file to containing counts of each depth of coverage", dest="hist_filename")
    parser.add_option("-m", "--max_loci", help="maximum number of loci to parse (for debugging)", default=sys.maxint, dest="max_loci", type="int")
    parser.add_option("-v", "--discovery", help="run discovery rather than genotyping calls", default=False, dest="do_discovery", action="store_true")
    parser.add_option("-e", "--evaluate", help="evaluate genotypes; requires a stats file and a histogram file", default=False, dest="evaluate_genotypes", action="store_true")
    parser.add_option("-d", "--debug", help="provide debugging output", default=False, dest="debug", action="store_true")
    

    (options, args) = parser.parse_args()

    #if len(args) < 1:
    #    usage(parser)
    #genotype_filename = args[0]

    if options.evaluate_genotypes:
        print "Evaluating genotypes"
        if options.hist_filename == None:
            sys.exit("Must provide -g histogram filename option")
        if options.stats_filename == None:
            sys.exit("Must provide -s stats fliname option")
        stats_from_hist(options.hist_filename, options.stats_filename)
    else:
        print "Creating performance tables from genotypes file"
        aggregate = aggregate_stats(options.genotype_filename, options.max_loci)
        stats_filename = options.genotype_filename+".stats"
        create_coverage_stats_table(aggregate, stats_filename, options.debug)



