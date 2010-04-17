lanes_file = open("/humgen/gsa-hpprojects/FHS/indexed/production/oneOffAnalyses/coverage_4_14/docs/fhs_squid_lanes_4_15.txt")
samples_file = open("/humgen/gsa-hpprojects/FHS/indexed/production/oneOffAnalyses/coverage_4_14/docs/fhs_squid_samples_4_15.txt")
lanes_header = lanes_file.readline().split("\t")
import os
import time
DEBUG = False
DEBUG_RECORDS = 50
# dumb TSV often saves with strings like "190,12"
def str2num(strn):
    if ( strn == "" or strn == "''" or strn == '""' ):
        return -1
    elif ( strn.startswith('"') ):
        return str2num(strn.split('"')[1])
    elif ( strn.find(",") > -1):
        return str2num(strn.replace(",",""))
    elif ( strn.find(".") > -1 ):
        return float(strn)
    else:
        try:
            return int(strn)
        except ValueError:
            print("Odd format for int: "+strn)
            exit()

def getCaptureDate(flow,lane,lib):
    __f = os.listdir("/seq/picard/"+flow+"/")[0]
    checkdir = "/seq/picard/"+flow+"/"+__f+"/"+lane+"/"+lib+"/"
    earlydate = None
    for file in os.listdir(checkdir):
        ctime = os.stat(checkdir+file)
        ctime = ctime[len(ctime)-1]
        if ( earlydate == None or ctime < earlydate ):
            earlydate = ctime
    return ctime

class LaneRecord:
    def __init__(self,line):
        spline = line.strip().split("\t")
        self.flowcell = spline[lanes_header.index("Flowcell")]
        self.lane_number = spline[lanes_header.index("Lane")]
        self.sample_id = spline[lanes_header.index("External ID")]
        self.lane_id = self.flowcell+"."+self.lane_number
        self.library = spline[lanes_header.index("Library")]
        self.aligned_reads = str2num(spline[lanes_header.index("AL_PF_HQ_ALIGNED_READS")])
        self.read_length = str2num(spline[lanes_header.index("AL_MEAN_READ_LENGTH")])
        self.dup_pct = str2num(spline[lanes_header.index("DUP_PERCENT_DUPLICATION")])
        self.hs_lib_size = str2num(spline[lanes_header.index("HS_LIBRARY_SIZE")])
        self.hs_pf_uq_reads = str2num(spline[lanes_header.index("HS_PF_UNIQUE_READS")])
        self.picard_snps = str2num(spline[lanes_header.index("SNP_TOTAL_SNPS")])
        try:
            self.ic_rd2_error_rate = str2num(spline[lanes_header.index("Lane IC PCT Mean RD2 Err Rate")])
        except IndexError:
            print(self.lane_id)
            self.ic_rd2_error_rate = 0.5
        self.date_run = getCaptureDate(self.flowcell,self.lane_number,self.library)

class LaneRecordAggregator:
    def __init__(self,by_type):
        self.by = by_type
        self.entries = dict()

    def inc(self,rec):
        if ( self.by == "sample" ):
            if ( rec.sample_id in self.entries.keys() ):
                self.entries[rec.sample_id].add(rec)
            else:
                self.entries[rec.sample_id] = set()
                self.entries[rec.sample_id].add(rec)
        elif ( self.by == "flowcell" ):
            if ( rec.flowcell in self.entries.keys() ):
                self.entries[rec.flowcell].add(rec)
            else:
                self.entries[rec.flowcell] = set()
                self.entries[rec.flowcell].add(rec)
        elif ( self.by == "lane_number" ):
            if ( rec.lane_number in self.entries.keys() ):
                self.entries[rec.lane_number].add(rec)
            else:
                self.entries[rec.lane_number] = set()
                self.entries[rec.lane_number].add(rec)
        elif ( self.by == "library" ):
            if ( rec.library in self.entries.keys() ):
                self.entries[rec.lane_number].add(rec)
            else:
                self.entries[rec.lane_number] = set()
                self.entries[rec.lane_number].add(rec)
        else:
            print("You use sample, flowcell, or lane_number. You should hit ctrl-c now.")

    def summaryTable(self):
        strOut = "\tNumBarcodes\tDate\tReads\tLength\tDup_Pct\tLib_Size\tHS_pf_unique_reads\trd2_error_rate\tsnps"
        for id in self.entries.keys():
            strOut += "\n"+id
            lanes = self.entries[id]
            avg_date = 0.0
            avg_reads = 0.0
            avg_length = 0.0
            avg_dup = 0.0
            avg_lib_size = 0.0
            avg_hs_pf = 0.0
            avg_err_rt = 0.0
            avg_snps = 0.0
            for lane in lanes:
                avg_date += lane.date_run
                avg_reads += lane.aligned_reads
                avg_length += lane.read_length
                avg_dup += lane.dup_pct
                avg_lib_size += lane.hs_lib_size
                avg_hs_pf += lane.hs_pf_uq_reads
                avg_err_rt += lane.ic_rd2_error_rate
                avg_snps += lane.picard_snps
            avg_date = avg_date/len(lanes)
            avg_reads = avg_reads/len(lanes)
            avg_length = avg_length/len(lanes)
            avg_dup = avg_dup/len(lanes)
            avg_lib_size = avg_lib_size/len(lanes)
            avg_hs_pf = avg_hs_pf/len(lanes)
            avg_err_rt = avg_err_rt/len(lanes)
            avg_snps = avg_snps/len(lanes)
            strOut += "\t"+"\t".join([str(len(lanes)),str(time.asctime(time.gmtime(int(round(avg_date))))),str(avg_reads),str(avg_length),str(avg_dup),str(avg_lib_size),str(avg_hs_pf),str(avg_err_rt),str(avg_snps)])
        return strOut

by_flowcell_metrics = LaneRecordAggregator("flowcell")
by_lane_metrics = set()
by_lane_number_metrics = LaneRecordAggregator("lane_number")
by_sample_metrics = LaneRecordAggregator("sample")
by_library_metrics = LaneRecordAggregator("library")
line_no = 0
for line in lanes_file.readlines():
    lane = LaneRecord(line)  
    by_lane_metrics.add(lane)
    by_flowcell_metrics.inc(lane)
    by_lane_number_metrics.inc(lane)
    by_sample_metrics.inc(lane)
    by_library_metrics.inc(lane)
    line_no += 1
    if ( DEBUG and line_no > DEBUG_RECORDS ):
        break
    if ( line_no % 100 == 0 ):
        print("Read: "+str(line_no)+" lines.")

out = open("flowcell_metrics.txt",'w')
out.write(by_flowcell_metrics.summaryTable())
out.close()
out = open("sample_metrics.txt",'w')
out.write(by_sample_metrics.summaryTable())
out.close()
out = open("lane_number_metrics.txt",'w')
out.write(by_lane_number_metrics.summaryTable())
out.close()
out = open("library_metrics.txt",'w')
out.write(by_library_metrics.summaryTable())
out.close()
out = open("per_lane_metrics.txt",'w')
out.write("\t"+"\t".join(["lane_id","date","reads","length","dups","lib_size","hs_reads_pf_uq","rt2_err","snps"]))
for lane in by_lane_metrics:
    out.write("\n"+"\t".join([str(lane.lane_id),str(time.asctime(time.gmtime(int(round(lane.date_run))))),str(lane.aligned_reads),str(lane.dup_pct),str(lane.hs_lib_size),str(lane.hs_pf_uq_reads),str(lane.ic_rd2_error_rate),str(lane.picard_snps)]))
out.close()
