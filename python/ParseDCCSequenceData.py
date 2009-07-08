#!/usr/bin/env python

import operator, FlatFileTable

class db_file:
    filename = ""

    def __init__(self, filenm):
        self.filename = filenm

    def count_fields(self, fixed_field, fixed_field_values, count_field):
        record_gen = FlatFileTable.record_generator(self.filename)

        counts = dict()
        #fixed_field_num = self.field_names[fixed_field]
        print count_field+" for "+fixed_field+" = "+" or ".join(fixed_field_values)

        for record in record_gen:
            if record[fixed_field] in fixed_field_values:
                #fixed_field_value = fields[fixed_field_num]
                count_field_num = record[count_field]
                if counts.has_key(count_field_num):
                    counts[count_field_num] += 1
                else:
                    counts[count_field_num] = 0

        for k,v in sorted(counts.items(), key=operator.itemgetter(1), cmp=lambda x,y: y-x ):
            print str(k)+"\t"+str(v)

    def count_bases(self, fixed_field_values): #, fixed_field_values):
        record_gen = FlatFileTable.record_generator(self.filename)

        base_count = 0
        #fixed_field_num = self.field_names[fixed_field]
        #print "For "+fixed_field+" = "+" or ".join(fixed_field_values)+":",
        print "For "+ " AND ".join( [one_ffv[0]+" = "+" OR ".join(one_ffv[1]) for one_ffv in fixed_field_values] )

        for record in record_gen:
            #if record[fixed_field] in fixed_field_values:
            if FlatFileTable.record_matches_values(record, fixed_field_values):
                try:
                    base_count += int(record["BASE_COUNT"])
                except ValueError:
                    pass

        print "%e bases" % base_count

if __name__ == "__main__":
    db = db_file("sequence.index")
    
    platforms = (("ILLUMINA",), ("AB SOLiD","SOLID","ABI_SOLID","AB SOLiD System 2.0"), ("LS454",))
    studies = (("1000Genomes Project Pilot 1",), ("1000Genomes Project Pilot 2",), ("1000Genomes Project Pilot 3",))

    for select_field, select_field_values in (): #(("INSTRUMENT_PLATFORM", platforms), ("STUDY_NAME", studies)):
        for count_field in ("CENTER_NAME", "STUDY_NAME", "INSTRUMENT_PLATFORM"):
            for select_field_value in select_field_values:
                db.count_fields(select_field, select_field_value, count_field)
                
        for select_field_value in select_field_values:
            db.count_bases(((select_field, select_field_value),))

        print

    for field1, value1 in zip(["INSTRUMENT_PLATFORM"]*len(platforms), platforms):
        for field2, value2 in zip(["STUDY_NAME"]*len(studies), studies):
            db.count_bases(((field1, value1), (field2, value2)))









