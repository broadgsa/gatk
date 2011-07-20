#include "MainTools.h"
#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "lookup/SerialQltout.h"

unsigned int MatchingEnd(look_align &la, vecbasevector &candidates, vecbasevector &ref) {
    //la.PrintParseable(cout);

    for (int i = 0; i < candidates.size(); i++) {
        look_align newla = la;

        if (newla.rc1) { candidates[i].ReverseComplement(); }
        newla.ResetFromAlign(newla.a, candidates[i], ref[la.target_id]);

        //newla.PrintParseable(cout, &candidates[i], &ref[newla.target_id]);
        //cout << newla.Errors() << " " << la.Errors() << endl;

        if (newla.Errors() == la.Errors()) {
            return i;
        }
    }

    //FatalErr("Query id " + ToString(la.query_id) + " had no matches.");

    return candidates.size() + 1;
}

int main(int argc, char **argv) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String(ALIGNS);
    CommandArgument_String(FASTB_END_1);
    CommandArgument_String(FASTB_END_2);
    CommandArgument_String(REFERENCE);

    CommandArgument_String(ALIGNS_END_1_OUT);
    CommandArgument_String(ALIGNS_END_2_OUT);
    EndCommandArguments;

    vecbasevector ref(REFERENCE);
    vecbasevector reads1(FASTB_END_1);
    vecbasevector reads2(FASTB_END_2);

    ofstream aligns1stream(ALIGNS_END_1_OUT.c_str());
    ofstream aligns2stream(ALIGNS_END_2_OUT.c_str());

    basevector bv;

    SerialQltout sqltout(ALIGNS);
    look_align la;
    while (sqltout.Next(la)) {
        vecbasevector candidates(2);
        candidates[0] = reads1[la.query_id];
        candidates[1] = reads2[la.query_id];

        unsigned int matchingend = MatchingEnd(la, candidates, ref);
        if (matchingend < 2) {
            bv = (matchingend == 0) ? reads1[la.query_id] : reads2[la.query_id];

            //la.PrintParseable(cout, &bv, &ref[la.target_id]);
            la.PrintParseable(((matchingend == 0) ? aligns1stream : aligns2stream), &bv, &ref[la.target_id]);
        }
    }

    aligns1stream.close();
    aligns2stream.close();

    return 0;
}
