#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::ostream;

#include <fragmentprocessor.hxx>
#include <recrec.hxx>

using namespace tkbio;

recombination_detector::recombination_detector(int coverage, float hetero_threshold){
    _coverage = coverage;
    _threshold = hetero_threshold;
}

void recombination_detector::process_fragments(const vector<recfragment*>& fragments,
                                               chromosome_seq const* chromosome,
                                               int start, int end, ostream& ost) throw (exception) {
    if (chromosome != NULL) {
        std::cout << chromosome->name() << ";" << chromosome->length() << "\t" << start << "\t" << end << "\t" << fragments.size() << std::endl;
    }
}
