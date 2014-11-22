/* 
The MIT Lincense

Copyright (c) 2014 Takaho A. Endo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
#include <set>
#include <string>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>

using namespace std;
//using namespace std::stringstream;

#include <gtf.hxx>
#include <tktools.hxx>

using namespace tkbio;
using namespace tktools::bio;


string* gtfexon::_features = NULL;
const int gtfexon::_num_features = 10;
namespace {
    string* generate_feature_array() {
        string* _features = new string[10];
        _features[gtfexon::OTHERS] = string(".");
        _features[gtfexon::EXON] = "exon";
        _features[gtfexon::CDS] = "CDS";
        _features[gtfexon::START_CODON] = "start_codon";
        _features[gtfexon::STOP_CODON] = "stop_codon";
        _features[gtfexon::UTR] = "UTR";
        _features[gtfexon::transcript] = "transcript";
        _features[gtfexon::gene] = "gene";
        _features[gtfexon::Selenocystein] = "Selenocystein";
        return _features;
    }
}

const string& gtfexon::get_feature(Feature _feature) {
    if (_features == NULL) {
        _features = generate_feature_array();//new string[_num_features];
        // _features[0] = string(".");
        // _features[1] = "exon";
        // _features[2] = "CDS";
        // _features[3] = "start_codon";
        // _features[4] = "stop_codon";
    }
    if (_feature < 0 || _feature >= _num_features) {
        return _features[0];
    } else {
        return _features[_feature];
    }
}

namespace {
    set<string> _others;
}

gtfexon::Feature gtfexon::encode_feature(const char* feature) {
    if (_features == NULL) {
        //cerr << "instanciate _features\n";
        _features = generate_feature_array();
    }
    for (int i = 1; i < _num_features; i++) {
        //cerr << i << " : " << _features[i] << endl;
        if (strncmp(feature, _features[i].c_str(), _features[i].size()) == 0) {
            return (Feature)i;
        }
    }
    if (_others.find(string(feature)) == _others.end()) {
        cerr << "new feature : " << feature << endl;
        _others.insert(string(feature));
    }
    return OTHERS;
}

gtfgene::gtfgene(const string& transcript_id, const string& name, const string& tss_id) {
    _name = name;
    _transcript_id = transcript_id;
    _tss_id = tss_id;
    _orientation = '\0';
    _position3 = _position5 = -1;
    _chromosome_code = -1;
}

string gtfgene::to_string() const {
    std::stringstream ss;
    ss << _name << "\t" << _transcript_id << "\t" << _chromosome << ":" << _orientation << ":" << _position5 << "-" << _position3;
    return ss.str();
}

void gtfgene::insert(const string& chromosome, char orientation, int pos5, int pos3, const char* feature) {
    _chromosome = chromosome;
    _chromosome_code = convert_chromosome_to_code(chromosome.c_str());
    _orientation = orientation;
    if (_position5 < 0 || _position5 > pos5) _position5 = pos5;
    if (_position3 < 0 || _position3 < pos3) _position3 = pos3;
    //cerr << feature << endl;
    gtfexon::Feature feature_id = gtfexon::encode_feature(feature);
    //cerr << feature_id << endl;
    if (feature_id == gtfexon::EXON
        || feature_id == gtfexon::UTR
        || feature_id == gtfexon::CDS) {
        _exons.push_back(gtfexon(pos5, pos3, feature_id));
    }
}

void gtfgene::determine_span() {
    if (_position5 < 0) {
        _position5 = numeric_limits<int>::max();
        _position3 = 0;
        for (vector<gtfexon>::const_iterator it = _exons.begin(); it != _exons.end(); it++) {
            int p5 = it->_pos5;
            int p3 = it->_pos3;
            if (p5 < _position5) _position5 = p5;
            if (p3 > _position3) _position3 = p3;
        }
    }
}

void gtfgene::insert(int pos5, int pos3, const char* feature) {
    _exons.push_back(gtfexon(pos5, pos3, feature));
    if (_position5 < 0 || _position5 > pos5) _position5 = pos5;
    if (_position3 < 0 || _position3 < pos3) _position3 = pos3;
}

pair<int,int> gtfgene::get_span() const {
    //if (_position5 < 0) {
    const_cast<gtfgene*>(this)->determine_span();
    //}
    return make_pair(_position5, _position3);
}

bool gtfgene::compare_position(const gtfgene* lhs, const gtfgene* rhs) {
    if (lhs->_chromosome_code != rhs->_chromosome_code) {
        return lhs->_chromosome_code < rhs->_chromosome_code;
    }
    return lhs->_position5 < rhs->_position5;
}

vector<gtfgene*> gtffile::get_tss() const {
    vector<gtfgene*> tss;
    for (map<string, vector<gtfgene*> >::const_iterator it = _tssid.begin();
         it != _tssid.end(); it++) {
        tss.push_back((it->second)[0]);
    }
    return tss;
}

vector<string> gtffile::get_symbols() const {
    vector<string> genes;
    for (map<string, vector<gtfgene*> >::const_iterator it = _symbol.begin();
         it != _symbol.end(); it++) {
        genes.push_back(it->first);
    }
    return genes;
}

vector<gtfgene*> gtffile::get_genes() const {
    vector<gtfgene*> genes;
    for (map<string, gtfgene*>::const_iterator it = _genes.begin();
         it != _genes.end(); it++) {
        genes.push_back(it->second);
    }
    return genes;
}

gtffile::gtffile() {
    _sorted = false;
}

void gtffile::add(const char* transcript_id, const char* name, const char* tss_id, const char* feature, const char* chromosome, char orientation, int position5, int position3) throw (exception) {
    map<string, gtfgene*>::iterator it = _genes.find(transcript_id);
    gtfgene* gene;
    if (it == _genes.end()) {
        gene = new gtfgene(transcript_id, name, tss_id);
        _genes[transcript_id] = gene;
        if (_tssid.find(tss_id) == _tssid.end()) {
            vector<gtfgene*> tmp_genes;
            tmp_genes.push_back(gene);
            _tssid[tss_id] = tmp_genes;
        } else {
            _tssid[tss_id].push_back(gene);
        }
        if (_symbol.find(name) == _symbol.end()) {
            vector<gtfgene*> tmp_genes;
            tmp_genes.push_back(gene);
            _symbol[name] = tmp_genes;
        } else {
            _symbol[name].push_back(gene);
        }
        _sorted = false;
        gene->insert(chromosome, orientation, position5, position3, feature);
    } else {
        gene = it->second;
        gene->insert(position5, position3, feature);
    }
}

void gtffile::sort_genes() {
    if (_sorted) return;
    _maximum_gene_span = 0;
    _sorted_genes.erase(_sorted_genes.begin(), _sorted_genes.end());
    for (map<string, gtfgene*>::iterator it = _genes.begin(); it != _genes.end(); it++) {
        gtfgene* g = it->second;
        int span = g->_position3 - g->_position5 + 1;
        if (span > _maximum_gene_span) { _maximum_gene_span = span; }
        _sorted_genes.push_back(g);
    }
    std::sort(_sorted_genes.begin(), _sorted_genes.end(), gtfgene::compare_position);
    _sorted = true;
}

gtffile* gtffile::load_gtf(const char* filename) throw (exception) {
    if (tktools::io::file_exists(filename) == false) {
        throw logic_error(string("cannot open ") + filename);
    }
    gtffile* data = new gtffile();
    ifstream fi(filename);
    while (!fi.eof()) {
        string line;
        getline(fi, line);
        int col = 0;
        int span = line.size();
        char* ptr = new char[line.size() + 1];
        memcpy(ptr, line.c_str(), line.size() + 1);
        int prev = 0;
        const char* chromosome = ptr;
        const char* feature = NULL;
        char ori = '\0';
        int p3 = -1, p5 = -1;
        for (int i = 0; i <= span; i++) {
            char c = ptr[i];
            if (c < ' ') {
                //cout << line.substr(prev, i - prev) << endl;
                //ptr[i] = '\0';
                //cerr << col << ":" << ptr + prev << endl;
                if (col == 0) { // chromosome
                    ptr[i] = '\0';
                } else if (col == 2) { // feature
                    ptr[i] = '\0';
                    feature = ptr + prev;
                    //cerr << "FEATURE====" << feature << endl;
                } else if (col == 3) { // 
                    p5 = atoi(ptr + prev);
                } else if (col == 4) {
                    p3 = atoi(ptr + prev);
                } else if (col == 6) {
                    ori = ptr[prev];
                    //cerr << "ORI: " << ori << endl;
                } else if (col == 8) {
                    //cerr << feature << endl;
                    //cerr << col << ":" << string(ptr + prev, i - prev) << endl;
                    const char* key = ptr + prev;// + i + 1;
                    const char* value = NULL;
                    const char* tss_id = NULL;
                    const char* transcript_id = NULL;
                    const char* name = NULL;
                    //cerr << "j=" << (prev) << "-" << span << endl;

                    for (int j = prev; j < span; j++) {
                        c = ptr[j];
                        if (c == ';') { // end of an item
                            key = value = NULL;
                        } else if (c == '"') { // value
                            if (value == NULL) {
                                value = ptr + j + 1;
                                //cerr << "value=" << value << endl;
                            } else {
                                ptr[j] = '\0';
                                if (key != NULL) {
                                    //cerr << key << " => " << value << endl;
                                    if (strncmp(key, "tss_id", 6) == 0) {
                                        tss_id = value;
                                        //cerr << "TSS_ID:" << tss_id << endl;
                                    } else if (strncmp(key, "gene_name", 9) == 0) {
                                        name = value;
                                        //cout << "GENE_NAME:" << name << endl;
                                    } else if (strncmp(key, "transcript_id", 13) == 0) {
                                        transcript_id = value;
                                        //cerr << "TRANSCRIPT_ID:" << transcript_id << endl;
                                    }
                                }
                            }
                        } else if (c > ' ') { // 
                            if (key == NULL) {
                                key = ptr + j;
                                //cerr << key << endl;
                            }
                        }
                    }

                    if (name != NULL && transcript_id != NULL) {
                        if (tss_id == NULL) {
                            tss_id = transcript_id;
                        }
                        data->add(transcript_id, name, tss_id, feature, chromosome, ori, p5, p3);
                    }
                    break;
                }
                col++;
                prev = i + 1;
            }
        }
        delete[] ptr;
    }
    fi.close();

    //cerr << "stat : " << data->_genes.size() << " genes, " << data->_tssid.size() << " TSSs, " << data->_symbol.size() << " symbols" << endl;
    return data;
}

gtffile::~gtffile() {
    for (map<string,gtfgene*>::iterator it = _genes.begin(); it != _genes.end(); it++) {
        delete it->second;
    }
}

bool gtfgene::contains_in_exon(int chromosome_code, int location) const {
    return contains_in_exon(chromosome_code, location, location);
}

bool gtfgene::contains_in_exon_debug(int ccode, int location) const {
    if (ccode != _chromosome_code || location< _position5 || location > _position3) {
        return false;
    }
    for (int i = 0; i < (int)_exons.size(); i++) {
        const gtfexon& e = _exons[i];
        if (e._feature == gtfexon::EXON) {
            cout << e.feature_name() << ":" << e._pos5 << "-" << e._pos3 << "  " << location << endl;
            if (e._pos5 <= location && location <= e._pos3) {
                return true;
            } else if (location > e._pos3) {
                return false;
            }
        }
    }
    return false;
}

bool gtfgene::contains(int chromosome_code, int position) const {
    if (chromosome_code != _chromosome_code || position < _position5 || position > _position3) {
        return false;
    }
    for (int i = 0; i < (int)_exons.size(); i++) {
        const gtfexon& e = _exons[i];
        if (e._pos5 <= position && position <= e._pos3) {
            return true;
        }
    }
    return false;
}

bool gtfgene::contains_in_exon(int chromosome_code, int start, int end) const {
    if (chromosome_code != _chromosome_code || end < _position5 || start > _position3) {
        return false;
    }
    for (int i = 0; i < (int)_exons.size(); i++) {
        const gtfexon& e = _exons[i];
        if (e._feature == gtfexon::EXON) {
            if (e._pos5 <= end && start <= e._pos3) {
                return true;
            }
        }
    }
    return false;
}


namespace {
    bool paircomp(const pair<int, int>& lhs, const pair<int, int>& rhs) {
        return lhs.first < rhs.first;
    }
}

vector<const gtfgene*> gtffile::find_genes(const string& chromosome, int start, int end) const {
    int ccode = convert_chromosome_to_code(chromosome.c_str());
    return find_genes(ccode, start, end);
}

const gtfgene* gtffile::get_container(int ccode, int pos) const {
    const_cast<gtffile*>(this)->sort_genes();
    int left = 0;
    int right = _sorted_genes.size();
    //cout << ccode << ":" << pos << " // " << _sorted_genes.size() << endl;
    if (right == 0) return NULL;
    int center = -1;
    for (;;) {
        center = (left + right) / 2;
        const gtfgene* g = _sorted_genes[center];
        //cout << center << " " << g->chromosome() << ":" << g->_orientation << ":" << g->_position5 << "-" << g->_position3 << " " << g->name() << " // " << ccode << ":" << pos << endl;
        if (ccode < g->_chromosome_code) {
            right = center;
        } else if (ccode > g->_chromosome_code) {
            left = center + 1;
        } else if (pos < g->_position5 - _maximum_gene_span) {
            right = center;
        } else if (pos > g->_position3 + _maximum_gene_span) {
            left = center + 1;
        } else {
            break;
        }
        if (left == right) {
            center = -1;
            break;
        }
    }
    if (center < 0) return NULL;
    for (int index = center; index >= 0; index--) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position3 + _maximum_gene_span < pos) {
            break;
        } else if (g->contains(ccode, pos)) {
            return g;
        }
    }
    for (int index = center + 1; index < (int)_sorted_genes.size(); index++) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position5 - _maximum_gene_span > pos) {
            break;
        } else if (g->contains(ccode, pos)) {
            return g;
        }
    }
    return NULL;
}

bool gtffile::contains(int ccode, int pos) const {
    return get_container(ccode, pos) != NULL;
}

vector<const gtfgene*> gtffile::find_genes(int ccode, int start, int end) const {
    vector<const gtfgene*> selected;
    //cout << ccode << ":" << start << "-" << end << endl;
    if (ccode < 0) {
        vector<const gtfgene*> selected;
        return selected;
    }

    if (end < 0) end = start + 1;
    const_cast<gtffile*>(this)->sort_genes();
    int left = 0; 
    int right = _sorted_genes.size();
    int pivot = -1;
    while (left < right) {
        int center = (left + right) >> 1;
        const gtfgene* g = _sorted_genes[center];
        if (ccode < g->_chromosome_code) {//chromosome) {
            right = center;
        } else if (ccode > g->_chromosome_code) {//chromosome) {
            left = center + 1;
        } else if (end < g->_position5 - _maximum_gene_span) {
            right = center;
        } else if (start > g->_position3 + _maximum_gene_span) {
            left = center + 1;
        } else {
            pivot = center;
            break;
        }
    }

    if (pivot < 0) {
        return selected;
    }
    for (int index = pivot; index >= 0; index--) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position3 + _maximum_gene_span < start) {
            break;
        } else if (start <= g->_position3 && g->_position5 <= end) {
            selected.push_back(const_cast<const gtfgene*>(g));
        }
    }
    for (int index = pivot + 1; index < (int)_sorted_genes.size(); index++) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position5 - _maximum_gene_span > end) {
            break;
        } else if (start <= g->_position3 and g->_position5 <= end) {
            selected.push_back(const_cast<const gtfgene*>(g));
        }
    }

    return selected;
}

bool gtffile::contains_in_exon(const string& chromosome, int location) const {
    int ccode = convert_chromosome_to_code(chromosome.c_str());
    const_cast<gtffile*>(this)->sort_genes();
    int left = 0; 
    int right = _sorted_genes.size();
    int pivot = -1;
    while (left < right) {
        int center = (left + right) >> 1;
        const gtfgene* g = _sorted_genes[center];
        if (g->_chromosome_code < ccode) {//chromosome) {
            left = center + 1;
        } else if (g->_chromosome_code > ccode) {//chromosome) {
            right = center;
        } else if (location < g->_position5 - _maximum_gene_span) {
            right = center;
        } else if (location > g->_position3 + _maximum_gene_span) {
            left = center + 1;
        } else {
            pivot = center;
            break;
        }
    }

    if (pivot < 0) {
        return false;
    }
    for (int index = pivot; index >= 0; index--) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position3 + _maximum_gene_span < location) {
            break;
        }
        if (g->contains_in_exon(ccode, location)) return true;
    }
    for (int index = pivot + 1; index < (int)_sorted_genes.size(); index++) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome_code != ccode || g->_position5 - _maximum_gene_span > location) {
            break;
        }
        if (g->contains_in_exon(ccode, location)) return true;
    }
    return false;
}

bool gtffile::contains_in_exon_debug(const string& chromosome, int location) const {
    int ccode = convert_chromosome_to_code(chromosome.c_str());
    const_cast<gtffile*>(this)->sort_genes();
    int left = 0; 
    int right = _sorted_genes.size();
    int pivot = -1;
    while (left < right) {
        int center = (left + right) >> 1;
        const gtfgene* g = _sorted_genes[center];
        if (g->_chromosome_code < ccode) {//chromosome) {
            left = center + 1;
        } else if (g->_chromosome_code > ccode) {//chromosome) {
            right = center;
        } else if (location < g->_position5 - _maximum_gene_span) {
            right = center;
        } else if (location > g->_position3 + _maximum_gene_span) {
            left = center + 1;
        } else {
            pivot = center;
            break;
        }
    }

    if (pivot < 0) {
        return false;
    }
    for (int index = pivot; index >= 0; index--) {
        const gtfgene* g = _sorted_genes[index];
        cout << g->transcript_id() << ":" << g->name() << endl;
        if (g->_chromosome_code != ccode || g->_position3 + _maximum_gene_span < location) {
            break;
        }
        if (g->contains_in_exon_debug(ccode, location)) return true;
    }
    for (int index = pivot + 1; index < (int)_sorted_genes.size(); index++) {
        const gtfgene* g = _sorted_genes[index];
        cout << g->transcript_id() << ":" << g->name() << endl;
        if (g->_chromosome_code != ccode || g->_position5 - _maximum_gene_span > location) {
            break;
        }
        if (g->contains_in_exon_debug(ccode, location)) return true;
    }
    return false;
}

vector<pair<int,int> > gtffile::get_exonregion(const string& chromosome, bool utr) const {
    int ccode = convert_chromosome_to_code(chromosome.c_str());
    vector<pair<int,int> > regions;
    if (ccode < 0) {
        return regions;
    }
    for (map<string,gtfgene*>::const_iterator it = _genes.begin(); it != _genes.end(); it++) {
        //for (int i = 0; i < (int)_genes.size(); i++) {
        const gtfgene* gene = it->second;//_genes[i];
        if (gene->_chromosome_code == ccode) {
            vector<gtfexon> exons = gene->exons();
            for (int j = 0; j < (int)exons.size(); j++) {
                if (!utr && exons[j].feature() != gtfexon::CDS) {
                    continue;
                }
                regions.push_back(make_pair(exons[j].position5(), exons[j].position3()));
            }
        }
    }
    if (regions.size() <= 1) {
        return regions;
    }

    // merge overlaps
    std::sort(regions.begin(), regions.end(), paircomp);
    vector<pair<int,int> > merged;
    int start = regions[0].first;
    int end = regions[1].second;
    for (int i = 1; i < (int)regions.size(); i++) {
        int s = regions[i].first;
        int e = regions[i].second;
        if (start < e && s <= end) {
            end = e < end ? end : e;
        } else {
            merged.push_back(make_pair(start, end));
            start = s;
            end = e;
        }
    }
    merged.push_back(make_pair(start, end));
    return merged;
}
