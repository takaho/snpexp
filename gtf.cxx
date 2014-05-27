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

using namespace std;

#include <gtf.hxx>

using namespace tkbio;

namespace {
    string* generate_feature_array() {
        string* _features = new string[5];
        _features[gtfexon::OTHERS] = string(".");
        _features[gtfexon::EXON] = "exon";
        _features[gtfexon::CDS] = "CDS";
        _features[gtfexon::START_CODON] = "start_codon";
        _features[gtfexon::STOP_CODON] = "stop_codon";
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
    if (_feature < 0 || _feature >= 4) {
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
        _features = generate_feature_array();
    }
    for (int i = 1; i < _num_features; i++) {
        if (strncmp(feature, _features[i].c_str(), _features[i].size()) == 0) {
            return (Feature)i;
        }
    }
    if (_others.find(string(feature)) == _others.end()) {
        cerr << "new feature : " << feature;
        _others.insert(string(feature));
    }
    return OTHERS;
}

string* gtfexon::_features = NULL;
const int gtfexon::_num_features = 5;

gtfgene::gtfgene(const string& transcript_id, const string& name, const string& tss_id) {
    _name = name;
    _transcript_id = transcript_id;
    _tss_id = tss_id;
    _orientation = '\0';
    _position3 = _position5 = -1;
}

void gtfgene::insert(const string& chromosome, char orientation, int pos5, int pos3, const char* feature) {
    _chromosome = chromosome;
    _orientation = orientation;
    if (_position5 < 0 || _position5 > pos5) _position5 = pos5;
    if (_position3 < 0 || _position3 < pos3) _position3 = pos3;
    _exons.push_back(gtfexon(pos5, pos3, feature));
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
    if (lhs->_chromosome != rhs->_chromosome) {
        return lhs->_chromosome < rhs->_chromosome;
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
    //cerr << "add " << transcript_id << " " << name << " " << tss_id << " " << feature << " " << chromosome << ":" << orientation << ":" << position5 << "-" << position3 << "   " << __LINE__ <<  endl;
    map<string, gtfgene*>::iterator it = _genes.find(transcript_id);
    gtfgene* gene;
    if (it == _genes.end()) {
        gene = new gtfgene(transcript_id, name, tss_id);
        _genes[transcript_id] = gene;
        if (_tssid.find(tss_id) == _tssid.end()) {
            //cerr << "new tss : " <<  tss_id << endl;
            vector<gtfgene*> tmp_genes;
            tmp_genes.push_back(gene);
            _tssid[tss_id] = tmp_genes;
            //gene);
        } else {
            _tssid[tss_id].push_back(gene);
        }
        if (_symbol.find(name) == _symbol.end()) {
            //cerr << "new symbol " << name << endl;
            vector<gtfgene*> tmp_genes;
            tmp_genes.push_back(gene);
            _symbol[name] = tmp_genes;
            //gene);
        } else {
            _symbol[name].push_back(gene);
        }
        _sorted = false;
        //if (_maximum_gene_span < span) _maximum_gene_span = span;
        // _symbol[name].push_back(gene);
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

gtffile* gtffile::load(const char* filename) throw (exception) {
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
                                        //cerr << "GENE_NAME:" << name << endl;
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

                    if (tss_id != NULL && name != NULL && transcript_id != NULL) {
                        //cout << chromosome << endl;
                        if (strncmp(chromosome, "chr", 3) == 0 || strncmp(chromosome, "Chr", 3) == 0) {
                            chromosome += 3;
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

namespace {
    bool paircomp(const pair<int, int>& lhs, const pair<int, int>& rhs) {
        return lhs.first < rhs.first;
    }
}

vector<const gtfgene*> gtffile::find_genes(string chromosome, int start, int end) const {
    if (chromosome.find("chr") == 0) {
        chromosome = chromosome.substr(3);
    }
    if (end < 0) end = start + 1;
    const_cast<gtffile*>(this)->sort_genes();
    int left = 0; 
    int right = _sorted_genes.size();
    int pivot = -1;
    while (left < right) {
        int center = (left + right) >> 1;
        const gtfgene* g = _sorted_genes[center];
        if (g->_chromosome < chromosome) {
            left = center + 1;
        } else if (g->_chromosome > chromosome) {
            right = center;
        } else if (end < g->_position5 - _maximum_gene_span) {
            right = center;
        } else if (start > g->_position3 + _maximum_gene_span) {
            left = center + 1;
        } else {
            pivot = center;
            break;
        }
    }

    vector<const gtfgene*> selected;
    if (pivot < 0) {
        return selected;
    }
    for (int index = pivot; index >= 0; index--) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome != chromosome || g->_position3 + _maximum_gene_span < start) {
            break;
        } else if (start <= g->_position3 && g->_position5 <= end) {
            selected.push_back(const_cast<const gtfgene*>(g));
        }
    }
    for (int index = pivot + 1; index < (int)_sorted_genes.size(); index++) {
        const gtfgene* g = _sorted_genes[index];
        if (g->_chromosome != chromosome || g->_position5 - _maximum_gene_span > end) {
            break;
        } else if (start <= g->_position3 and g->_position5 <= end) {
            selected.push_back(const_cast<const gtfgene*>(g));
        }
    }


    // for (map<string, gtfgene*>::const_iterator it = _genes.begin(); it != _genes.end(); it++) {
    //     const gtfgene* g = it->second;
    //     if (g->_chromosome == chromosome) {
    //         //cerr << chromosome << ":" << start << "-" << end << " <-> " << g->_chromosome << ":" << g->_position5 << "-" << g->_position3 << endl;
    //         if (start <= g->_position3 and g->_position5 <= end) {
    //             selected.push_back(const_cast<const gtfgene*>(g));
    //         }
    //     }
    // }
    return selected;
}

vector<pair<int,int> > gtffile::get_exonregion(string chromosome, bool utr) const {
    if (chromosome.find("chr") == 0) {
        chromosome = chromosome.substr(3);
    }
    //const vector<gtfgene*>& genes = gtf->get_genes();
    vector<pair<int,int> > regions;
    for (map<string,gtfgene*>::const_iterator it = _genes.begin(); it != _genes.end(); it++) {
        //for (int i = 0; i < (int)_genes.size(); i++) {
        const gtfgene* gene = it->second;//_genes[i];
        if (gene->chromosome() == chromosome) {
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
