#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>


//  Primitive MS/MS search engine
//  Copyright (C) 2006 Jacques Colinge

//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

//  Contact:
//   Dr. Jacques Colinge
//   CeMM
//   Lazarettgasse 19/3
//   A-1090 Vienna, Austria
//   http://www.cemm.at


// Simple C++ code to implement the PMF DB search algorithm. Some elementary classes are defined to
// collect and organize data. No real OOP! Compiled with g++ and tested under Linux. Read pmfDBSearch
// code first to understand this code more easily.

// No variable modifications are possible, cysteines are modified. Tryptic digestion only.

// Exercises: (1) Dancik score (Dancik et al., J Comp Biol, 1999);
//            (2) Improve the p-value estimation;
//            (3) Implement variable oxidation on methionines;
//            (4) Implement a simple user interface (web-based).


using namespace std;

// Global variables to store basic quantities
int verbose = 0;
enum {
    noScore = -300, nBest = 10
};
vector<double> aaMasses;
double waterMass, protonMass, hydrogenMass, oxygenMass, carbonMass, ammoniaMass, nitrogenMass;
double minMass = 200, maxMass = 3500;
double tol;


class PeptidePosition {
public:
    PeptidePosition() : first(0), last(0), length(0), nmc(0) {}

    PeptidePosition(unsigned f, unsigned l, unsigned len, unsigned n) {
        first = f;
        last = l;
        length = len;
        nmc = n;
    }

    ~PeptidePosition() {}

    unsigned first, last, length, nmc;

}; // class PeptidePosition


class Peak {
public:
    Peak() : moz(0.0), intensity(0.0), charge(0) {}

    Peak(const double m, const double i, const int c) {
        moz = m;
        intensity = i;
        charge = c;
    }

    ~Peak() {}

    double moz, intensity;
    int charge;

}; // class Peak


bool operator<(const Peak &a, const Peak &b) {
    return a.moz < b.moz;

} // operator<(const Peak&, const Peak&)


class ExpPeptide {
public:
    ExpPeptide() {
        description = "";
    }

    ~ExpPeptide() {}

    Peak parentIon;
    vector<Peak> fragments;
    string description;
    unsigned order;

}; // class ExpPeptide


bool operator<(const ExpPeptide &a, const ExpPeptide &b) {
    return a.parentIon.moz < b.parentIon.moz;

} // operator<(const ExpPeptide&, const ExpPeptide&)


class Mass2Pept {
public:
    Mass2Pept() {}

    Mass2Pept(const double m, ExpPeptide *ep) {
        mass = m;
        expPeptide = ep;
    }

    ~Mass2Pept() {};

    double mass;
    ExpPeptide *expPeptide;

}; // class Mass2Pept


bool operator<(const Mass2Pept &a, const Mass2Pept &b) {
    return a.mass < b.mass;

} // operator<(const Mass2Pept&, const Mass2Pept&)


class PeptideMatch {
public:
    PeptideMatch() {
        ac = "";
        sequence = "";
        pos = 0;
        score = 0.0;
        pValue = 1.0;
        expPeptide = 0;
        first = last = 0;
    }

    PeptideMatch(const string &a, const fstream::pos_type p, const string &s, const double sc, ExpPeptide *ep,
                 const unsigned f, const unsigned l) {
        ac = a;
        pos = p;
        score = sc;
        pValue = 1.0;
        sequence = s;
        expPeptide = ep;
        first = f;
        last = l;
    }

    ~PeptideMatch() {}

    double score, pValue;
    string ac, sequence;
    ExpPeptide *expPeptide;
    fstream::pos_type pos;
    unsigned first, last;

}; // class PeptideMatch


bool operator<(const PeptideMatch &a, const PeptideMatch &b) {
    return a.score < b.score;

} // operator<(const PeptideMatch&, const PeptideMatch&)


class ProteinMatch {
public:
    ProteinMatch() {
        ac = id = descr = sequence = "";
        pos = 0;
        score = 0.0;
        coverage = 0.0;
    }

    ProteinMatch(const string &a, const fstream::pos_type p) {
        ac = a;
        pos = p;
        score = 0.0;
        id = sequence = "";
        coverage = 0.0;
    }

    ~ProteinMatch() {}

    double score, coverage;
    string ac, id, descr, sequence;
    vector<int> peptMatch;
    fstream::pos_type pos;

}; // class ProteinMatch


bool operator<(const ProteinMatch &a, const ProteinMatch &b) {
    return a.score < b.score;

} // operator<(const ProteinMatch&, const ProteinMatch&)


class TheoreticalFragmentation {
public:
    TheoreticalFragmentation() {}

    ~TheoreticalFragmentation() {}

    string type;
    vector<double> moz;

}; // TheoreticalFragmentation


unsigned startWith(string &s, char *b) {
    unsigned l = strlen(b);
    if (s.substr(0, l) == b)
        return l;
    else
        return 0;

} // startWith


bool startWithANumber(string &s) {
    string::size_type pos = s.find_first_of("0123456789");
    return pos != string::npos;

} // startWithANumber


void readMgf(const string mgf, vector<ExpPeptide *> &expPeptides) {
    string title;
    double pepMoz, pepIntens, fragMoz, fragIntens;
    int charge;
    Peak frag;
    bool inIon = false;
    ExpPeptide *ep;
    unsigned order = 0;

    ifstream pkfile(mgf.c_str());
    string line;
    unsigned v;
    getline(pkfile, line);
    while (pkfile.good()) {
        if ((v = startWith(line, "BEGIN IONS"))) {
            ep = new ExpPeptide();
            inIon = true;
        } else if ((v = startWith(line, "TITLE="))) {
            title = line.substr(v);
            ep->description = title;
        } else if ((v = startWith(line, "PEPMASS="))) {
            istringstream in(line.substr(v));
            in >> pepMoz >> pepIntens;
            ep->parentIon.moz = pepMoz;
            ep->parentIon.intensity = pepIntens;
            charge = 0;
        } else if ((v = startWith(line, "CHARGE="))) {
            if (line.substr(v) == "1+")
                charge = 1;
            else if (line.substr(v) == "2+")
                charge = 2;
            else if (line.substr(v) == "3+")
                charge = 3;
            else if (line.substr(v) == "4+")
                charge = 4;
            else
                charge = 0;
        } else if (inIon && startWithANumber(line)) {
            istringstream in(line.substr(v));
            in >> fragMoz >> fragIntens;
            frag.moz = fragMoz;
            frag.intensity = fragIntens;
            frag.charge = 0;
            ep->fragments.push_back(frag);
        } else if (startWith(line, "END IONS")) {
            ep->parentIon.charge = charge;
            ep->order = order++;
            sort(ep->fragments.begin(), ep->fragments.end()); // fragment masses in ascending order
            expPeptides.push_back(ep);
            inIon = false;
        }

        getline(pkfile, line);
    }

} // readMgf


void readAAMasses(vector<double> &aaMasses) {
    aaMasses.assign('Z', -1.0e10);
    aaMasses['A'] = 71.03711;
    aaMasses['R'] = 156.10111;
    aaMasses['N'] = 114.04293;
    aaMasses['D'] = 115.02694;
    aaMasses['C'] = 103.00919;
    aaMasses['E'] = 129.04259;
    aaMasses['Q'] = 128.05858;
    aaMasses['G'] = 57.02146;
    aaMasses['H'] = 137.05891;
    aaMasses['I'] = 113.08406;
    aaMasses['L'] = 113.08406;
    aaMasses['K'] = 128.09496;
    aaMasses['M'] = 131.04049;
    aaMasses['F'] = 147.06841;
    aaMasses['P'] = 97.05276;
    aaMasses['S'] = 87.03203;
    aaMasses['T'] = 101.04768;
    aaMasses['W'] = 186.07931;
    aaMasses['Y'] = 163.06333;
    aaMasses['V'] = 99.06841;

    aaMasses['C'] += 57.02146; // iodoacetamide
    protonMass = 1.00728;
    hydrogenMass = 1.007825;
    oxygenMass = 15.994915;
    carbonMass = 12.000000;
    nitrogenMass = 14.003074;
    ammoniaMass = nitrogenMass + 3.0 * hydrogenMass;
    waterMass = 2.0 * hydrogenMass + oxygenMass;

} // readAAMasses


double getPeptideMass(const string &protein, const PeptidePosition &pos) {
    // Computes a peptide mass
    double mass = waterMass;
    for (unsigned i = pos.first; i <= pos.last; i++)
        mass += aaMasses[protein[i]];
    return mass;

} // getPeptideMass


double getPeptideMass(const string &peptide) {
    // Computes a peptide mass
    double mass = waterMass;
    for (unsigned i = 0; i < peptide.length(); i++)
        mass += aaMasses[peptide[i]];
    return mass;

} // getPeptideMass


void digestByTrypsin(const string &protein, unsigned nmc, vector<PeptidePosition> &pept) {
    // Digestion by trypsin with a maximum of nmc missed cleavages in a peptide. The
    // peptides are returned as a list of sequence position in the protein variable.

    // Forces size=0 in case previous positions were still in pept
    pept.resize(0);

    // Computes the position of peptides without missed cleavage
    unsigned previous = 0;
    for (unsigned i = 0; i < protein.length() - 1; i++)
        if (((protein[i] == 'K') || (protein[i] == 'R')) && (protein[i + 1] != 'P')) {
            pept.push_back(PeptidePosition(previous, i, i - previous + 1, 0));
            previous = i + 1;
        }
    pept.push_back(PeptidePosition(previous, protein.length() - 1, protein.length() - previous, 0));

    // Computes the position of peptides with up to nmc missed cleavages
    unsigned numPept = pept.size();
    for (unsigned i = 0; i < numPept - 1; i++)
        for (unsigned j = 1; (j <= nmc) && (i + j < numPept); j++)
            pept.push_back(PeptidePosition(pept[i].first, pept[i + j].last, pept[i + j].last - pept[i].first + 1, j));

} // digestByTrypsin


// Builds an index for retrieving experimental peptides when compared to theoretical masses
void buildMassIndex(vector<Mass2Pept> &massIndex, const vector<ExpPeptide *> &expPeptides) {
    Mass2Pept m2p;
    for (unsigned i = 0; i < expPeptides.size(); i++)
        if (expPeptides[i]->parentIon.charge == 0)
            for (int z = 2; z <= 3; z++) { // consider charges 2 and 3 only (90-95% of the data usually)
                m2p.mass = z * (expPeptides[i]->parentIon.moz - protonMass);
                m2p.expPeptide = expPeptides[i];
                massIndex.push_back(m2p);
            }
        else {
            m2p.mass = expPeptides[i]->parentIon.charge * (expPeptides[i]->parentIon.moz - protonMass);
            m2p.expPeptide = expPeptides[i];
            massIndex.push_back(m2p);
        }

    sort(massIndex.begin(), massIndex.end());

} // buildMassIndex


void computeTheoSpectrum(const string &peptSeq, vector<TheoreticalFragmentation> &theoSpectrum) {
    theoSpectrum.resize(9);
    for (unsigned i = 0; i < theoSpectrum.size(); i++)
        theoSpectrum[i].moz.resize(0);

    theoSpectrum[0].type = "a";
    theoSpectrum[1].type = "b";
    theoSpectrum[2].type = "b-NH3";
    theoSpectrum[3].type = "b-H2O";
    theoSpectrum[4].type = "b++";
    double m = hydrogenMass;
    bool waterLoss = false;
    bool ammoniaLoss = false;
    for (unsigned i = 0; i < peptSeq.length(); i++) {
        m += aaMasses[peptSeq[i]];
        if ((peptSeq[i] == 'S') || (peptSeq[i] == 'T'))
            waterLoss = true;
        if ((peptSeq[i] == 'N') || (peptSeq[i] == 'Q') || (peptSeq[i] == 'R'))
            ammoniaLoss = true;
        if ((i > 0) && (i < peptSeq.length() - 1)) {
            theoSpectrum[0].moz.push_back(m - carbonMass - oxygenMass);
            theoSpectrum[1].moz.push_back(m);
            theoSpectrum[2].moz.push_back(ammoniaLoss ? m - ammoniaMass : -1.0);
            theoSpectrum[3].moz.push_back(waterLoss ? m - waterMass : -1.0);
            theoSpectrum[4].moz.push_back((m + protonMass) * 0.5);
        } else
            for (unsigned i = 0; i <= 4; i++)
                theoSpectrum[i].moz.push_back(-1.0);
    }

    theoSpectrum[5].type = "y";
    theoSpectrum[6].type = "y-NH3";
    theoSpectrum[7].type = "y-H2O";
    theoSpectrum[8].type = "y++";
    m = waterMass + protonMass;
    waterLoss = ammoniaLoss = false;
    for (int i = peptSeq.length() - 1; i >= 0; i--) {
        m += aaMasses[peptSeq[i]];
        if ((peptSeq[i] == 'S') || (peptSeq[i] == 'T'))
            waterLoss = true;
        if ((peptSeq[i] == 'N') || (peptSeq[i] == 'Q') || (peptSeq[i] == 'R'))
            ammoniaLoss = true;
        if (i > 0) {
            theoSpectrum[5].moz.push_back(m);
            theoSpectrum[6].moz.push_back(ammoniaLoss ? m - ammoniaMass : -1.0);
            theoSpectrum[7].moz.push_back(waterLoss ? m - waterMass : -1.0);
            theoSpectrum[8].moz.push_back((m + protonMass) * 0.5);
        } else
            for (unsigned i = 5; i <= 8; i++)
                theoSpectrum[i].moz.push_back(-1.0);
    }

} // computeTheoSpectrum


void computeSpMatch(vector<vector<int> > &spMatch, vector<TheoreticalFragmentation> &theoSpectrum, ExpPeptide *ep) {
    spMatch.resize(theoSpectrum.size());
    for (unsigned frag = 0; frag < theoSpectrum.size(); frag++) {
        spMatch[frag].resize(0);
        for (unsigned i = 0; i < theoSpectrum[frag].moz.size(); i++) {
            int matchIndex = -1;
            if (theoSpectrum[frag].moz[i] != -1.0) {
                double theoMoz = theoSpectrum[frag].moz[i];
                if ((theoMoz >= minMass) && (theoMoz <= maxMass)) {
                    // Dichotomic search for the closest experimental moz
                    int left = 0;
                    int right = ep->fragments.size() - 1;
                    int middle;
                    while (right >= left) {
                        middle = (right + left) / 2;
                        if (ep->fragments[middle].moz >= theoMoz)
                            right = middle - 1;
                        if (ep->fragments[middle].moz <= theoMoz)
                            left = middle + 1;
                    }

                    int closest;
                    if (left - right == 2)
                        // Found the exact moz
                        closest = middle;
                    else if (right == -1)
                        // Search for a moz smaller than the smallest moz, return the smallest
                        closest = 0;
                    else if (left == int(theoSpectrum[frag].moz.size()))
                        // Search for a moz that is larger than the largest moz, return the largest
                        closest = theoSpectrum[frag].moz.size() - 1;
                    else
                        // Normal case
                        closest = (ep->fragments[left].moz - theoMoz < theoMoz - ep->fragments[right].moz) ? left
                                                                                                           : right;

                    double expMoz = ep->fragments[closest].moz;
                    if (fabs(theoMoz - expMoz) / (theoMoz + expMoz) * 2.0 <= tol)
                        // The closest experimental moz falls within tolerance tol, keep it
                        matchIndex = closest;
                }
            }
            spMatch[frag].push_back(matchIndex);
        }
    }

} // computeSpMatch


double computeSPCScore(const vector<vector<int> > &spMatch, const vector<TheoreticalFragmentation> &theoSpectrum,
                       ExpPeptide *ep, const string &peptSeq) {
    double score = 0.0;
    for (unsigned i = 0; i < spMatch.size(); i++)
        for (unsigned j = 0; j < spMatch[i].size(); j++)
            if (spMatch[i][j] != -1)
                score += 1.0;
    return score;

} // computeSPCScore


void
matchMSMSSpectra(const string &peptSeq, ExpPeptide *ep, const string &scoring, vector<PeptideMatch> &peptideMatches,
                 const string &ac, const fstream::pos_type pos, unsigned first, unsigned last) {
    // Theoretical fragmentation
    static vector<TheoreticalFragmentation> theoSpectrum;
    computeTheoSpectrum(peptSeq, theoSpectrum);

    // Spectrum match
    static vector<vector<int> > spMatch;
    computeSpMatch(spMatch, theoSpectrum, ep);

    // Score
    double score;
    if (scoring == "spc") {
        score = computeSPCScore(spMatch, theoSpectrum, ep, peptSeq);
        if (score >= 1.5 * peptSeq.length()) { // BE CAREFUL to set a threshold high enough to avoid memory full crashes
            PeptideMatch pm(ac, pos, peptSeq, score, ep, first, last);
            peptideMatches.push_back(pm);
        }
    }

} // matchMSMSSpectra


void searchFastaFile(vector<PeptideMatch> &peptideMatches, const string fasta, unsigned nmc,
                     const vector<Mass2Pept> &massIndex, const string &scoring) {
    vector<PeptidePosition> pept;

    // Processes the entire file of protein sequences
    string line, protein, ac, id;
    ifstream fastaFile(fasta.c_str());
    fstream::pos_type pos = fastaFile.tellg(), npos;
    getline(fastaFile, line);
    int v = line.find(' ');
    ac = line.substr(1, v - 1);
    id = line.substr(v + 5, line.find(' ', v + 5) - v - 5);
    while (fastaFile.good()) {
        protein = "";
        do {
            npos = fastaFile.tellg();
            getline(fastaFile, line);
            if (line[0] == '>')
                break;
            else
                protein += line;
        } while (fastaFile.good());

        // Digests and compares the peptides with experimental data
        digestByTrypsin(protein, nmc, pept);
        for (unsigned i = 0; i < pept.size(); i++) {
            double theoMass = getPeptideMass(protein, pept[i]);
            if ((theoMass >= minMass) && (theoMass <= maxMass)) {

                // Dichotomic search to find the closest experimental mass
                int left = 0;
                int right = massIndex.size() - 1;
                int middle;
                while (right >= left) {
                    middle = (right + left) / 2;
                    if (massIndex[middle].mass >= theoMass)
                        right = middle - 1;
                    if (massIndex[middle].mass <= theoMass)
                        left = middle + 1;
                }

                unsigned closest;
                if (left - right == 2)
                    // Found the exact mass
                    closest = middle;
                else if (right == -1)
                    // Search for a mass smaller than the smallest, return the smallest
                    closest = 0;
                else if (left == int(massIndex.size()))
                    // Search for a mass that is larger than the largest mass, return the largest
                    closest = massIndex.size() - 1;
                else
                    // Normal case
                    closest = (massIndex[left].mass - theoMass < theoMass - massIndex[right].mass) ? left : right;

                // Compares with close enough experimental masses
                string peptSeq = protein.substr(pept[i].first, pept[i].length);
                for (unsigned j = closest; (j >= 0) && (massIndex[j].mass >= theoMass * (1.0 - tol)); j--)
                    matchMSMSSpectra(peptSeq, massIndex[j].expPeptide, scoring, peptideMatches, ac, pos, pept[i].first,
                                     pept[i].last);
                for (unsigned j = closest + 1;
                     (j < massIndex.size()) && (massIndex[j].mass <= theoMass * (1.0 + tol)); j++)
                    matchMSMSSpectra(peptSeq, massIndex[j].expPeptide, scoring, peptideMatches, ac, pos, pept[i].first,
                                     pept[i].last);
            }
        }

        // Next protein sequence
        if (fastaFile.good()) {
            v = line.find(' ');
            ac = line.substr(1, v - 1);
            id = line.substr(v + 5, line.find(' ', v + 5) - v - 5);
            pos = npos;
        }
    }

} // searchFastaFile


void groupPeptideMatches(vector<ProteinMatch> &protMatches, const vector<PeptideMatch> &peptideMatches,
                         const string &fasta) {
    ProteinMatch pm;
    map<string, unsigned> acToProtMatch;
    for (unsigned i = 0; i < peptideMatches.size(); i++) {
        unsigned protIndex;
        if (acToProtMatch.find(peptideMatches[i].ac) == acToProtMatch.end()) {
            // First time we encounter this AC, create a new protein match
            pm.ac = peptideMatches[i].ac;
            pm.pos = peptideMatches[i].pos;
            protMatches.push_back(pm);
            protIndex = protMatches.size() - 1;
            acToProtMatch[peptideMatches[i].ac] = protIndex;
        } else
            protIndex = acToProtMatch[peptideMatches[i].ac];

        // Adds the peptide match numkber i to the list for the protein
        protMatches[protIndex].peptMatch.push_back(i);
    }

    // Compute protein coverages and scores, loads descriptions and sequences (descriptions start with \\DE=)
    ifstream fastaFile(fasta.c_str());
    vector<bool> covered;
    string sequence, line;
    map<string, double> already;
    for (unsigned i = 0; i < protMatches.size(); i++) {
        // Loads description and sequence
        fastaFile.seekg(protMatches[i].pos);
        getline(fastaFile, line);
        int v = line.find(' ');
        int w = line.find(' ', v + 5);
        protMatches[i].id = line.substr(v + 5, w - v - 5);
        v = line.find("\\DE=");
        protMatches[i].descr = line.substr(w + 5, line.length() - w - 5);
        sequence = "";
        do {
            getline(fastaFile, line);
            if (line[0] == '>')
                break;
            else
                sequence += line;
        } while (fastaFile.good());
        protMatches[i].sequence = sequence;

        // Computes coverage and score
        already.clear();
        covered.resize(sequence.length());
        for (unsigned j = 0; j < covered.size(); j++)
            covered[i] = false;
        for (unsigned j = 0; j < protMatches[i].peptMatch.size(); j++) {
            // Flags matched parts
            for (unsigned k = peptideMatches[protMatches[i].peptMatch[j]].first;
                 k <= peptideMatches[protMatches[i].peptMatch[j]].last; k++)
                covered[k] = true;
            // Distinct peptides highest scores
            if (already.find(peptideMatches[protMatches[i].peptMatch[j]].sequence) == already.end())
                already[peptideMatches[protMatches[i].peptMatch[j]].sequence] = peptideMatches[protMatches[i].peptMatch[j]].score;
            else if (already[peptideMatches[protMatches[i].peptMatch[j]].sequence] <
                     peptideMatches[protMatches[i].peptMatch[j]].score)
                already[peptideMatches[protMatches[i].peptMatch[j]].sequence] = peptideMatches[protMatches[i].peptMatch[j]].score;
        }
        // Coverage computation
        unsigned tot = 0;
        for (unsigned j = 0; j < covered.size(); j++)
            if (covered[j])
                tot++;
        protMatches[i].coverage = static_cast<double>(tot) / sequence.length();
        // Score computation
        protMatches[i].score = 0.0;
        for (map<string, double>::const_iterator si = already.begin(); si != already.end(); si++)
            protMatches[i].score += si->second;

        // Sorts peptide matches
        sort(protMatches[i].peptMatch.begin(), protMatches[i].peptMatch.end());
    }
    sort(protMatches.begin(), protMatches.end());

} // groupPeptideMatches


void outputProteinMatches(const vector<ProteinMatch> &protMatches, const vector<PeptideMatch> &peptideMatches) {
    cout << "Identified proteins:\n\n";
    for (int i = protMatches.size() - 1; i >= 0; i--) {
        cout << protMatches[i].score << "  " << protMatches[i].ac << " " << protMatches[i].id << " "
             << protMatches[i].descr << endl;
        cout << "    " << protMatches[i].peptMatch.size() << " peptides, coverage=" << protMatches[i].coverage << endl;
        for (int j = protMatches[i].peptMatch.size() - 1; j >= 0; j--)
            cout << "    " << peptideMatches[protMatches[i].peptMatch[j]].sequence << " "
                 << peptideMatches[protMatches[i].peptMatch[j]].score << " spectrum="
                 << peptideMatches[protMatches[i].peptMatch[j]].expPeptide->order << endl;
        cout << endl;
    }

} // outputProteinMatches


void usage() {
    cerr << "Usage: msmsDBSearch [options] -fasta <filename> -mgf <filename>\n\nOptions are:\n";
    cerr << "-help\n";
    cerr << "-verbose                         prints more details about the computation\n";
    cerr << "-fasta <filename>                a protein sequence sequence database in the fasta format\n";
    cerr << "-mgf <filename>                  experimental MS/MS data in the mgf format\n";
    cerr << "-nmc <int>                       maximum number of missed cleavages, default is 1\n";
    cerr << "-minmass <float>                 minimum mass considered in Daltons, default is 200\n";
    cerr << "-maxmass <float>                 maximum mass considered in Daltons, default 3500\n";
    cerr << "-delta <float>                   mass tolerance in ppm, default is 50\n";
    cerr << "-scoring <string>                scoring function name (spc)\n\n";
    exit(0);

} // usage


int main(int argc, char *argv[]) {
    string fasta, mgf;
    unsigned nmc = 1;
    string scoring = "spc";
    double delta = 800.0;

    readAAMasses(aaMasses);

    // Scans the command line
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-help") == 0)
            usage();
        else if ((strcmp(argv[i], "-fasta") == 0) && (i < argc - 1))
            fasta = string(argv[++i]);
        if ((strcmp(argv[i], "-mgf") == 0) && (i < argc - 1))
            mgf = string(argv[++i]);
        else if ((strcmp(argv[i], "-nmc") == 0) && (i < argc - 1))
            nmc = atoi(argv[++i]) >= 0 ? atoi(argv[i]) : 0;
        else if ((strcmp(argv[i], "-minmass") == 0) && (i < argc - 1))
            minMass = atof(argv[++i]);
        else if ((strcmp(argv[i], "-maxmass") == 0) && (i < argc - 1))
            maxMass = atof(argv[++i]);
        else if ((strcmp(argv[i], "-delta") == 0) && (i < argc - 1))
            delta = atof(argv[++i]);
        else if ((strcmp(argv[i], "-scoring") == 0) && (i < argc - 1))
            scoring = string(argv[++i]);
        else if (strcmp(argv[i], "-verbose") == 0)
            verbose = 1;
    }
    tol = delta / 1.0e6;

    // Reads experimental data
    cerr << "Reading...\n";
    vector<ExpPeptide *> expPeptides;
    readMgf(mgf, expPeptides);
    vector<Mass2Pept> massIndex;
    buildMassIndex(massIndex, expPeptides);

    // Searches the DB
    cerr << "Searching...\n";
    vector<PeptideMatch> peptideMatches;
    searchFastaFile(peptideMatches, fasta, nmc, massIndex, scoring);
    sort(peptideMatches.begin(), peptideMatches.end());

    // Groups peptide matches
    cerr << "Grouping...\n";
    vector<ProteinMatch> protMatches;
    groupPeptideMatches(protMatches, peptideMatches, fasta);

    // Print results
    outputProteinMatches(protMatches, peptideMatches);

    return 0;

} // main
