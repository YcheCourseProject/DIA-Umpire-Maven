#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


//  Primitive PMF search engine
//  Copyright (C) 2005, 2006 Jacques Colinge

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
// collect and organize data. No real OOP! Compiled with g++ and tested under Linux.

// No variable modifications are possible, cysteines are modified. Tryptic digestion only.

// Exercises: (1) Implements Profound (Zhang & Chait, Anal Chem, 2000; Zhang, et al., International Patent Application WO 00/73787 A1), RSA (Egelhofer, et al., Anal Chem, 2002), or OLAV-PMF (Magnin, et al., J Proteome Res, 2004) scores;
//            (2) Improve the p-value estimation;
//            (3) Implement variable oxidation on methionines;
//            (4) Implement a simple user interface (web-based).



using namespace std;

int verbose = 0;
enum {
    noScore = -300, nBest = 10
};
static const double reductionFactor = 0.25;


// Class used to represent peptides after digestion.
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


// Class used to represent one peak, either experimental or theoretical
class Peak {
public:
    Peak() : mass(0.0), intensity(0.0) {}

    Peak(double m, double i) {
        mass = m;
        intensity = i;
    }

    ~Peak() {}

    double mass, intensity;

}; // class Peak


bool operator<(const Peak &a, const Peak &b) {
    return a.mass < b.mass;

} // operator<(const Peak&, const Peak&)


// Class used to represent the match between one experimental peak and one theoretical peak
class MatchedPeaks {
public:
    MatchedPeaks() {
        theoPeak = Peak();
        expPeak = Peak();
    }

    MatchedPeaks(const Peak &t, const Peak &e) {
        theoPeak = t;
        expPeak = e;
    }

    ~MatchedPeaks() {}

    Peak theoPeak, expPeak;

}; // class MatchedPeaks


// Class used to represent the match of the experimental spectrum with one database entry
class DBMatch {
public:
    DBMatch() {
        ac = id = "";
        score = proteinMass = 0.0;
        pValue = 1.0;
        nbPept = nbTheoPept = 0;
    }

    DBMatch(string a, string i, double s, unsigned n, double pm, unsigned nt) {
        ac = a;
        score = s;
        pValue = 1.0;
        id = i;
        nbPept = n;
        proteinMass = pm;
        nbTheoPept = nt;
    }

    ~DBMatch() {}

    double score, pValue, proteinMass;
    string ac, id;
    unsigned nbPept, nbTheoPept;

}; // class DBMatch


bool operator<(const DBMatch &a, const DBMatch &b) {
    return a.score < b.score;

} // operator<(const DBMatch&, const DBMatch&)    


// Global variables to store basic masses and the MOWSE matrix
vector<double> aaMasses;
double waterMass, protonMass;
double minMass = 850, maxMass = 3500, delta = 50;
vector<vector<double> > M;


// Reads a MOWSE matrix computed by the program computeMOWSEMatrix
void readMOWSEMatrix(const string fName, vector<vector<double> > &M, const string &scoring) {
    if (scoring == "spc")
        return;

    ifstream mfile(fName.c_str());
    unsigned m, n;
    mfile >> m >> n;

    M.resize(m);
    for (unsigned i = 0; i < m; i++)
        M[i].resize(n);

    for (unsigned i = 0; i < M.size(); i++)
        for (unsigned j = 0; j < M[0].size(); j++)
            mfile >> M[i][j];

    // Checks for 0-values
    for (unsigned j = 0; j < M[0].size(); j++) {
        double min = 1.0;
        bool zeros = false;
        for (unsigned i = 0; i < M.size(); i++)
            if (M[i][j] == 0.0)
                zeros = true;
            else if (M[i][j] < min)
                min = M[i][j];
        for (unsigned i = 0; i < M.size(); i++)
            if (M[i][j] == 0.0)
                M[i][j] = min;
    }

    if (scoring == "pmowse") {
        // Converts into probabilities
        for (unsigned j = 0; j < M[0].size(); j++) {
            double sum = 0.0;
            for (unsigned i = 0; i < M.size(); i++)
                sum += M[i][j];
            for (unsigned i = 0; i < M.size(); i++)
                M[i][j] = log(M[i][j] / sum * M.size()); // log-odds with an equiprobable null-model
        }
    }

} // readMOWSEMatrix


void readAAMasses(vector<double> &aaMasses, double &water, double &proton) {
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
    water = 18.01056;
    proton = 1.00728;

} // readAAMasses


void readPkl(const string pkl, vector<Peak> &expSpectrum) {
    ifstream pkfile(pkl.c_str());
    string line;
    getline(pkfile, line);
    while (pkfile.good()) {
        int v = line.find(',');
        double m = atof(line.substr(0, v).c_str());
        double i = atof(line.substr(v + 1, line.find(',', v + 1) - v).c_str());
        expSpectrum.push_back(Peak(m, i));
        getline(pkfile, line);
    }

    // We assume the experimental masses to be sorted.
    // Otherwise: sort(expSpectrum.begin(), expSpectrum.end());

} // readPkl


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


double getProteinMass(const string &peptide) {
    // Computes a protein mass, accepts B, Z, X a.a.
    double mass = waterMass;
    for (unsigned i = 0; i < peptide.length(); i++)
        if ((peptide[i] == 'B') || (peptide[i] == 'X') || (peptide[i] == 'Z'))
            mass += 100.0; // average mass
        else
            mass += aaMasses[peptide[i]];
    return mass;

} // getProteinMass


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


void computeMatch(const double tol, const string &protein, const vector<PeptidePosition> &pept,
                  const vector<Peak> &expSpectrum, vector<MatchedPeaks> &matched, const double shift) {
    // Finds the experimental masses that fall within tol ppm to a theoretical mass
    matched.resize(0);
    for (unsigned i = 0; i < pept.size(); i++) {
        double theoMass = getPeptideMass(protein, pept[i]) + protonMass + shift;
        if ((theoMass >= minMass) && (theoMass <= maxMass)) {

            // Dichotomic search for the closest experimental mass
            int left = 0;
            int right = expSpectrum.size() - 1;
            int middle;
            while (right >= left) {
                middle = (right + left) / 2;
                if (expSpectrum[middle].mass >= theoMass)
                    right = middle - 1;
                if (expSpectrum[middle].mass <= theoMass)
                    left = middle + 1;
            }

            int closest;
            if (left - right == 2)
                // Found the exact mass
                closest = middle;
            else if (right == -1)
                // Search for a mass smaller than the smallest in expSpectrum, return the smallest
                closest = 0;
            else if (left == int(expSpectrum.size()))
                // Search for a mass that is larger than the largest mass in expSpectrum, return the largest
                closest = expSpectrum.size() - 1;
            else
                // Normal case
                closest = (expSpectrum[left].mass - theoMass < theoMass - expSpectrum[right].mass) ? left : right;

            double expMass = expSpectrum[closest].mass;
            if (fabs(theoMass - expMass) / (theoMass + expMass) * 2.0e6 <= tol)
                // The closest experimental mass falls within tolerance tol, keep it
                matched.push_back(MatchedPeaks(Peak(theoMass, 1.0), expSpectrum[closest]));
        }
    }

} // computeMatch


void matchMass(bool recalibration, const string &protein, const vector<PeptidePosition> &pept,
               const vector<Peak> &expSpectrum, vector<MatchedPeaks> &matched) {
    // Finds the experimental masses that fall with delta ppm of a theoretical mass
    computeMatch(delta, protein, pept, expSpectrum, matched, 0.0);

    if (recalibration && matched.size()) {
        // Determines a mass shift
        double shift = 0.0;
        for (unsigned i = 0; i < matched.size(); i++)
            shift += matched[i].expPeak.mass - matched[i].theoPeak.mass;
        shift /= matched.size();89
        computeMatch(reductionFactor * delta, protein, pept, expSpectrum, matched, shift);
    }

} // matchMass


double computeSPCScore(const string &protein, const vector<PeptidePosition> &pept,
                       const vector<Peak> &expSpectrum, vector<MatchedPeaks> &matched) {
    return static_cast<double>(matched.size());

} // computeSPCScore


double computeMOWSEScore(const string &protein, const vector<PeptidePosition> &pept,
                         const vector<Peak> &expSpectrum, vector<MatchedPeaks> &matched) {
    if (matched.size() == 0)
        return noScore;

    double proteinMass = getProteinMass(protein);
    double score = 50000.0 / proteinMass;
    unsigned j = static_cast<unsigned>(proteinMass / 10000.0);
    if (j >= M[0].size())
        j = M[0].size() - 1;
    for (unsigned k = 0; k < matched.size(); k++) {
        unsigned i = static_cast<unsigned>(matched[k].theoPeak.mass / 100);
        if (i >= M.size())
            i = M.size() - 1;
        score /= M[i][j];
    }

    if (score <= 1.0e-10)
        return noScore;
    else
        return log(score);

} // computeMOWSEScore


double computeProbMOWSEScore(const string &protein, const vector<PeptidePosition> &pept,
                             const vector<Peak> &expSpectrum, vector<MatchedPeaks> &matched) {
    if (matched.size() == 0)
        return noScore;

    double proteinMass = getProteinMass(protein);
    double score = log(50000.0 / proteinMass);
    unsigned j = static_cast<unsigned>(proteinMass / 10000.0);
    if (j >= M[0].size())
        j = M[0].size() - 1;
    for (unsigned k = 0; k < matched.size(); k++) {
        unsigned i = static_cast<unsigned>(matched[k].theoPeak.mass / 100);
        if (i >= M.size())
            i = M.size() - 1;
        score += M[i][j];
    }

    return score;

} // computeProbMOWSEScore


void searchFastaFile(vector<DBMatch> &dbMatch, const string fasta, unsigned nmc, const vector<Peak> &expSpectrum,
                     const string &scoring, bool recalibration) {
    vector<PeptidePosition> pept;
    vector<MatchedPeaks> matched;

    // Processes the entire file of protein sequences, WE ASSUME >AC ID ... as format
    string line, protein, ac, id;
    ifstream fastaFile(fasta.c_str());
    getline(fastaFile, line);
    int v = line.find(' ');
    ac = line.substr(1, v - 1);
    id = line.substr(v + 5, line.find(' ', v + 5) - v - 5);
    while (fastaFile.good()) {
        protein = "";
        do {
            getline(fastaFile, line);
            if (line[0] == '>')
                break;
            else
                protein += line;
        } while (fastaFile.good());

        // Digests and computes the score
        digestByTrypsin(protein, nmc, pept);
        matchMass(recalibration, protein, pept, expSpectrum, matched);
        double score;
        if (scoring == "mowse")
            score = computeMOWSEScore(protein, pept, expSpectrum, matched);
        else if (scoring == "pmowse")
            score = computeProbMOWSEScore(protein, pept, expSpectrum, matched);
        else if (scoring == "spc")
            score = computeSPCScore(protein, pept, expSpectrum, matched);
        else
            return;
        dbMatch.push_back(DBMatch(ac, id, score, matched.size(), getProteinMass(protein), pept.size()));

        if (fastaFile.good()) {
            v = line.find(' ');
            ac = line.substr(1, v - 1);
            id = line.substr(v + 5, line.find(' ', v + 5) - v - 5);
        }
    }

} // searchFastaFile


double gaussianTail(const double z) {
    //  Based upon algorithm 5666 for the error function, from:
    //  Hart, J.F. et al, 'Computer Approximations', Wiley 1968
    static double p[7] = {220.2068679123761,
                          221.2135961699311,
                          112.0792914978709,
                          33.91286607838300,
                          6.373962203531650,
                          0.7003830644436881,
                          0.03526249659989109};
    static double q[8] = {440.4137358247522,
                          793.8265125199484,
                          637.3336333788311,
                          296.5642487796737,
                          86.78073220294608,
                          16.06417757920695,
                          1.755667163182642,
                          0.08838834764831844};

    if (z > 37.)
        return 0.0;
    if (z < -37)
        return 1.0;

    double cutOff = 7.7071;
    double root2Pi = 2.506628274631001;
    double zabs = fabs(z);
    double expntl = exp(-.5 * zabs * zabs);
    double pdf = expntl / root2Pi;
    double tail;
    if (zabs < cutOff)
        tail = expntl *
               ((((((p[6] * zabs + p[5]) * zabs + p[4]) * zabs + p[3]) * zabs + p[2]) * zabs + p[1]) * zabs + p[0]) /
               (((((((q[7] * zabs + q[6]) * zabs + q[5]) * zabs + q[4]) * zabs + q[3]) * zabs + q[2]) * zabs + q[1]) *
                zabs + q[0]);

    else
        tail = pdf / (zabs + 1. / (zabs + 2. / (zabs + 3. / (zabs + 4. / (zabs + 0.65)))));

    return (z > 0.0) ? tail : 1.0 - tail;

} // gaussianTail


void outputResults(vector<DBMatch> &dbMatch, const string &scoring) {
    // Finds the best results
    sort(dbMatch.begin(), dbMatch.end());
    double maxScore = dbMatch.back().score;

    // Estimates p-values by fitting a distribution with all but the nBest best scores
    if ((scoring == "mowse") || (scoring == "pmowse")) {
        // Assumes a Gaussian distribution of the random scores
        double mean = 0.0, xx = 0.0;
        unsigned count = 0;
        for (unsigned i = 0; i < dbMatch.size() - nBest; i++)
            if (dbMatch[i].score != noScore) {
                register double v = dbMatch[i].score;
                mean += v;
                xx += v * v;
                count++;
            }
        mean /= static_cast<double>(count);
        xx /= static_cast<double>(count);
        if (xx > mean * mean) {
            double sd = sqrt(xx - mean * mean);
            int num = 0;
            double lastScore = maxScore;
            for (int i = dbMatch.size() - 1; (i >= 0) && ((num < nBest) || (dbMatch[i].score == lastScore)); i--)
                if (dbMatch[i].score != noScore) {
                    dbMatch[i].pValue = gaussianTail((dbMatch[i].score - mean) / sd);
                    num++;
                    lastScore = dbMatch[i].score;
                }
        }
    } else if (scoring == "spc") {
        // Assume an exponential distribution of the random scores
        double mean = 0.0;
        for (unsigned i = 0; i < dbMatch.size() - nBest; i++)
            mean += dbMatch[i].score;
        mean /= dbMatch.size();
        double lambda = 1.0 / mean;
        int num = 0;
        for (int i = dbMatch.size() - 1; (i >= 0) && ((i >= dbMatch.size() - nBest) || (dbMatch[i].score ==
                                                                                        dbMatch[dbMatch.size() -
                                                                                                nBest].score)); i--, num++)
            dbMatch[i].pValue = exp(-lambda * dbMatch[i].score);
    }

    // Prints the nBest first sequences
    int num = 0;
    double lastScore = maxScore;
    for (int i = dbMatch.size() - 1; (i >= 0) && ((num < nBest) || (dbMatch[i].score == lastScore)); i--)
        if (dbMatch[i].score != noScore) {
            cout << num + 1 << ": " << dbMatch[i].ac << "\t" << dbMatch[i].id << "\t" << dbMatch[i].score << "\t"
                 << dbMatch[i].nbPept << "\t" << dbMatch[i].pValue << "\t" << dbMatch[i].proteinMass << "\t"
                 << dbMatch[i].nbTheoPept << endl;
            lastScore = dbMatch[i].score;
            num++;
        }

    if (verbose)
        // Outputs the "random" score empirical distribution
        for (unsigned i = 0; i < dbMatch.size(); i++)
            if (dbMatch[i].score != noScore)
                cerr << i + 1 << "\t" << dbMatch[i].score << endl;

} // outputResults


void usage() {
    cerr << "Usage: pmfDBSearch [options] -fasta <filename> -pkl <filename>\n\nOptions are:\n";
    cerr << "-help\n";
    cerr << "-verbose                         prints more details about the computation\n";
    cerr << "-fasta <filename>                a protein sequence sequence database in the fasta format\n";
    cerr << "-pkl <filename>                  experimental PMF data in the pkl format\n";
    cerr << "-nmc <int>                       maximum number of missed cleavages, default is 1\n";
    cerr << "-minmass <float>                 minimum mass considered in Daltons, default is 850\n";
    cerr << "-maxmass <float>                 maximum mass considered in Daltons, default 3500\n";
    cerr << "-delta <float>                   mass tolerance in ppm, default is 50\n";
    cerr << "-scoring <string>                scoring function name (mowse|pmowse|spc)\n";
    cerr << "-recalibration                   triggers a posteriori recalibration\n\n";
    exit(0);

} // usage


int main(int argc, char *argv[]) {
    string fasta, pkl;
    unsigned nmc = 1;
    string scoring = "mowse";
    bool recalibration = false;

    readAAMasses(aaMasses, waterMass, protonMass);

    // Scans the command line
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-help") == 0)
            usage();
        else if ((strcmp(argv[i], "-fasta") == 0) && (i < argc - 1))
            fasta = string(argv[++i]);
        else if ((strcmp(argv[i], "-pkl") == 0) && (i < argc - 1))
            pkl = string(argv[++i]);
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
        else if (strcmp(argv[i], "-recalibration") == 0)
            recalibration = true;
    }

    // Reads experimental data
    vector<Peak> expSpectrum;
    readPkl(pkl, expSpectrum);

    // Searches the DB
    readMOWSEMatrix("MOWSE.txt", M, scoring);
    vector<DBMatch> dbMatch;
    searchFastaFile(dbMatch, fasta, nmc, expSpectrum, scoring, recalibration);

    // Print results
    outputResults(dbMatch, scoring);

    return 0;

} // main
