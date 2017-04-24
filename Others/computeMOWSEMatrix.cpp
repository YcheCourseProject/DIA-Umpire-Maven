#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//  MOWSE matrix computation
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


// Computation of the MOWSE matrix from a database of protein sequences in the fasta format. Compiled with g++ and tested
// under Linux.


using namespace std;


class PeptidePosition {
public:
    PeptidePosition() : first(-1), last(-1), length(-1), nmc(0) {}

    PeptidePosition(int f, int l, int len, int n) {
        first = f;
        last = l;
        length = len;
        nmc = n;
    }

    ~PeptidePosition() {}

    int first, last, length, nmc;

}; // class PeptidePosition


// Global variables to store basic masses
vector<double> aaMasses;
double waterMass, protonMass;


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


double getPeptideMass(const string &protein, const PeptidePosition &pos) {
    // Computes a peptide mass
    double mass = waterMass;
    for (int i = pos.first; i <= pos.last; i++) {
        if ((protein[i] == 'B') || (protein[i] == 'X') || (protein[i] == 'Z'))
            return 0.0;
        mass += aaMasses[protein[i]];
    }
    return mass;

} // getPeptideMass


double getPeptideMass(const string &peptide) {
    // Computes a peptide mass
    double mass = waterMass;
    for (int i = 0; i < peptide.length(); i++) {
        if ((peptide[i] == 'B') || (peptide[i] == 'X') || (peptide[i] == 'Z'))
            return 0.0;
        mass += aaMasses[peptide[i]];
    }
    return mass;

} // getPeptideMass


void digestByTrypsin(const string &protein, int nmc, vector<PeptidePosition> &pept) {
    // Digestion by trypsin with a maximum of nmc missed cleavages in a peptide. The
    // peptides are returned as a list of sequence position in the protein variable.

    // Forces size=0 in case previous positions were still in pept
    pept.resize(0);

    // Computes the position of peptides without missed cleavage
    int previous = 0;
    for (int i = 0; i < protein.length() - 1; i++)
        if (((protein[i] == 'K') || (protein[i] == 'R')) && (protein[i + 1] != 'P')) {
            pept.push_back(PeptidePosition(previous, i, i - previous + 1, 0));
            previous = i + 1;
        }
    pept.push_back(PeptidePosition(previous, protein.length() - 1, protein.length() - previous, 0));

    // Computes the position of peptides with up to nmc missed cleavages
    int numPept = pept.size();
    for (int i = 0; i < numPept - 1; i++)
        for (int j = 1; (j <= nmc) && (i + j < numPept); j++)
            pept.push_back(PeptidePosition(pept[i].first, pept[i + j].last, pept[i + j].last - pept[i].first + 1, j));

} // digestByTrypsin


void computeMatrix(vector<vector<unsigned> > &M, const string fasta, int nmc) {
    // Digests all the proteins of a fasta file and collect the frequencies of peptide masses.
    vector<PeptidePosition> pept;
    M.resize(40);
    for (unsigned i = 0; i < M.size(); i++)
        M[i].assign(30, 0);

    // Processes the entire file of protein sequences
    string line, protein;
    ifstream fastaFile(fasta.c_str());
    getline(fastaFile, line); // Skips first fasta header
    while (fastaFile.good()) {
        protein = "";
        do {
            getline(fastaFile, line);
            if (line[0] == '>')
                break;
            else
                protein += line;
        } while (fastaFile.good());

        double proteinMass = getPeptideMass(protein, PeptidePosition(0, protein.length() - 1, protein.length(), 0));
        digestByTrypsin(protein, nmc, pept);
        unsigned j = static_cast<unsigned>(proteinMass / 10000.0);
        for (unsigned k = 0; k < pept.size(); k++) {
            double mass = getPeptideMass(protein, pept[k]) + protonMass;
            unsigned i = static_cast<unsigned>(mass / 100.0);
            if (i >= M.size())
                i = M.size() - 1;
            if (j >= M[0].size())
                j = M[0].size() - 1;

            // Counts
            M[i][j]++;
        }
    }

    // Finds the max counts per protein mass range
    vector<unsigned> maxCount;
    for (unsigned j = 0; j < M[0].size(); j++) {
        unsigned max = M[0][j];
        for (unsigned k = 0; k < M.size(); k++)
            if (M[k][j] > max)
                max = M[k][j];
        maxCount.push_back(max);
    }

    // Prints the MOWSE matrix on standard output
    cout << M.size() << "\t" << M[0].size() << endl;
    for (unsigned i = 0; i < M.size(); i++) {
        for (unsigned j = 0; j < M[0].size(); j++) {
            if (j > 0)
                cout << "\t";
            cout << static_cast<double>(M[i][j]) / maxCount[j];
        }
        cout << endl;
    }

} // computeMatrix


int main(int argc, char *argv[]) {
    string fasta;
    int nmc = 1;

    readAAMasses(aaMasses, waterMass, protonMass);

    // Scans the command line
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-help") == 0) {
            cerr << "Usage: computeMOWSEMatrix -fasta <filename> -nmc <int>\n";
            exit(0);
        } else if ((strcmp(argv[i], "-fasta") == 0) && (i < argc - 1))
            fasta = string(argv[++i]);
        else if ((strcmp(argv[i], "-nmc") == 0) && (i < argc - 1))
            nmc = atoi(argv[++i]);
    }

    // Digests everything
    vector<vector<unsigned> > M;
    computeMatrix(M, fasta, nmc);

    return 0;

} // main
