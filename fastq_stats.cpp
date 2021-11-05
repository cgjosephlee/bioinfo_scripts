#include <zlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include "kseq.h"
using namespace std;
KSEQ_INIT(gzFile, gzread)

// https://github.com/rrwick/Filtlong/blob/f0cb3b34907b9fdf204fda3a36f61cf07b50f68f/src/read.cpp#L270
// https://github.com/lh3/seqtk/blob/7c04ce7898ad5909bd309c6ba3cd9c3bd0651f0e/seqtk.c#L1843
// https://lellansin.wordpress.com/2015/05/18/cc-使用-memoization-优化算法/
inline float qscore_to_quality(char qscore, float *memo) {
    int q = qscore - 33;
    if (memo[q] == -1.0)
        memo[q] = powf(10.0, -q / 10.0);
    return memo[q];
}

float get_mean_quality(char *qscores, float *memo) {
    int length = strlen(qscores);
    float sum = 0.0;
    for (int i = 0; i < length; ++i) {
        sum += qscore_to_quality(qscores[i], memo);
    }
    return -10.0f * log10f(sum / length);
}

inline void print_seq(FILE *fp, const kseq_t *seq) {
    if (seq->qual.l) {
        // fastq
        fputc('@', fp); fputs(seq->name.s, fp);
        if (seq->comment.l)
            fputc(' ', fp); fputs(seq->comment.s, fp);
        fputc('\n', fp);
        fputs(seq->seq.s, fp); fputc('\n', fp);
        fputc('+', fp); fputc('\n', fp);
        fputs(seq->qual.s, fp); fputc('\n', fp);
    } else {
        // fasta
        fputc('>', fp); fputs(seq->name.s, fp);
        if (seq->comment.l)
            fputc(' ', fp); fputs(seq->comment.s, fp);
        fputc('\n', fp);
        fputs(seq->seq.s, fp); fputc('\n', fp);
    }
}

vector<unsigned long> get_Nxx(
        unsigned long &seqNumber,
        unsigned long long &totLen,
        vector<unsigned long> &seqLengths
        ) {
    // assume sorted seqLengths
    vector<float> partitions {0.5, 0.9}; // N50, N90
    vector<double> thresholds(partitions.size());
    vector<unsigned long> Nxx(partitions.size());

    for (unsigned int i = 0; i < partitions.size(); ++i) {
        thresholds[i] = totLen * partitions[i];
    }

    unsigned long long sum = 0;
    unsigned int j = 0;
    for (unsigned int i = 0; i < seqNumber; ++i) {
        sum += seqLengths[i];
        if (sum > thresholds[j]) {
            Nxx[j] = seqLengths[i];
            ++j;
            if (j == partitions.size())
                break;
        }
    }
    return Nxx;
}

struct ReadStats {
    // length stat, count, base
    unsigned long lc500;
    unsigned long lc1000;
    unsigned long lc5000;
    unsigned long lc10000;
    unsigned long lc30000;
    unsigned long long lb500;
    unsigned long long lb1000;
    unsigned long long lb5000;
    unsigned long long lb10000;
    unsigned long long lb30000;
    // quality stat, count, base
    unsigned long qc5;
    unsigned long qc7;
    unsigned long qc10;
    unsigned long qc15;
    unsigned long qc20;
    unsigned long long qb5;
    unsigned long long qb7;
    unsigned long long qb10;
    unsigned long long qb15;
    unsigned long long qb20;
};

ReadStats get_ReadStats(
        unsigned long &seqNumber,
        unsigned long long &totLen,
        vector<unsigned long> &seqLengths,
        vector<float> &seqQuals
        ) {
    // assume original order
    bool hasQ = (seqLengths.size() != seqQuals.size()) ? false : true;

    ReadStats out;
    out.lc500 = out.lc1000 = out.lc5000 = out.lc10000 = out.lc30000 = 0;
    out.lb500 = out.lb1000 = out.lb5000 = out.lb10000 = out.lb30000 = 0;
    out.qc5 = out.qc7 = out.qc10 = out.qc15 = out.qc20 = 0;
    out.qb5 = out.qb7 = out.qb10 = out.qb15 = out.qb20 = 0;

    for (unsigned int i = 0; i < seqNumber; ++i) {
        if (seqLengths[i] >= 500)
            ++out.lc500, out.lb500 += seqLengths[i];
        if (seqLengths[i] >= 1000)
            ++out.lc1000, out.lb1000 += seqLengths[i];
        if (seqLengths[i] >= 5000)
            ++out.lc5000, out.lb5000 += seqLengths[i];
        if (seqLengths[i] >= 10000)
            ++out.lc10000, out.lb10000 += seqLengths[i];
        if (seqLengths[i] >= 30000)
            ++out.lc30000, out.lb30000 += seqLengths[i];
        if (hasQ) {
            if (seqQuals[i] >= 5.0f)
                ++out.qc5, out.qb5 += seqLengths[i];
            if (seqQuals[i] >= 7.0f)
                ++out.qc7, out.qb7 += seqLengths[i];
            if (seqQuals[i] >= 10.0f)
                ++out.qc10, out.qb10 += seqLengths[i];
            if (seqQuals[i] >= 15.0f)
                ++out.qc15, out.qb15 += seqLengths[i];
            if (seqQuals[i] >= 20.0f)
                ++out.qc20, out.qb20 += seqLengths[i];
        }
    }
    return out;
}

template <typename T> 
double median(vector<T> &vec, int size, bool is_sorted=true) {
    if (is_sorted) {
        // works for both ascending and descending
        if (size % 2 == 0)
            return (vec[size/2-1] + vec[size/2]) / 2.0;
        else 
            return vec[size/2];
    } else {
        // original vector will be modified
        // O(n)
        if (size % 2 == 0) {
            T m1, m2;
            nth_element(vec.begin(), vec.begin()+(size/2-1), vec.end());
            m1 = vec[size/2-1];
            m2 = *min_element(vec.begin()+(size/2), vec.end());
            return (m1 + m2) / 2.0;
        } else {
            nth_element(vec.begin(), vec.begin()+(size/2), vec.end());
            return vec[size/2];
        }
    }
}

void generate_output(
        FILE *fp,
        unsigned long &seqNumber,
        unsigned long long &totLen,
        vector<unsigned long> &seqLengths,
        vector<float> &seqQuals
        ) {
    bool hasQ = (seqLengths.size() != seqQuals.size()) ? false : true;

    double totLenGb = 0, avgLen = 0, medLen = 0;
    totLenGb = totLen / 1000000000.0;
    avgLen = totLen / (double)seqNumber;
    ReadStats Stat = get_ReadStats(seqNumber, totLen, seqLengths, seqQuals);

    // sort here
    sort(seqLengths.begin(), seqLengths.end(), greater<>());
    auto Nxx = get_Nxx(seqNumber, totLen, seqLengths);
    medLen = median(seqLengths, seqNumber);

    double totQual = 0, avgQual = 0, medQual = 0;
    if (hasQ) {
        totQual = accumulate(seqQuals.begin(), seqQuals.end(), 0.0);
        avgQual = totQual / (double)seqNumber;
        medQual = median(seqQuals, seqNumber, false);
    }

    fprintf(fp, "minimum length: %lu\n", seqLengths.back());
    fprintf(fp, "maximum length: %lu\n", seqLengths.front());
    fprintf(fp, "total seqs:     %lu\n", seqNumber);
    fprintf(fp, "total length:   %llu (%.2f Gb)\n", totLen, totLenGb);
    fprintf(fp, "N50:            %lu\n", Nxx[0]);
    fprintf(fp, "N90:            %lu\n", Nxx[1]);
    fprintf(fp, "avg. length:    %.2f\n", avgLen);
    fprintf(fp, "median length:  %.2f\n", medLen);
    if (hasQ) {
    fprintf(fp, "avg. quality:   %.2f\n", avgQual);
    fprintf(fp, "median quality: %.2f\n", medQual);
    }
    fprintf(fp, "\n");
    fprintf(fp, "L>=30000: %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.lc30000, 100.0*Stat.lc30000/seqNumber, Stat.lb30000, 100.0*Stat.lb30000/totLen);
    fprintf(fp, "L>=10000: %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.lc10000, 100.0*Stat.lc10000/seqNumber, Stat.lb10000, 100.0*Stat.lb10000/totLen);
    fprintf(fp, "L>=5000 : %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.lc5000, 100.0*Stat.lc5000/seqNumber, Stat.lb5000, 100.0*Stat.lb5000/totLen);
    fprintf(fp, "L>=1000 : %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.lc1000, 100.0*Stat.lc1000/seqNumber, Stat.lb1000, 100.0*Stat.lb1000/totLen);
    fprintf(fp, "L>=500  : %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.lc500, 100.0*Stat.lc500/seqNumber, Stat.lb500, 100.0*Stat.lb500/totLen);
    if (hasQ) {
    fprintf(fp, "\n");
    fprintf(fp, "Q>=20: %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.qc20, 100.0*Stat.qc20/seqNumber, Stat.qb20, 100.0*Stat.qb20/totLen);
    fprintf(fp, "Q>=15: %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.qc15, 100.0*Stat.qc15/seqNumber, Stat.qb15, 100.0*Stat.qb15/totLen);
    fprintf(fp, "Q>=10: %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.qc10, 100.0*Stat.qc10/seqNumber, Stat.qb10, 100.0*Stat.qb10/totLen);
    fprintf(fp, "Q>=7 : %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.qc7, 100.0*Stat.qc7/seqNumber, Stat.qb7, 100.0*Stat.qb7/totLen);
    fprintf(fp, "Q>=5 : %lu (%.2f%%)\t%llu (%.2f%%)\n", Stat.qc5, 100.0*Stat.qc5/seqNumber, Stat.qb5, 100.0*Stat.qb5/totLen);
    }
}

void generate_output_short(
        FILE *fp,
        unsigned long &seqNumber,
        unsigned long long &totLen,
        vector<unsigned long> &seqLengths,
        vector<float> &seqQuals
        ) {
    bool hasQ = (seqLengths.size() != seqQuals.size()) ? false : true;

    double totLenGb = 0, avgLen = 0, medLen = 0;
    // totLenGb = totLen / 1000000000.0;
    avgLen = totLen / (double)seqNumber;
    // ReadStats Stat = get_ReadStats(seqNumber, totLen, seqLengths, seqQuals);

    // sort here
    sort(seqLengths.begin(), seqLengths.end(), greater<>());
    auto Nxx = get_Nxx(seqNumber, totLen, seqLengths);
    medLen = median(seqLengths, seqNumber);

    double totQual = 0, avgQual = 0, medQual = 0;
    if (hasQ) {
        totQual = accumulate(seqQuals.begin(), seqQuals.end(), 0.0);
        avgQual = totQual / (double)seqNumber;
        medQual = median(seqQuals, seqNumber, false);
    }

    if (hasQ) {
        // total_base seq_num max min N50 N90 mean median meanQ medianQ
        fprintf(fp, "%llu\t%lu\t%lu\t%lu\t%lu\t%lu\t%.2f\t%.2f\t%.2f\t%.2f\n", totLen, seqNumber, seqLengths.front(), seqLengths.back(), Nxx[0], Nxx[1], avgLen, medLen, avgQual, medQual);
    } else {
        // total_base seq_num max min N50 N90 mean median
        fprintf(fp, "%llu\t%lu\t%lu\t%lu\t%lu\t%lu\t%.2f\t%.2f\n", totLen, seqNumber, seqLengths.front(), seqLengths.back(), Nxx[0], Nxx[1], avgLen, medLen);
    }
}

void generate_output_diff(
        FILE *fp,
        unsigned long nB, unsigned long n,
        unsigned long long totLenB, unsigned long long totLen
        ) {
    fprintf(fp, "\n");
    fprintf(fp, "%lu/%lu (%.2f%%) reads, ", n, nB, (double)n/nB*100.0);
    fprintf(fp, "%llu/%llu (%.2f%%) bases are preserved.\n", totLen, totLenB, (double)totLen/totLenB*100.0);
}

int main(int argc, char *argv[]) {
	gzFile fp;
	kseq_t *seq;
    FILE *fpDump = NULL, *fpStatBefore = NULL;
    unsigned long filterLen;
    float filterQual;
    bool ignoreQ = false, shortOutput = false, doFilter = false;
    int c;
    while ((c = getopt(argc, argv, "Al:q:B:D:s")) >= 0) {
        switch (c) {
            case 'A': ignoreQ = true; break;
            case 'l': doFilter = true; filterLen = atoi(optarg); break;
            case 'q': doFilter = true; filterQual = atof(optarg); break;
            case 'B': fpStatBefore = fopen(optarg, "w"); break;
            case 'D': fpDump = fopen(optarg, "w"); break;
            case 's': shortOutput = true; break;
        }
    }

	if (argc == optind) {
		printf("Usage: fastq_stats [options] <FASTX>\n");
        printf("Options:\n");
        printf("    -A       ignore qualities\n");
        printf("    -l NUM   minimum length\n");
        printf("    -q NUM   minimum mean quality\n");
        printf("    -B FILE  print read stats before filtering\n");
        printf("    -D FILE  dump per read stats to file\n");
        printf("    -s       short output\n");
		return 1;
	}
    // https://github.com/lh3/seqtk/blob/7c04ce7898ad5909bd309c6ba3cd9c3bd0651f0e/seqtk.c#L1410
    fp = optind < argc && strcmp(argv[optind], "-") ? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "Failed to open the input file/stream.\n");
		return 1;
	}
	seq = kseq_init(fp);

    // memoization
    float memo[100];
    for (int j = 0; j < 100; ++j)
        memo[j] = -1.0;

    // stats
	int r = 0;
    unsigned long n = 0;
    unsigned long long totLen = 0;
    vector<unsigned long> seqLengths;
    vector<float> seqQuals;
    float meanQ = 0.0;
    // stats before filtering
    unsigned long nB = 0;
    unsigned long long totLenB = 0;
    vector<unsigned long> seqLengthsB;
    vector<float> seqQualsB;
    // read seq
    if (doFilter) {
        while ((r = kseq_read(seq)) >= 0) {
            ++nB;
            totLenB += seq->seq.l;
            seqLengthsB.push_back(seq->seq.l);
            if (! ignoreQ && seq->qual.l) {
                meanQ = get_mean_quality(seq->qual.s, memo);
                seqQualsB.push_back(meanQ);
            }
            if (seq->seq.l >= filterLen && meanQ >= filterQual) {
                // pass
                print_seq(stdout, seq);
                ++n;
                totLen += seq->seq.l;
                seqLengths.push_back(seq->seq.l);
                if (! ignoreQ && seq->qual.l) {
                    seqQuals.push_back(meanQ);
                }
            }
            if (fpDump)
                fprintf(fpDump, "%s\t%lu\t%.6f\n", seq->name.s, seq->seq.l, meanQ);
        }
    } else {
        while ((r = kseq_read(seq)) >= 0) {
            ++n;
            totLen += seq->seq.l;
            seqLengths.push_back(seq->seq.l);
            if (! ignoreQ && seq->qual.l) {
                meanQ = get_mean_quality(seq->qual.s, memo);
                seqQuals.push_back(meanQ);
            }
            if (fpDump)
                fprintf(fpDump, "%s\t%lu\t%.6f\n", seq->name.s, seq->seq.l, meanQ);
        }
    }

	if (r != -1) fprintf(stderr, "ERROR: malformated FASTX\n");
	kseq_destroy(seq);
	gzclose(fp);
    if (fpDump) fclose(fpDump);

    // output stats
    // vectors might be modified and no longer preserve original order
    if (doFilter && fpStatBefore) {
        generate_output(fpStatBefore, nB, totLenB, seqLengthsB, seqQualsB);
        fclose(fpStatBefore);
    }

    if (shortOutput)
        generate_output_short(stderr, n, totLen, seqLengths, seqQuals);
    else
        generate_output(stderr, n, totLen, seqLengths, seqQuals);

    if (doFilter)
        generate_output_diff(stderr, nB, n, totLenB, totLen);

	return 0;
}
