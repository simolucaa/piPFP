#include <assert.h>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>

#include <oneapi/tbb/concurrent_hash_map.h>

extern "C" {
#include "utils.h"
}

using namespace std;
using namespace oneapi::tbb;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX - 1)
typedef uint32_t word_int_t;
// Note: we can probably raise this to UINT32_MAX-1, but first
// we need to remove the limitationon of 2^32-2 on the number
// of words in the parsing (see bwtparse.c)

// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
  string data_directory = "";
  int w = 4;                  // sliding window size and its default
  int p = 10;                 // modulus for establishing stopping w-tuples
  string outputFileName = ""; // output file name
  int th = 0;                 // number of helper threads
  int verbose = 0;            // verbosity level
};

// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int *window;
  int asize;
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot; // asize^(wsize-1) mod prime

  KR_window(int w) : wsize(w) {
    asize = 256;
    asize_pot = 1;
    for (int i = 1; i < wsize; i++)
      asize_pot =
          (asize_pot * asize) % prime; // ugly linear-time power algorithm
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for (int i = 0; i < wsize; i++)
      window[i] = 0;
    // init hash value and related values
    hash = tot_char = 0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers
    hash += (prime -
             (window[k] * asize_pot) % prime); // remove window[k] contribution
    hash = (asize * hash + c) % prime;         //  add char i
    window[k] = c;
    return hash;
  }
  // debug only
  string get_window() {
    string w = "";
    int k = (tot_char - 1) % wsize;
    for (int i = k + 1; i < k + 1 + wsize; i++)
      w.append(1, window[i % wsize]);
    return w;
  }

  ~KR_window() { delete[] window; }
};
// -----------------------------------------------------------

static void save_update_word(uint64_t hash, unordered_set<uint64_t> &freq);
static void mergeThreadHashMapToCollectiveMap(
    unordered_set<uint64_t> &threadMap,
    concurrent_hash_map<uint64_t, uint32_t> &collectiveMap);

#ifndef NOTHREADS
#include "piPFP_growth.hpp"
#endif

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void save_update_word(uint64_t hash, unordered_set<uint64_t> &freq) {
  freq.insert(hash);
}

static void mergeThreadHashMapToCollectiveMap(
    unordered_set<uint64_t> &threadMap,
    concurrent_hash_map<uint64_t, uint32_t> &collectiveMap) {

  concurrent_hash_map<uint64_t, uint32_t>::accessor a;
  for (auto &x : threadMap) {
    auto inserted = collectiveMap.insert(a, {x, 1});
    if (inserted == false) {
      a->second++;
    }
    a.release();
  }
}

uint32_t process_file(Args &arg, unordered_set<uint64_t> &currentHS,
                      concurrent_hash_map<uint64_t, uint32_t> &wordFreq) {
  uint32_t n_sequences = 0;
  uint64_t n_chars_parsed = 0;
  uint64_t n_ts = 0;
  KR_window krw(arg.w);
  const uint64_t prime = 27162335252586509;
  for (auto &file : filesystem::directory_iterator(arg.data_directory)) {
    ifstream f(file.path());
    if (!f.rdbuf()->is_open()) {
      perror(__func__);
      throw new std::runtime_error("Cannot open input file " +
                                   file.path().generic_string());
    }

    // main loop on the chars of the input file
    int c, pc = '\n', IN_HEADER = 1;
    uint64_t currentHash = Dollar;
    while ((c = f.get()) != EOF) {
      if (pc == '\n')
        IN_HEADER = (c == '>');
      if (c != '\n' && !IN_HEADER) {
        n_chars_parsed++;
        currentHash += (256 * currentHash + c) % prime;
        uint64_t hash = krw.addchar(c);
        if (hash % arg.p == 0) {
          n_ts++;
          save_update_word(currentHash, currentHS);
          currentHash = krw.hash;
        }
      } else {
        if (c == '\n' && IN_HEADER && currentHash != Dollar) {
          save_update_word(currentHash, currentHS);
          currentHash = Dollar;
        }
      }
      pc = c;
    }
    n_ts++;
    save_update_word(currentHash, currentHS);
    mergeThreadHashMapToCollectiveMap(currentHS, wordFreq);
    currentHS.clear();
    krw.reset();
    n_sequences++;
    f.close();
  }

  cout << "Number of characters parsed " << n_chars_parsed + n_sequences
       << endl;
  cout << "Number of sequences parsed " << n_sequences << endl;
  cout << "Number of words " << wordFreq.size() << endl;
  cout << "Number of trigger strings " << n_ts << endl;

  return n_sequences;
}

void print_help(char **argv, Args &args) {
  cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
       << "\t-w W\tsliding window size, def. " << args.w << endl
       << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
#ifndef NOTHREADS
       << "\t-t M\tnumber of helper threads, def. none " << endl
#endif
       << "\t-o\t output name for the histogram and growth files (.hist and "
          ".growth), def. <input filename>"
       << endl;
  exit(1);
}

void parseArgs(int argc, char **argv, Args &arg) {
  int c;
  extern char *optarg;
  extern int optind;

  puts("==== Command line:");
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  puts("");

  string sarg;
  while ((c = getopt(argc, argv, "p:w:ht:vo:")) != -1) {
    switch (c) {
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 'p':
      sarg.assign(optarg);
      arg.p = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'v':
      arg.verbose++;
      break;
    case 'o':
      arg.outputFileName.assign(optarg);
      break;
    case 'h':
      print_help(argv, arg);
      exit(1);
    case '?':
      cout << "Unknown option. Use -h for help." << endl;
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1) {
    arg.data_directory.assign(argv[optind]);
  } else {
    cout << "Invalid number of arguments" << endl;
    print_help(argv, arg);
  }
  // check algorithm parameters
  if (arg.w < 4) {
    cout << "Windows size must be at least 4\n";
    exit(1);
  }
  if (arg.th < 0) {
    cout << "Number of threads cannot be negative\n";
    exit(1);
  }
}

int main(int argc, char **argv) {
  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;

  uint32_t totSeqNumber = 0;
  vector<unordered_set<uint64_t>> wordFreq(arg.th ? arg.th : 1);
  concurrent_hash_map<uint64_t, uint32_t> cumulativeWordFreq;

  // ------------ parsing input file
  try {
    if (arg.th == 0) {
      cout << "Parsing input file\n";
      totSeqNumber = process_file(arg, wordFreq[0], cumulativeWordFreq);
    } else {
      cout << "Parsing input file using " << arg.th << " threads\n";
#ifndef NOTHREADS
      totSeqNumber = mt_process_file(arg, wordFreq, cumulativeWordFreq);
#endif
    }
  } catch (const std::bad_alloc &) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
  }
  cout << "Parsing took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";

  start_wc = time(NULL);
  // fill array
  auto time1 = chrono::high_resolution_clock::now();
  uint64_t *histogram = new uint64_t[totSeqNumber + 1]();
  for (auto &x : cumulativeWordFreq) {
    histogram[x.second]++;
  }
  // print to file the histogram
  ofstream histFile(arg.outputFileName + ".hist");
  for (uint32_t i = 1; i < totSeqNumber; i++) {
    histFile << histogram[i] << endl;
  }
  histFile.close();
  auto time2 = chrono::high_resolution_clock::now();
  cout << "Time to fill histogram: "
       << chrono::duration_cast<chrono::milliseconds>(time2 - time1).count()
       << " ms\n";

  // compute pangenome growth
  ofstream growthFile(arg.outputFileName + ".growth");
  double tot = 0;
  double n_fall_m = 0;
  for (uint32_t i = 0; i < totSeqNumber + 1; i++)
    tot += histogram[i];

  double *F = new double[totSeqNumber + 1];
  for (uint32_t i = 0; i < totSeqNumber + 1; i++)
    F[i] = 0;

  for (uint32_t m = 1; m <= totSeqNumber; m++) {
    double y = 0;
    n_fall_m += log(totSeqNumber - m + 1);
    for (uint32_t i = 1; i <= totSeqNumber - m; i++) {
      F[i] += log((double)totSeqNumber - (double)m - (double)i + 1);
      y += exp(log(histogram[i]) + F[i] - n_fall_m);
    }
    if (m > 1)
      growthFile << " ";
    growthFile << (tot - y);
  }
  cout << '\n' << flush;

  cout << "==== Elapsed time: " << difftime(time(NULL), start_main)
       << " wall clock seconds\n";
  return 0;
}
