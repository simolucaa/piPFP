#include <assert.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <oneapi/tbb/concurrent_unordered_map.h>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#ifdef GZSTREAM
#endif
extern "C" {
#include "utils.h"
#include "xerrors.h"
}

using namespace std;
using namespace oneapi::tbb;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
  string inputFileName = "";
  int w = 10;      // sliding window size and its default
  int p = 100;     // modulus for establishing stopping w-tuples
  int th = 0;      // number of helper threads
  int verbose = 0; // verbosity level
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
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
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

static void
save_update_word(uint64_t currentWordLength, uint64_t hash,
                 concurrent_unordered_map<uint64_t, uint64_t> &freq);

#ifndef NOTHREADS
#include "piPFP.hpp"
#endif

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void
save_update_word(uint64_t currentWordLength, uint64_t hash,
                 concurrent_unordered_map<uint64_t, uint64_t> &freq) {
  freq.insert({hash, currentWordLength});
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args &arg,
                      concurrent_unordered_map<uint64_t, uint64_t> &wordFreq) {
  // open a, possibly compressed, input file
  string fnam = arg.inputFileName;
  ifstream f(fnam);
  if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + fnam);
  }

  // main loop on the chars of the input file
  int c, pc = '\n', IN_HEADER = 1;
  KR_window krw(arg.w);
  uint64_t currentWordLength = 1;
  uint64_t currentHash = Dollar;
  const uint64_t prime = 27162335252586509;
  uint64_t parseWords = 0;
  while ((c = f.get()) != EOF) {
    if (pc == '\n')
      IN_HEADER = (c == '>');
    if (c != '\n' && !IN_HEADER) {
      currentWordLength++;
      currentHash += (256 * currentHash + c) % prime;
      uint64_t hash = krw.addchar(c);
      if (hash % arg.p == 0) {
        save_update_word(currentWordLength, currentHash, wordFreq);
        currentWordLength = arg.w;
        currentHash = krw.hash; // kr_hash(krw.get_window());
        parseWords++;
      }
    } else {
      if (c == '\n' && IN_HEADER) {
        currentWordLength++;
        currentHash += (256 * currentHash + c) % prime;
        uint64_t hash = krw.addchar(c);
        if (hash % arg.p == 0) {
          save_update_word(currentWordLength, currentHash, wordFreq);
          currentWordLength = arg.w;
          currentHash = krw.hash; // kr_hash(krw.get_window());
          parseWords++;
        }
      }
    }
    pc = c;
  }
  // virtually add w null chars at the end of the file and add the last word
  // in the dict word.append(arg.w,Dollar);
  for (int i = 0; i < arg.w; i++) {
    currentWordLength++;
    currentHash += (256 * currentHash + Dollar) % prime;
  }
  save_update_word(currentWordLength, currentHash, wordFreq);
  // close input and output files
  f.close();
  return parseWords;
}

void print_help(char **argv, Args &args) {
  cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
       << "\t-w W\tsliding window size, def. " << args.w << endl
       << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
#ifndef NOTHREADS
       << "\t-t M\tnumber of helper threads, def. none " << endl
#endif
       << "\t-c  \tdiscard redundant information" << endl
       << "\t-h  \tshow help and exit" << endl
       << "\t-s  \tcompute suffix array info" << endl;
#ifdef GZSTREAM
  cout << "If the input file is gzipped it is automatically extracted\n";
#endif
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
  while ((c = getopt(argc, argv, "p:w:sht:vc")) != -1) {
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
    arg.inputFileName.assign(argv[optind]);
  } else {
    cout << "Invalid number of arguments" << endl;
    print_help(argv, arg);
  }
  // check algorithm parameters
  if (arg.w < 4) {
    cout << "Windows size must be at least 4\n";
    exit(1);
  }
  if (arg.p < 10) {
    cout << "Modulus must be at leas 10\n";
    exit(1);
  }
#ifdef NOTHREADS
  if (arg.th != 0) {
    cout << "The NT version cannot use threads\n";
    exit(1);
  }
#else
  if (arg.th < 0) {
    cout << "Number of threads cannot be negative\n";
    exit(1);
  }
#endif
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
  // init sorted map counting the number of occurrences of each word
  concurrent_unordered_map<uint64_t, uint64_t> wordFreq;
  pair<uint64_t, uint64_t> totCharAndWord;

  // ------------ parsing input file
  try {
    if (arg.th == 0) {
      cout << "Parsing input file\n";
      totCharAndWord.second = process_file(arg, wordFreq);
    } else {
      cout << "Parsing input file using " << arg.th << " threads\n";
#ifdef NOTHREADS
      cerr << "Sorry, this is the no-threads executable and you requested "
           << arg.th << " threads\n";
      exit(EXIT_FAILURE);
#else
      totCharAndWord = mt_process_file(arg, wordFreq);
#endif
    }
  } catch (const std::bad_alloc &) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
  }
  // first report
  uint64_t totDWord = wordFreq.size();
  cout << "Found " << totDWord << " distinct words" << endl;
  cout << "Parsing took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";

  // -------------- second pass
  start_wc = time(NULL);
  uint64_t sumLen = 0;
  for (auto &x : wordFreq) {
    sumLen += x.second;
  }
  cout << "Sum of lenghts of dictionary words: " << sumLen << endl;
  cout << "Total number of words: " << totCharAndWord.second + 1 << endl;
  cout << "==== Elapsed time: " << difftime(time(NULL), start_main)
       << " wall clock seconds\n";
  return 0;
}
