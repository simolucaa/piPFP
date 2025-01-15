/* ******************************************************************************
 * newscan.cpp
 *
 * parsing algorithm for bwt construction of repetitive sequences based
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](https://doi.org/10.4230/LIPIcs.WABI.2018.2)
 *
 * Usage:
 *   newscan.x wsize modulus file
 *
 * Unless the parameter -c (compression rather than BWT construction,
 * see "Compression mode" below) is used the input file cannot contain
 * the characters 0x0, 0x1, 0x2 which are used internally.
 *
 * Since the i-th thread accesses the i-th segment of the input file
 * random access (fseek) must be possible. For gzipped inputs use
 * cnewscan.x which doesn't use threads but automatically extracts the
 * content from a gzipped input using the lz library.
 *
 * The parameters wsize and modulus are used to define the prefix free parsing
 * using KR-fingerprints (see paper)
 *
 *
 * *** BWT construction ***
 *
 * The algorithm computes the prefix free parsing of
 *     T = (0x2)file_content(0x2)^wsize
 * in a dictionary of words D and a parsing P of the input T in terms of
 * the dictionary words. Note that consecutive words in the parsing overlap
 * by wsize chars and that the words contains also the occurences of 0x2
 *
 * Let d denote the number of words in D and p the number of phrases in
 * the parsing P
 *
 * newscan.x outputs the following files:
 *
 * file.dict
 * the dictionary words in lexicographic order with a 0x1 at the end of
 * each word and a 0x0 at the end of the file. Size: |D| + d + 1 where
 * |D| is the sum of the word lengths
 *
 * file.occ
 * the number of occurrences of each word in lexicographic order.
 * We assume the number of occurrences of each word is at most 2^32-1
 * so the size is 4d bytes
 *
 * file.parse
 * the parse P with each word identified with its 1-based lexicographic
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes
 *
 * file.last
 * the character in position w+1 from the end for each word in the
 * parsing. Since there is a size-w overlapping between consecutive
 * words file.last[i] is also the character immediately preceeding the
 * word i+1 in the parsing
 * Size: p bytes
 *
 * file.sai (if option -s is given on the command line)
 * containing the ending position +1 of each parsed word in the original
 * text written using IBYTES bytes for each entry (IBYTES defined in utils.h)
 * Size: p*IBYTES
 *
 * The output of newscan.x must be processed by bwtparse, which invoked as
 *
 *    bwtparse file
 *
 * computes the BWT of file.parse and produces file.ilist of size 4p+4 bytes
 * containing, for each dictionary word in lexicographic order, the list
 * of BWT positions where that word appears (ie i\in ilist(w) <=> BWT[i]=w).
 * There is also an entry for the EOF word which is not in the dictionary
 * but is assumed to be the smallest word.
 *
 * In addition, bwtparse permutes file.last according to
 * the BWT permutation and generates file.bwlast such that file.bwlast[i] is
 * the last char of P[SA[i]-2] (if SA[i]==0 , BWT[i]=0=EOF and file.bwlast[i]=0,
 * if SA[i]==1, BWT[i]=P[0] and file.bwlast[i] is taken from P[n-1], the last
 * word in the parsing).
 *
 * If the option -s is given to bwtparse, it permutes file.sai according
 * to the BWT permutation and generate file.bwsai using again IBYTES
 * per entry.  file.bwsai[i] is the ending position+1 of BWT[i] in the
 * original text
 *
 * The output of bwtparse (the files .ilist .bwlast) together with the
 * dictionary itself (file .dict) and the number of occurrences
 * of each word (file .occ) are used to compute the final BWT by the
 * pfbwt algorithm.
 *
 * As an additional check to the correctness of the parsing, it is
 * possible to reconstruct the original file from the files .dict
 * and .parse using the unparse tool.
 *
 *
 *  *** Compression mode ***
 *
 * If the -c option is used, the parsing is computed for compression
 * purposes rather than for building the BWT. In this case the redundant
 * information (phrases overlaps and 0x2's) is not written to the output files.
 *
 * In addition, the input can contain also the characters 0x0, 0x1, 0x2
 * (ie can be any input file). The program computes a quasi prefix-free
 * parsing (with no overlaps):
 *
 *   T = w_0 w_1 w_2 ... w_{p-1}
 *
 * where each word w_i, except the last one, ends with a lenght-w suffix s_i
 * such that KR(s_i) mod p = 0 and s_i is the only lenght-w substring of
 * w_i with that property, with the possible exception of the lenght-w
 * prefix of w_0.
 *
 * In Compression mode newscan.x outputs the following files:
 *
 * file.dicz
 * containing the concatenation of the (distinct) dictionary words in
 * lexicographic order.
 * Size: |D| where |D| is the sum of the word lengths
 *
 * file.dicz.len
 * containing the lenght in bytes of the dictionary words again in
 * lexicographic order. Each lenght is represented by a 32 bit int.
 * Size: 4d where d is the number of distinct dictionary words.
 *
 * file.parse
 * containing the parse P with each word identified with its 1-based
 * lexicographic rank (ie its position in D). We assume the number of distinct
 * words is at most 2^32-1, so the size is 4p bytes.
 *
 * From the above three files it is possible to recover the original input
 * using the unparsz tool.
 *
 */
#include <assert.h>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/concurrent_hash_map.h>
#include <oneapi/tbb/parallel_for.h>

#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
}

using namespace std;
using namespace oneapi::tbb;
// using namespace __gnu_cxx;

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

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  // string str;
  occ_int_t len;
  occ_int_t occ;
  // word_int_t rank=0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
  string inputFileName = "";
  int w = 10;                 // sliding window size and its default
  int p = 100;                // modulus for establishing stopping w-tuples
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
save_update_word(Args &arg, uint64_t currenWordLength, uint64_t hash,
                 unordered_map<uint64_t, pair<uint32_t, uint32_t>> &freq);

// #ifndef NOTHREADS
#include "newscan_faster_growth2.hpp"
// #endif

// // compute 64-bit KR hash of a string
// // to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
// uint64_t kr_hash(string s) {
//     uint64_t hash = 0;
//     //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
//     const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 +
//     2**47 + 2**13) for(size_t k=0;k<s.size();k++) {
//       int c = (unsigned char) s[k];
//       assert(c>=0 && c< 256);
//       hash = (256*hash + c) % prime;    //  add char k
//     }
//     return hash;
// }

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void
save_update_word(Args &arg, uint64_t currentSequence, uint64_t hash,
                 unordered_map<uint64_t, pair<uint32_t, uint32_t>> &freq) {
  // size_t minsize = arg.w;
  // assert(pos==0 || currentWordLength > minsize);
  // if(currentWordLength <= minsize) return;

  // freq[hash] = currentWordLength;
  // freq.insert({hash, currentWordLength});
  /* #ifndef NOTHREADS */
  /*   xpthread_mutex_lock(&map_mutex, __LINE__, __FILE__); */
  /* #endif */
  auto inserted = freq.insert({hash, {currentSequence, 1}});
  if (inserted.second == false) {
    if (inserted.first->second.first != currentSequence) {
      inserted.first->second.first = currentSequence;
      inserted.first->second.second++;
    }
  }
  /* #ifndef NOTHREADS */
  /*   xpthread_mutex_unlock(&map_mutex, __LINE__, __FILE__); */
  /* #endif */
  // currentWordLength = minsize;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
/* uint64_t process_file( */
/*     Args &arg, */
/*     concurrent_unordered_map<uint64_t, concurrent_set<uint32_t>> &wordFreq) {
 */
/*   // open a, possibly compressed, input file */
/*   string fnam = arg.inputFileName; */
/* #ifdef GZSTREAM */
/*   igzstream f(fnam.c_str()); */
/* #else */
/*   ifstream f(fnam); */
/* #endif */
/*   if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams */
/*     perror(__func__); */
/*     throw new std::runtime_error("Cannot open input file " + fnam); */
/*   } */

/*   // open the 1st pass parsing file */
/*   FILE *g = NULL; // open_aux_file(arg.inputFileName.c_str(),EXTPARS0,"wb");
 */
/*   FILE *sa_file = NULL, *last_file = NULL; */
/*   // if(!arg.compress) { */
/*   //   // open output file containing the char at position -(w+1) of each
 * word */
/*   //   last_file = open_aux_file(arg.inputFileName.c_str(),EXTLST,"wb"); */
/*   //   // if requested open file containing the ending position+1 of each
 * word */
/*   //   if(arg.SAinfo) */
/*   //     sa_file = open_aux_file(arg.inputFileName.c_str(),EXTSAI,"wb"); */
/*   // } */

/*   // main loop on the chars of the input file */
/*   int c; */
/*   uint64_t pos = 0; // ending position +1 of previous word in the original
 * text, */
/*                     // used for computing sa_info */
/*   assert(IBYTES <= */
/*          sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
 */
/*   // init first word in the parsing with a Dollar char unless we are just */
/*   // compressing string word(""); if(!arg.compress) word.append(1,Dollar);
 * init */
/*   // empty KR window: constructor only needs window size */
/*   KR_window krw(arg.w); */
/*   uint64_t currentWordLength = 1; */
/*   uint64_t currentHash = Dollar; */
/*   const uint64_t prime = 27162335252586509; */
/*   uint64_t parseWords = 0; */
/*   uint32_t seqNumber = 0; */
/*   while ((c = f.get()) != EOF) { */
/*     // if(c<=Dollar && !arg.compress) { */
/*     //   // if we are not simply compressing then we cannot accept 0,1,or 2
 */
/*     //   cerr << "Invalid char found in input file. Exiting...\n"; exit(1);
 */
/*     // } */
/*     if (c == '\n') { */
/*       save_update_word(arg, seqNumber, currentHash, wordFreq, g, last_file,
 */
/*                        sa_file, pos); */
/*       currentWordLength = 1; // arg.w; */
/*       currentHash = Dollar;  // kr_hash(krw.get_window()); */
/*       parseWords++; */
/*       seqNumber++; */
/*     } else { */
/*       currentWordLength++; */
/*       currentHash += (256 * currentHash + c) % prime; */
/*       // word.append(1,c); */
/*       uint64_t hash = krw.addchar(c); */
/*       if (hash % arg.p == 0) { */
/*         // end of word, save it and write its full hash to the output file */
/*         // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" <<
 */
/*         // krw.get_window() << ">" <<  endl; */
/*         save_update_word(arg, seqNumber, currentHash, wordFreq, g, last_file,
 */
/*                          sa_file, pos); */
/*         currentWordLength = arg.w; */
/*         currentHash = kr_hash(krw.get_window()); */
/*         parseWords++; */
/*       } */
/*     } */
/*   } */
/*   // virtually add w null chars at the end of the file and add the last word
 * in */
/*   // the dict word.append(arg.w,Dollar); */
/*   for (int i = 0; i < arg.w; i++) { */
/*     currentWordLength++; */
/*     currentHash += (256 * currentHash + Dollar) % prime; */
/*   } */
/*   save_update_word(arg, seqNumber, currentHash, wordFreq, g, last_file,
 * sa_file, */
/*                    pos); */
/*   f.close(); */
/*   return seqNumber; */
/* } */

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b) { return *a <= *b; }

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
  // if (arg.p < 10) {
  //   cout << "Modulus must be at leas 10\n";
  //   exit(1);
  // }
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

class ProcessWordFreq {

  Args &arg;
  vector<unordered_map<uint64_t, pair<uint32_t, uint32_t>>> &wordFreq;
  concurrent_hash_map<uint64_t, pair<uint32_t, uint32_t>> &cumulativeWordFreq;

public:
  void operator()(const blocked_range<uint32_t> &r) const {
    concurrent_hash_map<uint64_t, pair<uint32_t, uint32_t>>::accessor a;
    for (uint32_t i = r.begin(); i != r.end(); i++) {
      for (auto &y : wordFreq[i]) {
        auto inserted = cumulativeWordFreq.insert(a, y);
        if (inserted == false) {
          a->second.second += y.second.second;
        }
        a.release();
      }
      // delete hash map
      std::unordered_map<uint64_t, pair<uint32_t, uint32_t>>().swap(
          wordFreq[i]);
    }
  }

  ProcessWordFreq(
      Args &arg,
      vector<unordered_map<uint64_t, pair<uint32_t, uint32_t>>> &wordFreq,
      concurrent_hash_map<uint64_t, pair<uint32_t, uint32_t>>
          &cumulativeWordFreq)
      : arg(arg), wordFreq(wordFreq), cumulativeWordFreq(cumulativeWordFreq) {}
};

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
  vector<unordered_map<uint64_t, pair<uint32_t, uint32_t>>> wordFreq(arg.th);
  uint32_t totSeqNumber = 0;
  uint64_t totChar;

  // ------------ parsing input file
  try {
    if (arg.th == 0) {
      cout << "Parsing input file\n";
      /* totSeqNumber = process_file(arg, wordFreq); */
    } else {
      cout << "Parsing input file using " << arg.th << " threads\n";
#ifdef NOTHREADS
      cerr << "Sorry, this is the no-threads executable and you requested "
           << arg.th << " threads\n";
      exit(EXIT_FAILURE);
#else
      totSeqNumber = mt_process_file(arg, wordFreq);
#endif
    }
  } catch (const std::bad_alloc &) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
  }
  // first report
  uint64_t totDWord = wordFreq.size();
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" << endl;
  cout << "Parsing took: " << difftime(time(NULL), start_wc)
       << " wall clock seconds\n";
  // check # distinct words
  if (totDWord > MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinc words (" << totDWord
         << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }

  // -------------- second pass
  start_wc = time(NULL);
  // create array of dictionary words
  // vector<const string *> dictArray;
  // dictArray.reserve(totDWord);

  // merge all hash maps into one
  auto time1 = chrono::high_resolution_clock::now();
  concurrent_hash_map<uint64_t, pair<uint32_t, uint32_t>> cumulativeWordFreq;
  concurrent_hash_map<uint64_t, pair<uint32_t, uint32_t>>::accessor a;
  parallel_for(blocked_range<uint32_t>(0, arg.th),
               ProcessWordFreq(arg, wordFreq, cumulativeWordFreq));

  /* for (auto &x : wordFreq) { */
  /*   for (auto &y : x) { */
  /*     auto inserted = cumulativeWordFreq.insert(y); */
  /*     if (inserted.second == false) { */
  /*       inserted.first->second.second += y.second.second; */
  /*     } */
  /*   } */
  /*   // delete hash map */
  /*   std::unordered_map<uint64_t, pair<uint32_t, uint32_t>>().swap(x); */
  /* } */
  // delete vector of hash maps
  std::vector<unordered_map<uint64_t, pair<uint32_t, uint32_t>>>().swap(
      wordFreq);
  auto time2 = chrono::high_resolution_clock::now();
  cout << "Time to merge: "
       << chrono::duration_cast<chrono::milliseconds>(time2 - time1).count()
       << " ms\n";

  // fill array
  time1 = chrono::high_resolution_clock::now();
  uint64_t sumLen = 0;
  // uint64_t totWord = 0;
  uint64_t *histogram = new uint64_t[totSeqNumber + 1]();
  for (auto &x : cumulativeWordFreq) {
    histogram[x.second.second]++;
  }
  // print to file the histogram
  ofstream histFile(arg.outputFileName + ".hist");
  for (uint32_t i = 1; i < totSeqNumber; i++) {
    histFile << histogram[i] << endl;
  }
  histFile.close();
  time2 = chrono::high_resolution_clock::now();
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

  // assert(dictArray.size()==totDWord);
  // cout << "Sum of lenghts of dictionary words: " << sumLen << endl;
  // cout << "Total number of words: " << totCharAndWord.second + 1 << endl;
  cout << "==== Elapsed time: " << difftime(time(NULL), start_main)
       << " wall clock seconds\n";
  return 0;
}

// // sort dictionary
//   sort(dictArray.begin(), dictArray.end(),pstringCompare);
//   // write plain dictionary and occ file, also compute rank for each hash
//   cout << "Writing plain dictionary and occ file\n";
//   writeDictOcc(arg, wordFreq, dictArray);
//   dictArray.clear(); // reclaim memory
//   cout << "Dictionary construction took: " <<
//   difftime(time(NULL),start_wc)
//   << " wall clock seconds\n";

//   // remap parse file
//   start_wc = time(NULL);
//   cout << "Generating remapped parse file\n";
//   remapParse(arg, wordFreq);
//   cout << "Remapping parse file took: " << difftime(time(NULL),start_wc)
//   << " wall clock seconds\n";
