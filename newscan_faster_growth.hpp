extern "C" {
#include "xerrors.h"
}
#include <assert.h>
#include <fstream>
#include <iostream>
#include <oneapi/tbb/concurrent_set.h>
#include <oneapi/tbb/concurrent_unordered_map.h>
#include <stdint.h>
using namespace std;
using namespace oneapi::tbb;

pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

// struct shared via mt_parse
typedef struct {
  concurrent_unordered_map<uint64_t, concurrent_set<uint32_t>>
      *wordFreq;               // shared dictionary
  Args *arg;                   // command line input
  long start, end;             // input
  long skipped, parsed, words; // output
  FILE *parse, *last, *sa;
  uint32_t startingSeqNumber;
} mt_data;

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
  uint64_t hash = 0;
  // const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
  const uint64_t prime =
      27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
  for (size_t k = 0; k < s.size(); k++) {
    int c = (unsigned char)s[k];
    assert(c >= 0 && c < 256);
    hash = (256 * hash + c) % prime; //  add char k
  }
  return hash;
}

void *mt_parse(void *dx) {
  // extract input data
  mt_data *d = (mt_data *)dx;
  Args *arg = d->arg;
  concurrent_unordered_map<uint64_t, concurrent_set<uint32_t>> *wordFreq =
      d->wordFreq;

  if (arg->verbose > 1)
    printf("Scanning from %ld, size %ld\n", d->start, d->end - d->start);

  // open input file
  ifstream f(arg->inputFileName);
  if (!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }

  // prepare for parsing
  f.seekg(d->start); // move to the beginning of assigned region
  KR_window krw(arg->w);
  uint64_t currentWordLength = 1;
  uint64_t currentHash = Dollar;
  const uint64_t prime = 27162335252586509;
  int c; // string word = "";

  // there is some parsing to do
  uint64_t pos = d->start; // ending position+1 in text of previous word
  d->parsed = 0;
  // if(pos>0) pos+= d->skipped+ arg->w;  // or 0 for the first word
  while ((c = f.get()) != EOF) {
    if (c <= Dollar) {
      // if we are not simply compressing then we cannot accept 0,1,or 2
      cerr << "Invalid char found in input file. Exiting...\n";
      exit(1);
    }
    // word.append(1,c);
    if (c == '\n') {
      save_update_word(*arg, d->startingSeqNumber, currentHash, *wordFreq,
                       d->parse, d->last, d->sa, pos);
      currentWordLength = 1; // arg->w;
      currentHash = Dollar;  // kr_hash(krw.get_window());
      d->words++;
      d->parsed++;
      d->startingSeqNumber++;
      if (d->start + d->parsed >= d->end) {
        /* cout << "Thread " << d->start */
        /*      << " with startingSeqNumber = " << d->startingSeqNumber */
        /*      << " has finished" << endl; */
        f.close();
        return NULL;
      }
    } else {
      currentWordLength++;
      currentHash += (256 * currentHash + c) % prime;
      uint64_t hash = krw.addchar(c);
      d->parsed++;
      if (hash % arg->p == 0 && d->parsed > arg->w) {
        // end of word, save it and write its full hash to the output file
        // pos is the ending position+1 of previous word and is updated in the
        // next call
        // save_update_word(*arg,word,*wordFreq,d->parse,d->last,d->sa,pos);
        save_update_word(*arg, d->startingSeqNumber, currentHash, *wordFreq,
                         d->parse, d->last, d->sa, pos);
        currentWordLength = arg->w;
        // currentHash = hash;
        currentHash = kr_hash(krw.get_window());
        d->words++;
        if (d->start + d->parsed >= d->end + arg->w) {
          // cout << "Thread " << d->start << " with startingSeqNumber = " <<
          // d->startingSeqNumber << " has finished" << endl;
          f.close();
          return NULL;
        }
      }
    }
  }
  // end of file reached
  // virtually add w null chars at the end of the file and add the last word in
  // the dict word.append(arg->w,Dollar);
  for (int i = 0; i < arg->w; i++) {
    currentWordLength++;
    currentHash += (256 * currentHash + Dollar) % prime;
  }
  save_update_word(*arg, currentWordLength, currentHash, *wordFreq, d->parse,
                   d->last, d->sa, pos);
  // close input file and return
  f.close();
  /* cout << "Thread " << d->start */
  /*      << " with startingSeqNumber = " << d->startingSeqNumber */
  /*      << " has finished" << endl; */
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
uint32_t mt_process_file(
    Args &arg,
    concurrent_unordered_map<uint64_t, concurrent_set<uint32_t>> &wf) {
  // get input file size
  ifstream f(arg.inputFileName, std::ifstream::ate);
  if (!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + arg.inputFileName);
  }
  long size = f.tellg();

  // rewind file and store starting line positions
  f.seekg(0);
  vector<uint64_t> line_pos;
  line_pos.push_back(0);
  int c;
  while ((c = f.get()) != EOF) {
    if (c == '\n') {
      line_pos.push_back(f.tellg());
      line_pos.back()++;
    }
  }

  f.close();

  // distribute the work among threads based on line positions
  vector<uint64_t> ranges(arg.th + 1);
  for (int i = 0; i < arg.th; i++) {
    ranges[i] = i * (line_pos.size() / arg.th);
    /* cout << "Range " << i << " starts at " << ranges[i] << endl; */
  }
  ranges[arg.th] = line_pos.size() - 1;

  // prepare and execute threads
  assert(arg.th > 0);
  pthread_t t[arg.th];
  mt_data td[arg.th];
  for (int i = 0; i < arg.th; i++) {
    td[i].wordFreq = &wf;
    td[i].arg = &arg;
    // td[i].start = i*(size/arg.th); // range start
    // td[i].end = (i+1==arg.th) ? size : (i+1)*(size/arg.th); // range end
    td[i].start = line_pos[ranges[i]];
    td[i].end = i == arg.th - 1 ? size : line_pos[ranges[i + 1]];
    td[i].startingSeqNumber = ranges[i];

    /* cout << "Thread " << i << " starts at (" << ranges[i] << ") " <<
     * td[i].start */
    /*      << " and ends at (" << ranges[i + 1] << ") " << td[i].end << endl;
     */

    assert(td[i].end <= size);
    // open the 1st pass parsing file
    // td[i].parse =
    // open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,i,"wb");
    td[i].last = td[i].sa = NULL;
    xpthread_create(&t[i], NULL, &mt_parse, &td[i], __LINE__, __FILE__);
  }

  // wait for the threads to finish (in order) and close output files
  long tot_char = 0;
  uint64_t tot_word = 0;
  for (int i = 0; i < arg.th; i++) {
    xpthread_join(t[i], NULL, __LINE__, __FILE__);
    if (arg.verbose) {
      cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:";
      cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" << td[i].words
           << endl;
    }
    // close thread-specific output files
    // fclose(td[i].parse);
    // if(td[i].last) fclose(td[i].last);
    // if(td[i].sa) fclose(td[i].sa);
    // tot_word += td[i].words;
    // if(td[i].words>0) {
    //   // extra check
    //   assert(td[i].parsed>arg.w);
    //   tot_char += td[i].parsed - (i!=0? arg.w: 0); //parsed - overlapping
    // }
    // else assert(i>0); // the first thread must produce some words
    /* cout << "Thread " << i */
    /*      << " has td[i].startingSeqNumber = " << td[i].startingSeqNumber */
    /*      << endl; */
  }
  // assert(tot_char==size);
  return td[arg.th - 1].startingSeqNumber;
}
