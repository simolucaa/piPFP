extern "C" {
#include "xerrors.h"
}
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

// struct shared via mt_parse
typedef struct {
  concurrent_unordered_map<uint64_t, uint64_t> *wordFreq; // shared dictionary
  Args *arg;                                              // command line input
  long skipped, parsed, words;                            // output
  long true_start, true_end, start, end;                  // input
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
  concurrent_unordered_map<uint64_t, uint64_t> *wordFreq = d->wordFreq;

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
  uint64_t currentWordLength = 0;
  uint64_t currentHash = 0;
  const uint64_t prime = 27162335252586509;
  int c; // string word = "";
  d->skipped = d->parsed = d->words = 0;
  int pc = '\n';
  int IN_HEADER = 1;
  if (d->true_start == 0) {
    currentWordLength = 1;
    currentHash = Dollar;
  } else { // reach the next breaking window
    while ((c = f.get()) != EOF) {
      if (pc == '\n')
        IN_HEADER = (c == '>');
      if (c != '\n' && !IN_HEADER) {
        d->skipped++;
        if (d->true_start + d->skipped == d->true_end + arg->w) {
          f.close();
          return NULL;
        }
        currentWordLength++;
        currentHash += (256 * currentHash + c) % prime;
        // word.append(1,c);
        uint64_t hash = krw.addchar(c);
        if (hash % arg->p == 0 && d->skipped >= arg->w)
          break;
      }
      pc = c;
    }
    if (c == EOF) {
      f.close();
      return NULL;
    } // reached EOF without finding a breaking point nothing to do
    d->parsed = arg->w;   // the kr-window is part of the next word
    d->skipped -= arg->w; // ... so w less chars have been skipped
    // word.erase(0,word.size() - arg->w);// keep only the last w chars
    currentWordLength = arg->w;
    currentHash = krw.hash; // kr_hash(krw.get_window());
  }

  // there is some parsing to do
  while ((c = f.get()) != EOF) {
    if (pc == '\n')
      IN_HEADER = (c == '>');
    if (c != '\n' && !IN_HEADER) {
      currentWordLength++;
      currentHash += (256 * currentHash + c) % prime;
      uint64_t hash = krw.addchar(c);
      d->parsed++;
      if (hash % arg->p == 0 && d->parsed > arg->w) {
        save_update_word(currentWordLength, currentHash, *wordFreq);
        currentWordLength = arg->w;
        currentHash = krw.hash; // kr_hash(krw.get_window());
        d->words++;
        if (d->true_start + d->skipped + d->parsed >= d->true_end + arg->w) {
          f.close();
          return NULL;
        }
      }
    } else {
      if (c == '\n' && IN_HEADER && d->parsed > arg->w) {
        currentWordLength++;
        currentHash += (256 * currentHash + c) % prime;
        uint64_t hash = krw.addchar(c);
        /* d->parsed++; */
        if (hash % arg->p == 0 && d->parsed > arg->w) {
          save_update_word(currentWordLength, currentHash, *wordFreq);
          currentWordLength = arg->w;
          currentHash = krw.hash; // kr_hash(krw.get_window());
          d->words++;
          if (d->true_start + d->skipped + d->parsed >= d->true_end + arg->w) {
            f.close();
            return NULL;
          }
        }
      }
    }
    pc = c;
  }
  // end of file reached
  // virtually add w null chars at the end of the file and add the last word in
  // the dict word.append(arg->w,Dollar);

  // append '\n'
  cout << "End of file" << endl;
  c = '\n';
  currentWordLength++;
  currentHash += (256 * currentHash + c) % prime;
  uint64_t hash = krw.addchar(c);
  if (hash % arg->p == 0 && d->parsed > arg->w) {
    save_update_word(currentWordLength, currentHash, *wordFreq);
    currentWordLength = arg->w;
    currentHash = krw.hash; // kr_hash(krw.get_window());
    d->words++;
  }
  for (int i = 0; i < arg->w; i++) {
    currentWordLength++;
    currentHash += (256 * currentHash + Dollar) % prime;
  }
  save_update_word(currentWordLength, currentHash, *wordFreq);
  // close input file and return
  f.close();
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
pair<uint64_t, uint64_t>
mt_process_file(Args &arg, concurrent_unordered_map<uint64_t, uint64_t> &wf) {
  assert(arg.th > 0);
  pthread_t t[arg.th];
  mt_data td[arg.th];
  // scan file for start positions and execute threads
  FILE *fp = fopen(arg.inputFileName.c_str(), "r");
  if (fp == NULL) {
    throw new std::runtime_error("Cannot open input file " + arg.inputFileName);
  }
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  rewind(fp);
  vector<size_t> th_sts(arg.th);
  vector<size_t> true_starts(arg.th);
  vector<size_t> true_ends(arg.th);
  th_sts[0] = 0;
  for (int i = 1; i < arg.th; ++i) {
    th_sts[i] = (size_t)(size / arg.th) * i;
  }
  int IN_HEADER = 1;
  size_t true_pos = 0, file_pos = 0;
  int j = 0, pc = 0, c = 0;
  // this loop scans the fasta file in order to properly divvy it up
  // for the threads, so they they don't accidently start in a ">" header.
  // As soon as a proper start and end position has been found, execute the
  // thread
  while (((c = fgetc(fp)) != EOF)) {
    if (j == arg.th)
      break;
    if (pc == '\n')
      IN_HEADER = (c == '>');
    if (file_pos == th_sts[j]) {
      if (IN_HEADER)
        th_sts[j]++;
      else {
        true_starts[j] = true_pos;
        if (j) {
          true_ends[j - 1] = true_pos;
          // prepare and execute thread j-1
          td[j - 1].wordFreq = &wf;
          td[j - 1].arg = &arg;
          td[j - 1].true_start = true_starts[j - 1];
          td[j - 1].true_end = true_ends[j - 1];
          td[j - 1].start = th_sts[j - 1]; // range start
          td[j - 1].end = (j == arg.th) ? size : th_sts[j];
          xpthread_create(&t[j - 1], NULL, &mt_parse, &td[j - 1], __LINE__,
                          __FILE__);
        }
        ++j;
        // check if previous thread spilled over computed start position of
        // next thread only possible in rare situations.
        if (j && j < arg.th && th_sts[j - 1] >= th_sts[j]) {
          th_sts[j] = file_pos + 1;
        }
      }
    }
    if (!IN_HEADER && c != '\n')
      ++true_pos;
    pc = c;
    ++file_pos;
  }
  assert(j == arg.th);
  // execute the last thread
  true_ends[j - 1] = size;
  td[j - 1].wordFreq = &wf;
  td[j - 1].arg = &arg;
  td[j - 1].true_start = true_starts[j - 1];
  td[j - 1].true_end = true_ends[j - 1];
  td[j - 1].start = th_sts[j - 1]; // range start
  td[j - 1].end = size;
  xpthread_create(&t[j - 1], NULL, &mt_parse, &td[j - 1], __LINE__, __FILE__);
  fclose(fp);

  // wait for the threads to finish (in order) and close output files
  uint64_t tot_word = 0;
  for (int i = 0; i < arg.th; i++) {
    xpthread_join(t[i], NULL, __LINE__, __FILE__);
    tot_word += td[i].words;
  }
  /* assert(tot_char==size); */
  return {size, tot_word};
  /* // get input file size */
  /* ifstream f(arg.inputFileName, std::ifstream::ate); */
  /* if (!f.is_open()) { */
  /*   perror(__func__); */
  /*   throw new std::runtime_error("Cannot open input file " +
   * arg.inputFileName); */
  /* } */
  /* long size = f.tellg(); */
  /* f.close(); */

  /* // prepare and execute threads */
  /* assert(arg.th > 0); */
  /* pthread_t t[arg.th]; */
  /* mt_data td[arg.th]; */
  /* for (int i = 0; i < arg.th; i++) { */
  /*   td[i].wordFreq = &wf; */
  /*   td[i].arg = &arg; */
  /*   td[i].start = i * (size / arg.th); // range start */
  /*   td[i].end = */
  /*       (i + 1 == arg.th) ? size : (i + 1) * (size / arg.th); // range end */

  /*   assert(td[i].end <= size); */
  /*   xpthread_create(&t[i], NULL, &mt_parse, &td[i], __LINE__, __FILE__); */
  /* } */

  /* // wait for the threads to finish (in order) and close output files */
  /* long tot_char = 0; */
  /* uint64_t tot_word = 0; */
  /* for (int i = 0; i < arg.th; i++) { */
  /*   xpthread_join(t[i], NULL, __LINE__, __FILE__); */
  /*   if (arg.verbose) { */
  /*     cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:"; */
  /*     cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" <<
   * td[i].words */
  /*          << endl; */
  /*   } */
  /*   tot_word += td[i].words; */
  /*   if (td[i].words > 0) { */
  /*     // extra check */
  /*     assert(td[i].parsed > arg.w); */
  /*     tot_char += td[i].parsed - (i != 0 ? arg.w : 0); // parsed -
   * overlapping */
  /*   } else */
  /*     assert(i > 0); // the first thread must produce some words */
  /* } */
  /* /1* assert(tot_char==size); *1/ */
  /* return {size, tot_word}; */
}
