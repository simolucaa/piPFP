#include <cstdint>
#include <filesystem>
#include <sys/types.h>
extern "C" {
#include "xerrors.h"
}
#include <assert.h>
#include <fstream>
#include <iostream>
#include <stdint.h>
using namespace std;

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>

// producer-consumer variables
static std::queue<string> dataQueue;
static std::mutex mtx;
static std::condition_variable cv;
static std::atomic<std::uint32_t> consumer_stop_counter;
static bool endOfProducer = false;

// struct shared via mt_parse
typedef struct {
  vector<ska::bytell_hash_set<uint64_t>> *wordFreq; // shared dictionary
  table_type *cumulativeWordFreq;
  Args *arg;                   // command line input
  long skipped, parsed, words; // output
  uint32_t startingSeqNumber;
  uint32_t nth_thread;
  string fileName;
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
  ska::bytell_hash_set<uint64_t> *wordFreq = &d->wordFreq->data()[d->nth_thread];
  table_type::handle_type handle = d->cumulativeWordFreq->get_handle();

  // open input file
  ifstream f(d->fileName);
  if (!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->data_directory);
  }

  // prepare for parsing
  KR_window krw(arg->w);
  const uint64_t prime = 27162335252586509;

  // main loop on the chars of the input file
  int c, pc = '\n', IN_HEADER = 1;
  uint64_t currentHash = Dollar;
  while ((c = f.get()) != EOF) {
    if (pc == '\n')
      IN_HEADER = (c == '>');
    if (c != '\n' && !IN_HEADER) {
      currentHash += (256 * currentHash + c) % prime;
      uint64_t hash = krw.addchar(c);
      if (hash % arg->p == 0) {
        save_update_word(currentHash, *wordFreq);
        currentHash = krw.hash;
      }
    } else {
      if (c == '\n' && IN_HEADER && currentHash != Dollar) {
        save_update_word(currentHash, *wordFreq);
        currentHash = Dollar;
      }
    }
    pc = c;
  }
  save_update_word(currentHash, *wordFreq);
  mergeThreadHashMapToCollectiveMap(*wordFreq, handle);
  wordFreq->clear();
  f.close();
  return NULL;
}

static uint64_t sequenceCount = 0;

void producer(const std::string &data_directory, const uint32_t numConsumers) {
  consumer_stop_counter = 0;

  uint32_t LIMIT = numConsumers;
  for (auto &file : filesystem::directory_iterator(data_directory)) {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock,
            [numConsumers] { return dataQueue.size() <= 2 * numConsumers; });
    dataQueue.emplace(string(file.path()));
    cv.notify_one();
    sequenceCount++;
  }

  endOfProducer = true;

  std::cout << "Producer finished reading " << sequenceCount << " sequences."
            << std::endl;
  while (consumer_stop_counter < numConsumers) {
    cv.notify_all();

    // optional: mitigate effects of busy wait
    std::this_thread::yield();
  }
  /* std::cout << "Producer finished." << std::endl; */
  return;
}

void consumer(const uint32_t id, mt_data arg) {
  while (true) {
    std::unique_lock<std::mutex> lock(mtx);
    while (dataQueue.empty()) {
      cv.wait(lock);

      // If we got woken up and the dataQueue is empty, that means that we exit.
      // We hold the lock after calling .wait(), so it is safe to access the
      // dataQueue.
      if (dataQueue.empty() && endOfProducer) {
        // std::cout << "consumer_stop_counter: " << consumer_stop_counter << "\n";
        consumer_stop_counter++;
        // std::cout << "Consumer " << id << " finished all the computation."
        //           << std::endl;
        return;
      }
    }
    string data = dataQueue.front();
    dataQueue.pop();
    // std::cout << "Consumer " << id << " consumed: " << data.first <<
    // std::endl;
    lock.unlock();

    arg.fileName = data;
    mt_parse(&arg);

    cv.notify_all();
  }
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
uint32_t mt_process_file(Args &arg, vector<ska::bytell_hash_set<uint64_t>> &wf,
                         table_type &cwf) {

  thread producerThread(producer, arg.data_directory, arg.th);
  vector<thread> consumers;
  consumers.reserve(arg.th);
  cout << "Parsing input file using " << arg.th << " threads\n";
  mt_data td[arg.th];
  cout << "Parsing input file\n";
  for (int i = 0; i < arg.th; i++) {
    // cout << "Creating thread " << i << endl;
    td[i].wordFreq = &wf;
    td[i].cumulativeWordFreq = &cwf;
    td[i].arg = &arg;
    td[i].nth_thread = i;

    consumers.push_back(thread(consumer, i, td[i]));
    // cout << "Thread " << i << " created" << endl;
  }

  producerThread.join();

  for (int i = 0; i < arg.th; i++) {
    consumers[i].join();
  }

  cv.notify_all();

  // cout << "All threads finished." << endl;
  return sequenceCount;
}
