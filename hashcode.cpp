
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <random>
#include <sstream>
#include <cassert>

#include "roaring.hh"
#include "roaring.c"
#include "hungarian.h"
#include "hungarian.c"

using namespace std;

class Problem {
 public:
  static Problem read(string filename);
  void readSolution(string filename);
  void writeSolution(string filename) const;

  int nbImages() const { return imageTags_.size(); }
  int nbVerticalImages() const;
  int nbTags() const { return nbTags_; }
  int nbTagReferences() const { return nbTagReferences_; }
  int nbSlides() const { return solution_.size(); }
  int objective() const;

  void initNaive();
  void initGreedy();
  void reoptimizeLocal();
  void reoptimizeAssignment();

  void report() const;
  void check() const;

 private:
  Roaring tags(pair<int, int> slide) const;
  int distance(pair<int, int> s1, pair<int, int> s2) const;
  int distance(const Roaring &tag1, const Roaring &tag2) const;

  int bestNextImageAny() const;
  int bestNextImageVertical() const;
  int randomNextImageAny() const;
  int randomNextImageVertical() const;

  void try2Opt();
  void tryExchange();
  void tryVerticalExchange();
  void optimizeVerticalAssignment();

 private:
  Problem()
    : nbTags_(0)
    , nbTagReferences_(0)
    , rgen_(random_device()())
    {}

  int nbTags_;
  int nbTagReferences_;
  vector<Roaring> imageTags_;
  vector<bool> imageVertical_;

  vector<pair<int, int> > solution_;
  unordered_set<int> solutionImages_;
  mutable mt19937 rgen_;
};

Problem Problem::read(string filename) {
  ifstream s(filename);
  int nbImgs;
  s >> nbImgs;

  Problem pb;
  unordered_map<string, int> tagIds;
  pb.imageTags_.resize(nbImgs);
  pb.imageVertical_.resize(nbImgs);
  for (int i = 0; i < nbImgs; ++i) {
    string orient;
    int nbImgTags;
    s >> orient;
    s >> nbImgTags;
    pb.imageVertical_[i] = orient == "V";
    for (int j = 0; j < nbImgTags; ++j) {
      string tag;
      s >> tag;
      auto insPair = tagIds.emplace(tag, pb.nbTags_);
      if (insPair.second)
        ++pb.nbTags_;
      int tagId = insPair.first->second;
      pb.imageTags_[i].add(tagId);
    }
    pb.imageTags_[i].runOptimize();
    pb.nbTagReferences_ += pb.imageTags_[i].cardinality();
  }

  return pb;
}

void Problem::readSolution(string filename) {
  solution_.clear();
  solutionImages_.clear();
  ifstream s(filename);
  string l;
  getline(s, l);
  stringstream ss1(l);
  int nbSlides;
  ss1 >> nbSlides;
  for (int i = 0; i < nbSlides; ++i) {
    getline(s, l);
    stringstream ss2(l);
    int f = -1;
    int s = -1;
    ss2 >> f;
    if (ss2.good())
      ss2 >> s;
    assert(f >= 0);
    solution_.emplace_back(f, s);
  }
}

void Problem::writeSolution(string filename) const {
  ofstream s(filename);
  s << solution_.size() << endl;
  for (pair<int, int> slide : solution_) {
    s << slide.first;
    if (slide.second >= 0)
      s << " " << slide.second;
    s << endl;
  }
}

int Problem::nbVerticalImages() const {
  int nb = 0;
  for (bool v : imageVertical_) {
    nb += v;
  }
  return nb;
}

void Problem::report() const {
  cout << nbImages() << " images" << endl;
  cout << nbVerticalImages() << " vertical images" << endl;
  cout << nbTags() << " tags" << endl;
  cout << nbTagReferences() << " tag references" << endl;
}

Roaring Problem::tags(pair<int, int> slide) const {
  if (slide.second >= 0)
    return imageTags_[slide.first] | imageTags_[slide.second];
  else
    return imageTags_[slide.first];
}

int Problem::distance(pair<int, int> s1, pair<int, int> s2) const {
  return distance(tags(s1), tags(s2));
}

int Problem::distance(const Roaring &tags1, const Roaring &tags2) const {
  int l1 = tags1.and_cardinality(tags2);
  int l2 = tags1.andnot_cardinality(tags2);
  int l3 = tags2.andnot_cardinality(tags1);
  return min(l1, min(l2, l3));
}

int Problem::objective() const {
  int obj = 0;
  for (int i = 0; i + 1 < (int) solution_.size(); ++i) {
    obj += distance(solution_[i], solution_[i+1]);
  }
  return obj;
}

void Problem::initNaive() {
  solution_.clear();
  solutionImages_.clear();
  int lastVSlide = -1;
  for (int i = 0; i < nbImages(); ++i) {
    if (imageVertical_[i]) {
      if (lastVSlide < 0) {
        lastVSlide = solution_.size();
        solution_.emplace_back(i, -1);
      }
      else {
        solution_[lastVSlide].second = i;
        lastVSlide = -1;
      }
    }
    else {
      solution_.emplace_back(i, -1);
    }
  }
  for (int i = 0; i < nbImages(); ++i) {
    solutionImages_.insert(i);
  }
}

void Problem::initGreedy() {
  solution_.clear();
  solutionImages_.clear();
  while (solutionImages_.size() < nbImages()) {
    if (solution_.empty()) {
      solution_.emplace_back(randomNextImageAny(), -1);
    }
    else if (imageVertical_[solution_.back().first] && solution_.back().second < 0) {
      if (solution_.size() > 1)
        solution_.back().second = bestNextImageVertical();
      else
        solution_.back().second = randomNextImageVertical();
    }
    else {
      solution_.emplace_back(bestNextImageAny(), -1);
    }
    solutionImages_.insert(solution_.back().first);
    if (solution_.back().second >= 0)
      solutionImages_.insert(solution_.back().second);
    if (solutionImages_.size() % 1000 == 0)
      cout << solutionImages_.size() << " images in the solution" << endl;
  }
}

int Problem::bestNextImageAny() const {
  int bestObj = -1;
  int best = -1;
  Roaring prevTags = tags(solution_.back());
  for (int i = 0; i < nbImages(); ++i) {
    if (solutionImages_.count(i)) continue;
    int obj = 0;
    if (!solution_.empty()) {
      obj = distance(prevTags, imageTags_[i]);
    }
    if (obj > bestObj) {
      bestObj = obj;
      best = i;
    }
  }
  return best;
}

int Problem::bestNextImageVertical() const {
  int bestObj = -1;
  int best = -1;
  pair<int, int> prevPair = solution_[solution_.size() - 1];
  Roaring prevTags = tags(prevPair);
  for (int i = 0; i < nbImages(); ++i) {
    if (solutionImages_.count(i)) continue;
    if (!imageVertical_[i]) continue;
    pair<int, int> nextPair = make_pair(solution_.back().first, i);
    int obj = distance(prevTags, tags(nextPair));
    if (obj > bestObj) {
      bestObj = obj;
      best = i;
    }
  }
  return best;
}

int Problem::randomNextImageAny() const {
  vector<int> imgs;
  for (int i = 0; i < nbImages(); ++i)
    imgs.push_back(i);
  shuffle(imgs.begin(), imgs.end(), rgen_);
  for (int i : imgs) {
    if (!solutionImages_.count(i))
      return i;
  }
  return -1;
}

void Problem::optimizeVerticalAssignment() {
  for (pair<int, int> &slide : solution_) {
    if (slide.second >= 0 && bernoulli_distribution()(rgen_))
      swap(slide.first, slide.second);
  }

  //cout << "Creating assignment problem" << endl;
  int subrangeSize = 20000;
  int start = uniform_int_distribution<int>(0, nbSlides() - 1)(rgen_);
  vector<int> allSlides;
  for (int i = 0; i < subrangeSize && i < nbSlides(); ++i)
    allSlides.push_back( (start + i) % nbSlides());
  shuffle(allSlides.begin(), allSlides.end(), rgen_);

  const int maxSlides = 5000;
  set<int> usedSet;
  for (int s : allSlides) {
    if (usedSet.size() > maxSlides) break;
    if (usedSet.count(s-1)) continue;
    if (usedSet.count(s+1)) continue;
    usedSet.insert(s);
  }

  vector<int> used(usedSet.begin(), usedSet.end());
  int **costMat;
  costMat = new int*[used.size()];
  for (int i = 0; i < used.size(); ++i) {
    costMat[i] = new int[used.size()];
  }
  // Solve an assignment problem
  for (int i = 0; i < used.size(); ++i) {
    int s1 = used[i];
    for (int j = 0; j < used.size(); ++j) {
      int s2 = used[j];
      int cost = 0;
      pair<int, int> slide (solution_[s2].first, solution_[s1].second);
      if (s2 > 0) cost += distance(slide, solution_[s2-1]);
      if (s2 + 1 < solution_.size()) cost += distance(slide, solution_[s2+1]);
      costMat[i][j] = cost;
    }
  }

  //cout << "Solving assignment problem" << endl;
  hungarian_problem_t p;
  hungarian_init(&p, costMat, used.size(), used.size(), HUNGARIAN_MODE_MAXIMIZE_UTIL);
  hungarian_solve(&p);

  vector<pair<int, int> > newSolution = solution_;
  for (int i = 0; i < used.size(); ++i) {
    int dest = -1;
    int sumTot = 0;
    for (int j = 0; j < used.size(); ++j) {
      sumTot += p.assignment[i][j];
      if (p.assignment[i][j])
        dest = j;
    }
    newSolution[used[dest]].second = solution_[used[i]].second;
    assert (sumTot == 1);
  }
  solution_ = newSolution;

  hungarian_free(&p);
  for (int i = 0; i < used.size(); ++i) {
    delete[] costMat[i];
  }
  delete[] costMat;
  //cout << "Solved assignment problem" << endl;
}

int Problem::randomNextImageVertical() const {
  vector<int> imgs;
  for (int i = 0; i < nbImages(); ++i)
    imgs.push_back(i);
  shuffle(imgs.begin(), imgs.end(), rgen_);
  for (int i : imgs) {
    if (imageVertical_[i] && !solutionImages_.count(i))
      return i;
  }
  return -1;
}

void Problem::try2Opt() {
  int e1 = uniform_int_distribution<int>(0, nbSlides()-2)(rgen_);
  //int e2 = uniform_int_distribution<int>(0, nbSlides()-2)(rgen_);
  int offset = uniform_int_distribution<int>(1, 100)(rgen_);
  int e2 = (e1 + offset) % (nbSlides() - 1);
  if (e1 == e2) return;
  assert (e1 + 1 < (int) solution_.size());
  assert (e2 + 1 < (int) solution_.size());
  int oldObj = distance(solution_[e1], solution_[e1+1]) + distance(solution_[e2], solution_[e2+1]);
  int newObj = distance(solution_[e1], solution_[e2]) + distance(solution_[e1+1], solution_[e2+1]);
  int gain = newObj - oldObj;
  if (gain > 0) {
    //cout << "Improvement found: " << gain << endl;
    reverse(solution_.begin() + min(e1, e2) + 1, solution_.begin() + max(e1, e2) + 1);
  }
  else if (gain == 0 && bernoulli_distribution(0.1)(rgen_)) {
    //cout << "Plateau found " << endl;
    reverse(solution_.begin() + min(e1, e2) + 1, solution_.begin() + max(e1, e2) + 1);
  }
}

void Problem::tryExchange() {
  int subrangeSize = 500;
  int e1 = uniform_int_distribution<int>(0, nbSlides()-1)(rgen_);
  //int e2 = uniform_int_distribution<int>(0, nbSlides()-2)(rgen_);
  int offset = uniform_int_distribution<int>(2, subrangeSize)(rgen_);
  int e2 = (e1 + offset) % nbSlides();
  if (abs(e1 - e2) <= 1) return;
  Roaring tags1 = tags(solution_[e1]);
  Roaring tags2 = tags(solution_[e2]);
  int oldObj = 0;
  int newObj = 0;
  if (e1 > 0) {
    Roaring prev = tags(solution_[e1-1]);
    oldObj += distance(tags1, prev);
    newObj += distance(tags2, prev);
  }
  if (e2 > 0) {
    Roaring prev = tags(solution_[e2-1]);
    newObj += distance(tags1, prev);
    oldObj += distance(tags2, prev);
  }
  if (e1 + 1 < nbSlides()) {
    Roaring prev = tags(solution_[e1+1]);
    oldObj += distance(tags1, prev);
    newObj += distance(tags2, prev);
  }
  if (e2 + 1 < nbSlides()) {
    Roaring prev = tags(solution_[e2+1]);
    newObj += distance(tags1, prev);
    oldObj += distance(tags2, prev);
  }
  int gain = newObj - oldObj;
  if (gain >= 0) {
    //cout << "Improvement found: " << gain << endl;
    swap(solution_[e1], solution_[e2]);
  }
}

void Problem::tryVerticalExchange() {
  int subrangeSize = 500;
  int e1 = uniform_int_distribution<int>(0, nbSlides()-1)(rgen_);
  //int e2 = uniform_int_distribution<int>(0, nbSlides()-2)(rgen_);
  int offset = uniform_int_distribution<int>(2, subrangeSize)(rgen_);
  int e2 = (e1 + offset) % nbSlides();
  if (abs(e1 - e2) <= 1) return;
  if (solution_[e1].second < 0 || solution_[e2].second < 0) return;

  pair<int, int> before1 = solution_[e1];
  pair<int, int> before2 = solution_[e2];
  bool first1 = bernoulli_distribution()(rgen_);
  bool first2 = bernoulli_distribution()(rgen_);
  pair<int, int> after1 = solution_[e1];
  pair<int, int> after2 = solution_[e2];
  if (first1 && first2)   swap(after1.first, after2.first);
  if (first1 && !first2)  swap(after1.first, after2.second);
  if (!first1 && first2)  swap(after1.second, after2.first);
  if (!first1 && !first2) swap(after1.second, after2.second);

  Roaring tb1 = tags(before1);
  Roaring tb2 = tags(before2);
  Roaring ta1 = tags(after1);
  Roaring ta2 = tags(after2);

  int oldObj = 0;
  int newObj = 0;
  if (e1 > 0) {
    Roaring prev = tags(solution_[e1-1]);
    oldObj += distance(tb1, prev);
    newObj += distance(ta1, prev);
  }
  if (e2 > 0) {
    Roaring prev = tags(solution_[e2-1]);
    oldObj += distance(tb2, prev);
    newObj += distance(ta2, prev);
  }
  if (e1 + 1 < nbSlides()) {
    Roaring prev = tags(solution_[e1+1]);
    oldObj += distance(tb1, prev);
    newObj += distance(ta1, prev);
  }
  if (e2 + 1 < nbSlides()) {
    Roaring prev = tags(solution_[e2+1]);
    oldObj += distance(tb2, prev);
    newObj += distance(ta2, prev);
  }
  int gain = newObj - oldObj;
  if (gain > 0) {
    //cout << "Improvement found: " << gain << endl;
    solution_[e1] = after1;
    solution_[e2] = after2;
  }
}

void Problem::reoptimizeLocal() {
  for (int i = 0; i < 1000000; ++i) {
    try2Opt();
    tryVerticalExchange();
    if (bernoulli_distribution(0.1)(rgen_))
        tryExchange();
  }
}

void Problem::reoptimizeAssignment() {
  for (int i = 0; i < 1; ++i) {
    optimizeVerticalAssignment();
  }
}

void Problem::check() const {
  unordered_set<int> used;
  for (pair<int, int> s : solution_) {
    assert (s.first >= 0 && s.first < nbImages());
    assert (s.second < nbImages());
    if (s.second < 0) {
      assert (!imageVertical_[s.first]);
      assert (!used.count(s.first));
      used.insert(s.first);
    }
    else {
      assert (imageVertical_[s.first] && imageVertical_[s.second]);
      assert (!used.count(s.first));
      used.insert(s.first);
      assert (!used.count(s.second));
      used.insert(s.second);
    }
  }
  assert (used.size() == nbImages());
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "Usage: hashcode inputFile [outputFile]" << endl;
    return 1;
  }

  Problem pb = Problem::read(argv[1]);
  //pb.report();
  if (argc >= 3) {
    pb.readSolution(argv[2]);
    pb.check();
  }
  //pb.initNaive();
  //pb.initGreedy();
  cout << "Initial objective value: " << pb.objective() << endl;

  pb.reoptimizeLocal();
  pb.reoptimizeAssignment();
  pb.check();

  cout << "Objective value: " << pb.objective() << endl;
  cout << endl;

  if (argc >= 3) {
    pb.writeSolution(argv[2]);
  }
}
