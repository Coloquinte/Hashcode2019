
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "roaring.hh"
#include "roaring.c"

using namespace std;

class Problem {
 public:
  static Problem read(string filename);
  void writeSolution(string filename) const;

  int nbImages() const { return imageTags_.size(); }
  int nbVerticalImages() const;
  int nbTags() const { return nbTags_; }
  int nbTagReferences() const { return nbTagReferences_; }
  int objective() const;

  void initSolution();
  void report();

 private:
  Roaring tags(pair<int, int> slide) const;
  int distance(pair<int, int> s1, pair<int, int> s2) const;

 private:
  Problem()
    : nbTags_(0)
    , nbTagReferences_(0)
    {}
  int nbTags_;
  int nbTagReferences_;
  vector<Roaring> imageTags_;
  vector<bool> imageVertical_;
  vector<pair<int, int> > solution_;
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

void Problem::writeSolution(string filename) const {
  ofstream s(filename);
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

void Problem::report() {
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
  Roaring tags1 = tags(s1);
  Roaring tags2 = tags(s2);
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

void Problem::initSolution() {
  vector<pair<int, int> > slides;
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
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "Usage: hashcode inputFile [outputFile]" << endl;
    return 1;
  }

  Problem pb = Problem::read(argv[1]);
  pb.report();
  pb.initSolution();

  cout << "Objective value: " << pb.objective() << endl;

  if (argc >= 3) {
    pb.writeSolution(argv[2]);
  }
}
