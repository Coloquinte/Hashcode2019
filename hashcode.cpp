
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

class Problem {
 public:
  static Problem read(string filename);
  void writeSolution(string filename) const;
  void report();

 private:
  Problem()
    : nbTags_(0)
    , nbTagReferences_(0)
    {}
  int nbTags_;
  int nbTagReferences_;
  vector<unordered_set<int> > imageTags_;
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
      pb.imageTags_[i].insert(tagId);
    }
    pb.nbTagReferences_ += pb.imageTags_[i].size();
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

void Problem::report() {
  cout << imageTags_.size() << " images" << endl; 
  cout << nbTags_ << " tags" << endl; 
  cout << nbTagReferences_ << " tag references" << endl; 
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "Usage: hashcode inputFile [outputFile]" << endl;
    return 1;
  }

  Problem pb = Problem::read(argv[1]);
  pb.report();

  if (argc >= 3) {
    pb.writeSolution(argv[2]);
  }
}
