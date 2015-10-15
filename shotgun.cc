#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

// if at any point the edit distance is known to be D or greater, where D > max_edits,
// then this function returns -D
int edit_distance (const string& x, const string& y, vector<vector<int> >& ed_matrix, int max_edits)
{
  const int xsize = x.size();
  const int ysize = y.size();
  const int min_edits_by_len = abs (xsize - ysize);
  if (min_edits_by_len > max_edits)
    return -min_edits_by_len;
  ed_matrix[0][0] = 0;
  for (int i = 1; i <= xsize; ++i) {
    int row_min_sc = max_edits + 1;
    for (int j = 1; j <= ysize; ++j)
      {
	const int match_sc = ed_matrix[i-1][j-1] + (x[i-1] == y[j-1] ? 0 : 1);
	const int ins_sc = ed_matrix[i-1][j] + 1;
	const int del_sc = ed_matrix[i][j-1] + 1;

	const int gap_sc = min (ins_sc, del_sc);
	const int cell_sc = min (match_sc, gap_sc);
	ed_matrix[i][j] = cell_sc;
	row_min_sc = min (row_min_sc, cell_sc);
      }
    if (row_min_sc > max_edits)
      return -row_min_sc;
  }

  const int dist = ed_matrix[xsize][ysize];

  // distance should never be zero...
  if (dist == 0 && x != y)
    {
      cerr << "Error: calculated zero edit distance for non-identical words " << x << " and " << y << "\n";
      cerr << '*';
      for (int i = 1; i <= xsize; ++i)
	cerr << '\t' << x[i-1];
      cerr << '\n';
      for (int j = 0; j <= ysize; ++j)
	{
	  cerr << (j==0 ? '*' : y[j-1]);
	  for (int i = 0; i <= xsize; ++i)
	    cerr << '\t' << ed_matrix[i][j];
	  cerr << '\n';
	}
      exit(1);
    }

  return dist;
}

typedef vector<map<int,int> > EDCache;
EDCache edit_distance_cache;
int get_edit_distance (const vector<string>& words, int x, int y, vector<vector<int> >& ed_matrix, int max_edits) {
  map<int,int>::const_iterator ed_iter = edit_distance_cache[x].find(y);
  if (ed_iter != edit_distance_cache[x].end()) {
    const int cached_ed = ed_iter->second;
    if (cached_ed >= 0 || (-cached_ed > max_edits))
      return cached_ed;
  }

  const int ed = edit_distance (words[x], words[y], ed_matrix, max_edits);
  edit_distance_cache[x][y] = ed;

  return ed;
}

void get_neighbors (const vector<string>& words, int src, int max_edits, vector<vector<int> >& ed_matrix, vector<int>& nbr, vector<int>& dist)
{
  nbr.clear();
  dist.clear();
  for (int dest = 0; dest < (int) words.size(); ++dest)
    if (dest != src)
      {
	const int ed = get_edit_distance (words, src, dest, ed_matrix, max_edits);
	if (ed >= 0 && ed <= max_edits)
	  {
	    nbr.push_back (dest);
	    dist.push_back (ed);
	  }
      }
}

void show_string (ostream& out, const vector<int>& sentence, const vector<string>& words)
{
  for (int n = 0; n < (int) sentence.size(); ++n)
    out << words[sentence[n]] << ' ';
  out << '\n';
}

void show_string (const vector<int>& sentence, const vector<string>& words)
{
  show_string (cout, sentence, words);
}

vector<int> evolve_sentence (const vector<int>& ancestor, int edits, double max_edits_per_letter, const vector<string>& words, vector<vector<int> >& ed_matrix)
{
  const int len = ancestor.size();
  vector<int> remaining_edits (len);
  vector<int> descendant (ancestor);
  for (int n = 0; n < len; ++n)
    remaining_edits[n] = max (1, (int) (.5 + max_edits_per_letter * (double) words[descendant[n]].size()));
  bool first_edit = true;
  while (edits > 0)
    {
      // log
      show_string (cerr, descendant, words);
      cerr << " (" << edits << " edits left)\n";
      // evolve
      vector<int> pos, new_word, distance, nbr, nbr_dist;
      vector<double> weight;
      double total_weight = 0;
      for (int n = 0; n < len; ++n)
	{
	  get_neighbors (words, descendant[n], min (edits, remaining_edits[n]), ed_matrix, nbr, nbr_dist);
	  if (first_edit) {
	    //	    cerr << nbr.size() << " neighbors of '" << words[descendant[n]] << "'\n";
	    if (nbr.size() == 0)
	      cerr << "Warning: stuck on '" << words[descendant[n]] << "', no neighbors within " << remaining_edits[n] << " edits\n";
	  }
	  const double weight_n = 1. / (double) nbr.size();
	  for (int k = 0; k < (int) nbr.size(); ++k)
	    {
	      pos.push_back (n);
	      new_word.push_back (nbr[k]);
	      distance.push_back (nbr_dist[k]);
	      weight.push_back (weight_n);
	      total_weight += weight_n;
	    }
	}

      if (pos.size() == 0)
	{
	  cout << "Ran out of possibilities -- string:\n";
	  show_string (descendant, words);
	  exit(1);
	}

      double mutation_weight = drand48() * total_weight;
      int mutation;
      for (mutation = 0; mutation_weight > 0; ++mutation)
	mutation_weight -= weight[mutation];

      descendant[pos[mutation]] = new_word[mutation];
      remaining_edits[pos[mutation]] -= distance[mutation];
      edits -= distance[mutation];

      first_edit = false;
    }

  return descendant;
}

int main (int argc, char** argv)
{
  if (argc != 6)
    {
      cerr << "Usage: " << argv[0] << " [dictionary file] [mean edits per word, per read] [max edits per letter, per read] [average coverage] [read length, in words] [source file]\n";
      return 1;
    }

  vector<string> word;
  map<string,int> word2int;
  string s;
  ifstream infile (argv[1]);
  int max_len = 0;
  while (!infile.eof())
    {
      infile >> s;
      if (!infile.eof())
	{
	  std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);  // transform to upper case
	  word2int[s] = word.size();
	  word.push_back (s);
	  if (s.size() > max_len)
	    max_len = s.size();
	}
    }
  edit_distance_cache = EDCache (word.size());


  // TODO: update all this
  
  
  vector<int> seed;
  for (int n = 0; n+5 < argc; ++n)
    {
      string s (argv[n+5]);
      std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);  // transform to upper case
      if (word2int.find(s) == word2int.end())
	{
	  cerr << "Can't find word " << s << " in dictionary file\n";
	  exit(1);
	}
      seed.push_back (word2int[s]);
      //      cerr << " [Word " << n+1 << ": " << s << " = " << seed.back() << "]\n";
    }

  const double edits_per_word = atof (argv[2]);
  const int edits_per_branch = (int) (.5 + edits_per_word * (double) seed.size());
  const double max_edits_per_letter = atof (argv[3]);
  const double coverage = atof (argv[4]);
  const int read_length = atoi (argv[5]);

  vector<vector<int> > ed_matrix (max_len + 1, vector<int> (max_len + 1, 0));
  ed_matrix[0][0] = 0;
  for (int i = 1; i <= max_len; ++i)
    ed_matrix[i][0] = ed_matrix[0][i] = i;

  evolve_subtree (seed, tree_depth, string(), edits_per_branch, max_edits_per_letter, word, ed_matrix);

  return 0;
}
