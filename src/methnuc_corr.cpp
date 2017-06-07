/*    methnuc_corr: calculate correlation between methylation and nucleosome
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith
 *
 *    Authors: Liz Ji, Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License as
 *    published by the Free Software Foundation, either version 3 of
 *    the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with this program.  If not, see
 *    <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <tr1/cmath>
#include <numeric>
#include <algorithm>

#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::sort;
using std::min;
using std::sqrt;



struct NMcouple {
  NMcouple() : chrom("chr0"), start(0), end(0), nu(0), tme(0), n_nu(0),
               n_me(0), score(0) {}
  string chrom;
  size_t start;
  size_t end;
  double nu;
  double tme;
  size_t n_nu;
  size_t n_me;
  double score;

};


std::istream& 
operator>>(std::istream& in, NMcouple &ns) {
  string line;
  if (!getline(in, line))
    return in; 

  std::istringstream iss(line);
  if (!(iss >> ns.chrom >> ns.start >> ns.end >> ns.nu))
    throw SMITHLABException("bad mapped read file");
  return in;
}

static size_t
loop_rshift(const size_t idx, const size_t l) {
  return(min(idx + 1, l-1));
}


template <class T> static double
corr(const vector<T> &x, const vector<T> &y) {

  const T n = x.size();

  const T x_sum = accumulate(x.begin(), x.end(), 0.0);
  const T y_sum = accumulate(y.begin(), y.end(), 0.0);

  const T ip = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

  const T x_ss = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const T y_ss = std::inner_product(y.begin(), y.end(), y.begin(), 0.0);

  const T p1 = ip - x_sum*y_sum / n;
  const T p2 = std::sqrt(x_ss - x_sum*x_sum / n);
  const T p3 = std::sqrt(y_ss - y_sum*y_sum / n);
  const T cor = p1 / (p2*p3);
  //if (cor > 1 || cor < -1) {
  //  cerr << p1 << '\t' << p2 << '\t' << p3 << '\t' << cor << endl;
  //}

  return cor;
}


static void
load_regions(const string &regions_file,
             vector<GenomicRegion> &regions) {

  std::ifstream in(regions_file.c_str());
  if (!in)
    throw SMITHLABException("bad regions file: " + regions_file);

  GenomicRegion r;
  while (in >> r)
    regions.push_back(r);

  if (!check_sorted(regions))
    throw SMITHLABException("regions file not sorted: " + regions_file);
}




static bool
site_precedes_region(const MSite &s, const NMcouple &t) {
  return (s.chrom < t.chrom ||
          (s.chrom == t.chrom && s.pos < t.start));
}


static bool
site_in_region(const MSite &s, const NMcouple &t) {
  return (s.chrom == t.chrom && s.pos >= t.start &&
           s.pos <= t.end);
}


static bool
region_precedes_target(const vector<GenomicRegion> &regions,
                       const size_t region_idx,
                       const NMcouple &t) {
  return (region_idx < regions.size() &&
          (regions[region_idx].get_chrom() < t.chrom ||
           (regions[region_idx].get_chrom() == t.chrom &&
            regions[region_idx].get_end() < t.start)));
}


static bool
region_overlaps_target(const vector<GenomicRegion> &regions,
                       const size_t region_idx,
                       const NMcouple &t) {
  // !!! ASSUMES REGION DOES NOT PRECEDE SITE
  return (region_idx < regions.size() &&
          regions[region_idx].get_chrom() == t.chrom &&
          regions[region_idx].get_end() >= t.start &&
          regions[region_idx].get_start() <= t.end);
}


static void
process_window(const vector<NMcouple> &nss, NMcouple &w,
              vector<NMcouple> &wins, const size_t min_reads,
              const size_t min_sites) {
  const size_t n = nss.size();
  vector<double> x;
  vector<double> y;
  // go through the nucleosome bins and filter out
  for (size_t i = 0; i < n; i++) {
    if (nss[i].n_me >= min_reads) {
      x.push_back(nss[i].nu);
      y.push_back(nss[i].tme / nss[i].n_me);
    }
  }
  if (x.size() >= min_sites) {
    w.score = corr(x, y);
    wins.push_back(w);
  }
}


static bool
window_allowed(const bool exclude_regions,
              const vector<GenomicRegion> &regions,
              const NMcouple &window, size_t &region_idx) {
  while (region_precedes_target(regions, region_idx, window))
    ++region_idx;
  const bool overlap = region_overlaps_target(regions, region_idx, window);
  return (exclude_regions ? !overlap : overlap);
}



int main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    size_t window = 2000;
    size_t bin = 100;
    size_t min_reads = 1;
    size_t min_sites = 10;

    string nucleo_file;
    string regions_file;

    bool VERBOSE = false;

    bool sliding = false;
    bool exclude_regions = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<meth-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("nucle", '\0', "nucleosome occupancy signal file", true,
                      nucleo_file);
    opt_parse.add_opt("window", 'w',
                      "window for computing correlation (default: 1000)",
                      false , window);
    opt_parse.add_opt("bin", 'b', "bin size (default: 10)",
                      false , bin);
    opt_parse.add_opt("reads", 'r', "min reads for nucbin (default: 1)",
                      false , min_reads);
    opt_parse.add_opt("g", 's', "min sites for correlation (default: 10)",
                      false , min_sites);
    opt_parse.add_opt("regions", '\0', "regions file", false,
                      regions_file);
    opt_parse.add_opt("slide", 'S', "sliding window (default: false)",
                      false, sliding);
    opt_parse.add_opt("exclude", 'E', "exclude regions (default: include)",
                      false, exclude_regions);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpg_file_name = leftover_args.front();
    /**********************************************************************/

    const size_t filesize = get_filesize(cpg_file_name);

    if (VERBOSE)
      cerr << "min_sites : " << min_sites << endl
           << "window : " << window << endl
           << "bin : " << bin << endl
           << "input_file : " << cpg_file_name << endl
           << "input_size : " << filesize << endl;

    std::ifstream in(cpg_file_name.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + cpg_file_name);
    
    std::ifstream inns(nucleo_file.c_str());
    if (!inns)
      throw SMITHLABException("bad nucleosome file: " + nucleo_file);

    MSite s;
    s.chrom = "chr0";
    s.pos = 0;
    vector<NMcouple> nss;
    NMcouple nu;

    vector<GenomicRegion> regions;
    if (!regions_file.empty()) {
      load_regions(regions_file, regions);
      if (VERBOSE)
        cerr << "load regions : " << regions.size() << endl;
    }

    vector<NMcouple> wins;
    NMcouple w;
    size_t left_most = 0;
    size_t right_most = 0;
    size_t region_idx = 0;
    while (!inns.eof() && inns >> nu) {
      // load CpGs
      while (site_precedes_region(s, nu) && !in.eof()) {
        in >> s;
      }
      while(site_in_region(s, nu) && !in.eof()) {
        nu.tme += s.meth;
        nu.n_me ++;
        in >> s; 
      }

      if (nu.chrom != w.chrom || nu.end > w.end) {
        // open/expand new window
        if (!nss.empty() &&
            (regions_file.empty() ||
             window_allowed(exclude_regions, regions, w, region_idx))) {
          w.end = nss[right_most].end; 
          process_window(nss, w, wins, min_reads, min_sites);
        }
        
        if (nu.chrom != w.chrom) {
          // enter new chromosome
          nss.clear();
          nss.push_back(nu);
          left_most = 0;
          right_most = 0;
          w.chrom = nu.chrom;
          w.start = nu.start;
          w.end = nu.start + window;
        } else {
          // enter new window
          if (sliding) {
            // sliding window: replace the first nu with the new one
            nss[left_most] = nu;
            right_most = left_most;
            left_most = loop_rshift(left_most, nss.size()); 
            w.start = nss[left_most].start; 
            w.end = nu.end;
          } else {
            nss.clear();
            nss.push_back(nu);
            left_most = 0;
            right_most = 0;
            w.start = nu.start;
            w.end = nu.start + window;
          }
        }
      } else {
        nss.push_back(nu);
        right_most++;
      }
    }
    
    if (!nss.empty() &&
        (regions_file.empty() ||
          window_allowed(exclude_regions, regions, w, region_idx))) {
      w.end = nss[right_most].end; 
      process_window(nss, w, wins, min_reads, min_sites);
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
  
    for (size_t i = 1; i < wins.size(); ++i)
      out << wins[i].chrom << '\t' << wins[i].start << '\t'
          << wins[i].end << '\t' << wins[i].score << endl;
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
