/*    methdiff: calculate methylation difference in nucleosome vs.
 *              flanks
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
#include <stdlib.h>
#include <gsl/gsl_histogram.h>


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
using std::max;
using std::sqrt;



struct NMcouple {
  NMcouple() : chrom("chr0"), start(0), end(0), tnu(0), tme(0),
               n_nu(1), n_me(0), max_me(0), min_me(0) {}
  string chrom;
  size_t start;
  size_t end;
  double tnu;
  double tme;
  size_t n_nu;
  size_t n_me;
  double max_me;
  double min_me;
};


static NMcouple
get_linker(const bool left, const NMcouple &t, const size_t flank) {
  const size_t lowest = 0;
  NMcouple linker;
  linker.chrom = t.chrom;
  linker.start = left ? max(t.start - flank, lowest) : t.end;
  linker.end = left ? t.start : t.end + flank;
  return linker;
}


static void
trim_linker(const bool left, NMcouple &linker, const NMcouple &nb) {
  if (left) {
    linker.start = max(linker.start, nb.end);
  } else {
    linker.end = min(linker.end, nb.start);
  }
}


std::istream& 
operator>>(std::istream& in, NMcouple &ns) {
  string line;
  if (!getline(in, line))
    return in; 
  
  string junk;
  std::istringstream iss(line);
  if (!(iss >> ns.chrom >> ns.start >> ns.end >> junk >> ns.tnu >> junk))
    throw SMITHLABException("bad mapped read file");
  ns.n_nu = 1;
  return in;
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
    string histfile;
    size_t flank = 40;
    size_t min_reads = 1;
    size_t min_sites = 1;
    double hbin = 0.01;

    string nucleo_file;
    string regions_file;

    bool VERBOSE = false;

    bool intact_linker = false;
    bool allow_nuc_overlap = false;
    bool exclude_regions = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<meth-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("hist", 'h', "Name of histogram file",
                      false , histfile);
    opt_parse.add_opt("nucle", '\0', "nucleosome peak file", true,
                      nucleo_file);
    opt_parse.add_opt("flank", 'f', "linker length of nucleosome (default: 40)",
                      false , flank);
    opt_parse.add_opt("reads", 'r', "min reads for counting a CpG site (default: 1)",
                      false , min_reads);
    opt_parse.add_opt("min_sites", 's', "min CpG sites (default: 1)",
                      false , min_sites);
    opt_parse.add_opt("hbin", 'H', "histogram bin size (default: 0.01)",
                      false , hbin);
    opt_parse.add_opt("regions", '\0', "regions file", false,
                      regions_file);
    opt_parse.add_opt("intact", 'I',
                      "keep intact flanking linker no trim (defulat: false)",
                      false, intact_linker);
    opt_parse.add_opt("overlap", 'O',
                      "allow and merge overlapping peaks (default: false)",
                      false, allow_nuc_overlap);
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
      cerr << "min_reads: " << min_reads << endl
           << "min_sites : " << min_sites << endl
           << "flank : " << flank << endl
           << "input_file : " << cpg_file_name << endl
           << "input_size : " << filesize << endl
           << "keep overlap nucleosomes: " << allow_nuc_overlap << endl
           << "keep intact linkers: " << intact_linker << endl << endl;

    std::ifstream in(cpg_file_name.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + cpg_file_name);
    
    std::ifstream inns(nucleo_file.c_str());
    if (!inns)
      throw SMITHLABException("bad nucleosome file: " + nucleo_file);


    vector<GenomicRegion> regions;
    if (!regions_file.empty()) {
      load_regions(regions_file, regions);
      if (VERBOSE)
        cerr << "load regions : " << regions.size() << endl << endl;
    }

    NMcouple nu;
    vector<NMcouple> nss;
    
    while (!inns.eof() && inns >> nu) {
      if (!nss.empty() && nu.start <= nss.back().end) {
        // detected overlapping nucleosome
        nss.back().end = nu.end;
        nss.back().tnu += nu.tnu;
        nss.back().n_nu++;
      } else {
        nss.push_back(nu);
      }
    }
    if (VERBOSE)
      cerr << "load nucleosomes : " << nss.size() << endl << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    
    vector<double> mdiff;  
    const size_t n = nss.size();
    MSite s;
    s.chrom = "chr0";
    s.pos = 0;
    size_t region_idx = 0;
    for (size_t i = 0; i < n; ++i) {
      if (!allow_nuc_overlap && nss[i].n_nu > 1)
        continue;
      if (!regions_file.empty() &&
          !window_allowed(exclude_regions, regions, nss[i], region_idx))
        continue;

      // locate linker
      NMcouple ll = get_linker(true, nss[i], flank);
      NMcouple rl = get_linker(false, nss[i], flank);
      if (!intact_linker && i > 0)
        trim_linker(true, ll, nss[i-1]);
      if (!intact_linker && i < n-1)
        trim_linker(false, rl, nss[i+1]);
      // assign CpGs
      while (site_precedes_region(s, ll) && !in.eof()) {
        in >> s;
      }
      while(site_in_region(s, ll) && !in.eof()) {
        // in left linker
        if (s.n_reads > min_reads) {
          ll.tme += s.meth;
          ll.n_me ++;
          ll.min_me = ll.n_me > 0 ? min(ll.min_me, s.meth) : 0;
        }
        in >> s; 
      }
      while(site_in_region(s, nss[i]) && !in.eof()) {
        // in nucleosome
        if (s.n_reads > min_reads) {
          nss[i].tme += s.meth;
          nss[i].n_me ++;
          nss[i].max_me = max(nss[i].max_me, s.meth);
        }
        in >> s; 
      }
      while(site_in_region(s, rl) && !in.eof()) {
        // in right linker
        if (s.n_reads > min_reads) {
          rl.tme += s.meth;
          rl.n_me ++;
          rl.min_me = rl.n_me > 0 ? min(rl.min_me, s.meth) : 0;
        }
        in >> s; 
      }

      if (ll.n_me > min_sites && nss[i].n_me > min_sites &&
          rl.n_me > min_sites) {
        const double lm = ll.tme / ll.n_me;
        const double nm = nss[i].tme / nss[i].n_me;
        const double rm = rl.tme / rl.n_me;
        out << nss[i].chrom << '\t' << nss[i].start << '\t' << nss[i].end
            << '\t' << nss[i].tnu / nss[i].n_nu << '\t'
            << lm << '\t' << nm << '\t' << rm << '\t'
            << ll.min_me << '\t' << nss[i].max_me << '\t' << rl.min_me << endl;
        mdiff.push_back(nm - (lm + rm) / 2);
      }
    }
    
    FILE * phist = fopen(histfile.c_str(), "w");
    gsl_histogram * h = gsl_histogram_alloc(2/hbin);
    gsl_histogram_set_ranges_uniform (h, -1, 1);

    if (VERBOSE)
      cerr << "compute histogram: " << hbin << endl << endl;
    
    for(size_t i = 0; i < mdiff.size(); ++i) {
      gsl_histogram_increment(h, mdiff[i]);
    }
    gsl_histogram_fprintf (phist, h, "%g", "%g");
    //gsl_histogram_fwrite(phist, h);
    gsl_histogram_free(h);

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
