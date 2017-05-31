/*    phasogram: generate phasogram from mapped reads
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
#include <numeric>
#include <algorithm>

#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::min;

struct ReadPos {
  string chrom;
  size_t start;
  size_t end;
  char strand;
};

std::istream& 
operator>>(std::istream& in, ReadPos &rp) {
  string line;
  if (!getline(in, line))
    return in; 

  string name, seq, scr;
  double score;
  std::istringstream iss(line);
  
  if (!(iss >> rp.chrom >> rp.start >> rp.end >> name
        >> score >> rp.strand >> seq >> scr))
    throw SMITHLABException("bad mapped read file");
  
  if (rp.strand == '-') {
    const size_t start = rp.end;
    rp.end = rp.start;
    rp.start = start;
  }
  return in;
}




static size_t
dist_to_idx(const size_t phase, const size_t start, const size_t bin) {
  return (phase - start) / bin;
}



static bool
region_precedes_site(const vector<GenomicRegion> &regions,
                     const size_t region_idx, const ReadPos &rp) {
  return (region_idx < regions.size() &&
          (regions[region_idx].get_chrom() < rp.chrom ||
           (regions[region_idx].get_chrom() == rp.chrom &&
            regions[region_idx].get_end() <= rp.start)));
}


static bool
region_contains_site(const vector<GenomicRegion> &regions,
                     const size_t region_idx, const ReadPos &rp) {
  // !!! ASSUMES REGION DOES NOT PRECEDE SITE
  return (region_idx < regions.size() &&
          regions[region_idx].get_chrom() == rp.chrom &&
          regions[region_idx].get_start() <= rp.start &&
          regions[region_idx].get_end() > rp.end);
}

/*
static size_t
boundary_next_to_site(const vector<GenomicRegion> &regions,
                      size_t &region_idx, const ReadPos &rp) {
  size_t b = rp.start;
  while (region_precedes_site(regions, region_idx, rp))
    ++region_idx;
  const bool contained = region_contains_site(regions, region_idx, rp);
  if (contained)
    b = regions[region_idx].get_end()-1;
  else
    if (region_idx < regions.size() &&
        regions[region_idx].get_chrom() == rp.chrom)
      b = regions[region_idx].get_start()-1;
    else b = std::numeric_limits<size_t>::max(); 

  return b;
}
*/

static void
process_chrom(const vector<ReadPos> &reads,
              vector<size_t> &phg, const size_t max_dist,
              const size_t start, const size_t bin) {
  
  const size_t L = phg.size(); 
  for (size_t i = 0; i < reads.size(); ++i) {
    size_t pos_limit = reads[i].start + max_dist;
    size_t j = i + 1;
    while (j < reads.size() &&
           min(reads[j].start, reads[j].end) <= pos_limit) {
        if (reads[j].strand == reads[i].strand) {
          if (reads[j].start >= reads[i].start) {
            const size_t phase = reads[j].start - reads[i].start;
            const size_t idx = dist_to_idx(phase, start, bin);
            if (idx < L)
              phg[idx] += 1;
          } else {
            cerr << "Wrongly sorted read pairs:" << endl;
            cerr << reads[i].chrom << "\t" << reads[i].start << "\t"
                 << reads[i].end << "\t" << reads[i].strand << endl;
            cerr << reads[j].chrom << "\t" << reads[j].start << "\t"
                 << reads[j].end << "\t" << reads[j].strand << endl;
          }
        }
      ++j;
    }
  }
}



static void
process_chrom_adj(const vector<ReadPos> &reads,
                  vector<size_t> &phg, const size_t max_dist,
                  const size_t start, const size_t bin) {
 
  const size_t L = phg.size(); 
  for (size_t i = 0; i < reads.size(); ++i) {
    size_t pos_limit = reads[i].start + max_dist;
    if (i < (reads.size() - 1) && reads[i+1].start <= pos_limit &&
        reads[i+1].strand == reads[i].strand) {
      if (reads[i+1].start >= reads[i].start) {
        const size_t phase = reads[i+1].start - reads[i].start;
        const size_t idx = dist_to_idx(phase, start, bin);
        if (idx < L)
          phg[idx] += 1;
      } else {
        cerr << "Wrongly sorted read pairs:" << endl;
        cerr << reads[i].chrom << "\t" << reads[i].start << "\t"
             << reads[i].end << "\t" << reads[i].strand << endl;
        cerr << reads[i+1].chrom << "\t" << reads[i+1].start << "\t"
             << reads[i+1].end << "\t" << reads[i+1].strand << endl;
      }
    }
  }
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
read_allowed(const bool exclude_regions,
             const vector<GenomicRegion> &regions,
             const ReadPos &rp, size_t &region_idx) {
  while (region_precedes_site(regions, region_idx, rp))
    ++region_idx;
  const bool contained = region_contains_site(regions, region_idx, rp);
  return (exclude_regions ? !contained : contained);
}


int main(int argc, const char **argv) {

  try {

    string outfile;
    size_t start = 1;
    size_t end = 2000;
    size_t bin = 10;
    size_t limit = 2000;

    string regions_file;
    bool adjacent = false;
    bool exclude_regions = false;

    bool VERBOSE = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<mr-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("start", 's', "Phasogram start (default: 1)",
                      false , start);
    opt_parse.add_opt("end", 'e', "Phasogram end (default: 2000)",
                      false , end);
    opt_parse.add_opt("bin", 'b', "Phasogram bin size (default: 10)",
                      false , bin);
    opt_parse.add_opt("limit", 'l',
                      "Maximum distance of pairs of reads (default: 2000)",
                      false , limit);
    opt_parse.add_opt("adjacent", 'A',
                      "only consider adjacent reads (default: F)",
                      false , adjacent);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , VERBOSE);
    opt_parse.add_opt("regions", '\0', "regions file", false, regions_file);
    opt_parse.add_opt("exclude", 'E', "exclude regions (default include)",
                      false, exclude_regions);

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
    const string mr_file_name = leftover_args.front();
    /**********************************************************************/

    const size_t filesize = get_filesize(mr_file_name);

    if (VERBOSE)
      cerr << "phasogram_start : " << start << endl
           << "phasogram_end : " << end << endl
           << "phasogram_bin : " << bin << endl
           << "max_dist : " << limit << endl
           << "input_file : " << mr_file_name << endl
           << "input_size : " << filesize << endl;

    std::ifstream in(mr_file_name.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + mr_file_name);

    size_t L = (end - start) / bin + 1;
    vector<size_t> phg(L, 0);
    
    ReadPos rp;
    vector<ReadPos> reads; 
    string prev_chrom;

    vector<GenomicRegion> regions;
    if (!regions_file.empty()) {
      load_regions(regions_file, regions);
      if (VERBOSE)
        cerr << "load regions : " << regions.size() << endl;
    }
    size_t region_idx = 0;

    while (in >> rp) {
      if (rp.chrom != prev_chrom) {
        if (adjacent)
          process_chrom_adj(reads, phg, limit, start, bin);
        else
          process_chrom(reads, phg, limit, start, bin);
        reads.clear();
      }

      if (regions_file.empty() ||
          read_allowed(exclude_regions, regions, rp, region_idx)) {
        reads.push_back(rp);
      }
      prev_chrom.swap(rp.chrom);
      
    }

    if (adjacent)
      process_chrom_adj(reads, phg, limit, start, bin);
    else
      process_chrom(reads, phg, limit, start, bin);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    size_t pos = start;
    for (size_t i = 0; i < L; ++i) {
      out << pos << '\t' << phg[i] << endl;
      pos = pos + bin;
    }
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
