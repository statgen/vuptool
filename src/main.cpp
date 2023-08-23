/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <vector>
#include <list>
#include <cstdint>
#include <string>

#include <savvy/reader.hpp>
#include <shrinkwrap/istream.hpp>

#include "inv_norm.hpp"
#include "getopt_wrapper.hpp"

typedef float signal_t;
enum class method_t { filter, mask, regress };

std::vector<std::string> split_string_to_vector(const std::string& in, char delim)
{
  std::vector<std::string> ret;
  std::string::const_iterator d;
  std::string token;
  auto s = in.begin();
  auto e = in.end();
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(s, d);
    s = d != e ? d + 1 : d;
  }
  ret.emplace_back(s,d);
  return ret;
}

class vuptool_args : public getopt_wrapper
{
private:
  std::vector<option> long_options_;
  std::string vcf_path_;
  std::string methy_path_;
  std::string manifest_path_;
  std::string output_path_ = "/dev/stdout";
  double filter_threhold_ = 0.05;
  method_t method_ = method_t::regress;
  bool version_ = false;
  bool help_ = false;
public:
  vuptool_args() :
    getopt_wrapper("Usage: vuptool [opts ...] methy_table.bed.gz genotypes.{sav,bcf,vcf.gz} epic_array_manifest.tsv.gz",
      {
        {"method", "<string>", 'm', "filter, mask, or regress (default: regress)"},
        {"help", "", 'h', "Print usage"},
        {"version", "", 'v', "Print version"},
        {"output", "<file>", 'o', "Output path (default: /dev/stdout)"},
        {"filter-threshold", "<float>", 'p', "Probes will be excluded when proportion of samples with variation under the probe >= this value (--method=filter only; default: 0.05)"}
      })
  {
  }

  const std::string& vcf_path() const { return vcf_path_; }
  const std::string& methy_path() const { return methy_path_; }
  const std::string& manifest_path() const { return manifest_path_; }
  const std::string& output_path() const { return output_path_; }


  double filter_threshold() const { return filter_threhold_; }
  method_t method() const { return method_; }
  bool version_is_set() const { return version_; }
  bool help_is_set() const { return help_; }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opt_string_.c_str(), long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'm':
        {
          std::string method_str(optarg ? optarg : "");
          if (method_str == "regress" || method_str == "r")
            method_ = method_t::regress;
          else if (method_str == "mask" || method_str == "m")
            method_ = method_t::mask;
          else if (method_str == "filter" || method_str == "f")
            method_ = method_t::filter;
          else
            return std::cerr << "Error: invalid --method\n", false;
        }
        break;
      case 'h':
        help_ = true;
        return true;
      case 'v':
        version_ = true;
        return true;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'p':
      {
        filter_threhold_ = std::atof(optarg ? optarg : "");
        break;
      }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 3)
    {
      methy_path_ = argv[optind];
      vcf_path_ = argv[optind + 1];
      manifest_path_ = argv[optind + 2];
    }
    else if (remaining_arg_count < 3)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

class methy_t
{
private:
  std::string chrom_;
  std::int64_t pos_;
  std::int64_t pos_end_;
  std::string id_;
  std::vector<signal_t> signals_;
  std::vector<bool> mask_;
  std::list<savvy::compressed_vector<std::int8_t>> predictor_list_;
  bool is_reverse_complement_;
public:
  const std::string& chrom() const { return chrom_; }
  const std::vector<signal_t>& signals() const { return signals_; }

  const std::int64_t start_pos(std::int64_t probe_length) const
  {
    return is_reverse_complement_ ? pos_ : pos_ - probe_length + 2;
  }

  const std::int64_t end_pos(std::int64_t probe_length) const
  { 
    return is_reverse_complement_ ? pos_ + probe_length : pos_ + 2;
  }  

  const std::int64_t distance(const savvy::site_info& var, std::int64_t probe_length) const
  {
    std::int64_t s = start_pos(probe_length);
    std::int64_t e = end_pos(probe_length);
    if (var.pos() < s)
      return var.pos() - s;
    else if (var.pos() > e)
      return var.pos() - e;
    return 0;
  }

  void mask_signals()
  {
    for (std::size_t i = 0; i < mask_.size(); ++i)
    {
      if (mask_[i])
        signals_[i] = -std::numeric_limits<signal_t>::infinity();
    }
  }

  void update_mask(const savvy::compressed_vector<std::int8_t>& gts)
  {
    std::size_t stride = gts.size() / signals_.size();
    for (auto a = gts.begin(); a != gts.end(); ++a)
    {
      mask_[a.offset() / stride] = true;
    }
  }

  double mask_proportion() const
  {
    return double(std::count(mask_.begin(), mask_.end(), true)) / mask_.size();
  }

  void add_predictor(const savvy::compressed_vector<std::int8_t>& gts)
  {
    predictor_list_.push_back(gts);
  }

  void regress_out_genotypes()
  {
    if (predictor_list_.size())
    {

    }
    throw std::runtime_error("TODO: " + std::to_string(__LINE__));
  }
  
  void inv_norm()
  {
    inverse_normalize(signals_);
  }

  static methy_t deserialize(const std::string& line, const std::unordered_set<std::string>& rev_comp_ids)
  {
    methy_t ret;
    auto vec = split_string_to_vector(line, '\t');
    if (vec.size() >= 5)
    {
      ret.chrom_ = vec[0];
      ret.pos_ = std::atoll(vec[1].c_str()) + 1;
      ret.pos_end_ = std::atoll(vec[2].c_str());
      ret.id_ = vec[3];
      ret.is_reverse_complement_ = rev_comp_ids.find(ret.id_) != rev_comp_ids.end();
      ret.signals_.reserve(vec.size() - 4);
      for (auto it = vec.begin() + 4; it != vec.end(); ++it)
        ret.signals_.push_back(std::atof(it->c_str()));
      ret.mask_.resize(ret.signals_.size(), false);
    }
    return ret;
  }

  static bool serialize(const methy_t& m, std::ostream& ofs, const std::string& mask_str)
  {
    ofs << m.chrom_ << "\t" << (m.pos_ - 1) << "\t" << m.pos_end_ << "\t" << m.id_;
    for (auto it = m.signals_.begin(); it != m.signals_.end(); ++it)
    {
      ofs.put('\t');
      if (std::isfinite(*it))
        ofs << *it;
      else
        ofs << mask_str;
    }
    return ofs.put('\n').good();
  }
};

int main(int argc, char** argv)
{

  vuptool_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "vuptool v" << VUPTOOL_VERSION << std::endl;
    return EXIT_SUCCESS;
  }


  std::string mask_string = "NA";
  bool inv_norm = false;
  std::int64_t probe_length = 50;
  method_t method = method_t::mask;
  double filter_threshold = 0.05;

  shrinkwrap::istream methy_file(argv[1]);
  savvy::reader vcf(argv[2]);
  std::ifstream probe_direction_map(argv[3]);
  std::ofstream output_file(argc > 4 ? argv[4] : "/dev/stdout", std::ios::binary);

  std::string line;
  std::unordered_set<std::string> reverse_complement_ids;
  while (std::getline(probe_direction_map, line))
  {
    std::size_t first_tab = line.find('\t');
    if (first_tab >= line.size())
      return std::cerr << "Error: malformed probe direction map file (" << argv[3] << ")\n", EXIT_FAILURE;
    std::string status = line.substr(first_tab + 1);
    if (status != "0")
      reverse_complement_ids.insert(line.substr(0, first_tab));
  }

  if (!std::getline(methy_file, line))
    return std::cerr << "Error: empty BED file (" << argv[1] << ")\n", EXIT_FAILURE;

  output_file << line << std::endl;
  std::vector<std::string> bed_sample_ids = split_string_to_vector(line, '\t');
  if (bed_sample_ids.size() < 5)
    return std::cerr << "Error: BED file must have at least 5 columns\n", EXIT_FAILURE;

  bed_sample_ids.erase(bed_sample_ids.begin(), bed_sample_ids.begin() + 4);
 
  auto vcf_sample_ids = vcf.subset_samples({bed_sample_ids.begin(), bed_sample_ids.end()});
  if (vcf_sample_ids != bed_sample_ids)
    return std::cerr << "Error: all sample IDs in BED file must exist in VCF file and in the same order\n", EXIT_FAILURE;

  std::list<methy_t> methy;
  while (std::getline(methy_file, line))
  {
    methy.push_back(methy_t::deserialize(line, reverse_complement_ids));
    if (methy.back().signals().size() != bed_sample_ids.size())
      return std::cerr << "Error: Invalid number of columns in  BED file\n", EXIT_FAILURE;
  }

  savvy::genomic_region reg(methy.front().chrom(), 
    std::max<std::int64_t>(0, methy.front().start_pos(probe_length)),
    std::max<std::int64_t>(0, methy.back().end_pos(probe_length)));

  if (methy.front().chrom() != methy.back().chrom())
    return std::cerr << "Error: multiple chromosomes in a BED file is not yet supported\n", EXIT_FAILURE;

  vcf.reset_bounds(reg);

  
  savvy::compressed_vector<std::int8_t> gts;
  savvy::variant rec;
  while (methy.size() && vcf >> rec)
  {
    while (methy.size() && methy.front().distance(rec, probe_length) > 0)
    {
      if (method != method_t::filter || methy.front().mask_proportion() < filter_threshold)
      {
        if (method == method_t::regress)
          methy.front().regress_out_genotypes();
        else if (method == method_t::mask)
          methy.front().mask_signals();
        if (inv_norm)
          methy.front().inv_norm();
        methy_t::serialize(methy.front(), output_file, mask_string);
      }

      methy.pop_front();
    }

    for (auto it = methy.begin(); it != methy.end() && it->distance(rec, probe_length) == 0; ++it)
    {
      rec.get_format("GT", gts);
      if (method == method_t::mask || method == method_t::filter)
      {
        it->update_mask(gts);
      }
      else
      {
        it->add_predictor(gts);
      }
    }
  }

  while (methy.size())  
  {
    if (inv_norm)
      methy.front().inv_norm();
    methy_t::serialize(methy.front(), output_file, mask_string);
    methy.pop_front();
  }

  return EXIT_SUCCESS;
}
