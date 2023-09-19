/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>
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

#include <xtensor/xtensor.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor-blas/xlinalg.hpp>

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
  std::string vcf_path_;
  std::string methy_path_;
  std::string manifest_path_;
  std::string output_path_ = "/dev/stdout";
  std::string mask_code_ = "NA";
  double filter_threhold_ = 0.05;
  method_t method_ = method_t::regress;
  bool inv_norm_ = false;
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
        {"filter-threshold", "<float>", 'p', "Probes will be excluded when proportion of samples with variation under the probe >= this value (--method=filter only; default: 0.05)"},
        {"mask-code", "<string>", 'c', "Character sequence used to denote missing/masked values when using --method=mask (default: NA)"},
        {"inv-norm", "", 'i', "Inverse-normalize output"},
      })
  {
  }

  const std::string& vcf_path() const { return vcf_path_; }
  const std::string& methy_path() const { return methy_path_; }
  const std::string& manifest_path() const { return manifest_path_; }
  const std::string& output_path() const { return output_path_; }

  const std::string& mask_code() const { return mask_code_; }
  double filter_threshold() const { return filter_threhold_; }
  method_t method() const { return method_; }
  bool inv_norm() const { return inv_norm_; }
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
      case 'h':
        help_ = true;
        return true;
      case 'v':
        version_ = true;
        return true;
      case 'c':
        mask_code_ = optarg ? optarg : "";
        break;
      case 'i':
        inv_norm_ = true;
        break;
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

class manifest_entry
{
public:
  std::string chrom;
  std::int64_t cpg_beg_0 = 0;
  short map_flag = 0;
  short type = 0;
  manifest_entry() {}
  manifest_entry(std::string c, std::int64_t p, short m, short t) :
    chrom(std::move(c)), cpg_beg_0(p), map_flag(m), type(t)
  {}
};

class methy_t
{
private:
  std::string chrom_;
  std::int64_t start_;
  std::int64_t end_;
  std::string id_;
  std::vector<signal_t> signals_;
  std::int64_t affected_region_start_;
  std::int64_t affected_region_end_;
  std::vector<bool> mask_;
  std::list<savvy::compressed_vector<std::int8_t>> predictor_list_;
public:
  const std::string& chrom() const { return chrom_; }
  const std::vector<signal_t>& signals() const { return signals_; }
  const std::int64_t start_pos() const { return affected_region_start_; }
  const std::int64_t end_pos() const { return affected_region_end_; }

  const std::int64_t distance(const savvy::site_info& var, std::size_t alt_idx) const
  {
    auto var_length = std::max<std::int64_t>(var.alts()[alt_idx].size(), var.ref().size());
    std::int64_t s = start_pos();
    std::int64_t e = end_pos();
    if (var.pos() + var_length - 1 < s)
      return (var.pos() + var_length - 1) - s;
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
    assert(gts.size() == signals_.size());
    predictor_list_.push_back(gts);
  }

  void regress_out_genotypes_additive()
  {
    using namespace xt;
    using namespace xt::linalg;

    if (predictor_list_.size())
    {
      auto y = xt::adapt(signals_, {signals_.size()});
      xtensor<double, 2> x = zeros<double>({predictor_list_.size() + 1, predictor_list_.front().size()});
      auto it = predictor_list_.begin();
      for (std::size_t i = 0; i < predictor_list_.size(); ++i,++it)
      {
        assert(it->size() == x.shape(1));
        for (auto jt = it->begin(); jt != it->end(); ++jt)
        {
          x(i, jt.offset()) = *jt;
        }
      }
      xt::row(x, predictor_list_.size()) = xt::ones<double>({x.shape(1)});

      auto pbetas = dot(dot(pinv(dot(x, transpose(x))), x), y);
      y -= dot(transpose(x), pbetas);
    }

  }

  void regress_out_genotypes()
  {
    using namespace xt;
    using namespace xt::linalg;

    if (predictor_list_.size())
    {
      auto y = xt::adapt(signals_, {signals_.size()});
      xtensor<double, 2> x = zeros<double>({predictor_list_.size() * 2 + 1, predictor_list_.front().size()});
      auto it = predictor_list_.begin();
      for (std::size_t i = 0; i < predictor_list_.size(); ++i,++it)
      {
        assert(it->size() == x.shape(1));
        for (auto jt = it->begin(); jt != it->end(); ++jt)
        {
          x(i * 2 + int(*jt == 2), jt.offset()) = 1;
        }
      }
      xt::row(x, predictor_list_.size() * 2) = xt::ones<double>({x.shape(1)});

      auto pbetas = dot(dot(pinv(dot(x, transpose(x))), x), y);
      y -= dot(transpose(x), pbetas);
    }

  }
  
  void inv_norm()
  {
    inverse_normalize(signals_);
  }

  static methy_t deserialize(const std::string& line, const std::unordered_map<std::string, manifest_entry>& manifest)
  {
    methy_t ret;
    auto vec = split_string_to_vector(line, '\t');
    if (vec.size() >= 5)
    {
      ret.chrom_ = vec[0];
      ret.start_ = std::atoll(vec[1].c_str());
      ret.end_ = std::atoll(vec[2].c_str());
      ret.id_ = vec[3];
      auto res = manifest.find(ret.id_);
      if (res != manifest.end())
      {
        if (res->second.map_flag == 0)
        {
          ret.affected_region_end_ = res->second.cpg_beg_0 + 2;
          ret.affected_region_start_ = ret.affected_region_end_ - 50 + 1; // - 50pb probe + 1 for inclusive end
          if (res->second.type == 1)
            ++(ret.affected_region_end_); // add 1 extension base
          else
            --(ret.affected_region_start_); // shift 1 for extension base
        }
        else
        {
          ret.affected_region_start_ = res->second.cpg_beg_0 + 1;
          ret.affected_region_end_ = ret.affected_region_start_ + 50 - 1; // + 50pb probe - 1 for inclusive end
          if (res->second.type == 1)
            --(ret.affected_region_start_); // subtract 1 extension base
          else
            ++(ret.affected_region_end_); // shift 1 for extension base
        }
      }

      ret.signals_.reserve(vec.size() - 4);
      for (auto it = vec.begin() + 4; it != vec.end(); ++it)
        ret.signals_.push_back(std::atof(it->c_str()));
      ret.mask_.resize(ret.signals_.size(), false);
    }
    return ret;
  }

  static bool serialize(const methy_t& m, std::ostream& ofs, const std::string& mask_str)
  {
    ofs << m.chrom_ << "\t" << m.start_ << "\t" << m.end_ << "\t" << m.id_;
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

std::unordered_map<std::string, manifest_entry> parse_epic_manifest(const std::string& file_path)
{
  shrinkwrap::istream manifest_file(file_path);
  std::string line;
  std::unordered_map<std::string, manifest_entry> manifest;
  int id_col = -1, chrom_col = -1, cpg_beg_col = -1, map_flag_col = -1, type_col = -1;
  int min_cols = 5;
  while (std::getline(manifest_file, line))
  {
    auto fields = split_string_to_vector(line, '\t');
    if (fields.size() < min_cols)
    {
      return std::cerr << "Error: malformed EPIC manifest file\n", manifest;
    }
    else
    {
      if (id_col < 0)
      {
        std::unordered_map<std::string, int> col_map;
        for (std::size_t i = 0; i < fields.size(); ++i)
          col_map[fields[i]] = i + 1;

        id_col = col_map["Probe_ID"];
        if (id_col == 0)
          return std::cerr << "Error: Probe_ID missing from manifest\n", manifest;
        chrom_col = col_map["CpG_chrm"];
        if (chrom_col == 0)
          return std::cerr << "Error: CpG_chrm missing from manifest\n", manifest;
        cpg_beg_col = col_map["CpG_beg"];
        if (cpg_beg_col == 0)
          return std::cerr << "Error: CpG_beg missing from manifest\n", manifest;
        map_flag_col = col_map["mapFlag_A"];
        if (map_flag_col == 0)
          return std::cerr << "Error: mapFlag_A missing from manifest\n", manifest;
        type_col = col_map["type"];
        if (type_col == 0)
          return std::cerr << "Error: type missing from manifest\n", manifest;

        min_cols = std::max(min_cols, id_col);
        min_cols = std::max(min_cols, chrom_col);
        min_cols = std::max(min_cols, cpg_beg_col);
        min_cols = std::max(min_cols, map_flag_col);
        min_cols = std::max(min_cols, type_col);
        --id_col; --chrom_col; --cpg_beg_col; --map_flag_col; --type_col;
      }
      else
      {
        manifest_entry e;
        e.chrom = fields[chrom_col];
        e.cpg_beg_0 = std::atoll(fields[cpg_beg_col].c_str());
        e.map_flag = std::atol(fields[map_flag_col].c_str());
        e.type = fields[type_col] == "II" ? 2 : 1;
        manifest.emplace(fields[id_col], std::move(e));
      }
    }
  }

  return manifest;
}

void write_methy(methy_t& m, std::ostream& output_file, const vuptool_args& args)
{
  if (args.method() != method_t::filter || m.mask_proportion() < args.filter_threshold())
  {
    if (args.method() == method_t::regress)
      m.regress_out_genotypes();
    else if (args.method() == method_t::mask)
      m.mask_signals();
    if (args.inv_norm())
      m.inv_norm();
    methy_t::serialize(m, output_file, args.mask_code());
  }
}

int main(int argc, char** argv)
{

  using namespace xt;
  using namespace xt::linalg;
  xtensor<double, 1> y = {0.2, 0.8, 0.9, 0.85};
  xtensor<double, 2> X = {
    {1., 0.1, 0.3, 0.2, 0.002, 0.15},
    {1., 0.7, 0.7, 0.5, 0.77, 0.59},
    {1., 0.5, 0.72, 0.6, 0.71, 0.7},
    {1., 0.9, 0.71, 0.66, 0.68, 0.8}
  };
  double lambda = 1.;
  auto I = diag(ones<double>({X.shape(1)}));
  xtensor<double, 1> beta = dot(pinv(dot(transpose(X), X) + lambda * I), dot(transpose(X), y));

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

  auto manifest = parse_epic_manifest(args.manifest_path());
  if (manifest.empty())
    return std::cerr << "Error: failed to parse EPIC array manifest\n", EXIT_FAILURE;

  shrinkwrap::istream methy_file(args.methy_path());

  savvy::reader vcf(args.vcf_path());
  if (!vcf)
    return std::cerr << "Error: failed to open VCF file\n", EXIT_FAILURE;
  std::ofstream output_file(args.output_path(), std::ios::binary);

  std::string line;
  if (!std::getline(methy_file, line))
    return std::cerr << "Error: empty BED file\n", EXIT_FAILURE;

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
    methy.push_back(methy_t::deserialize(line, manifest));
    if (methy.back().signals().size() != bed_sample_ids.size())
      return std::cerr << "Error: Invalid number of columns in  BED file\n", EXIT_FAILURE;
  }

  savvy::genomic_region reg(methy.front().chrom(),
    std::max<std::int64_t>(0, methy.front().start_pos()),
    std::max<std::int64_t>(0, methy.back().end_pos()));

  if (methy.front().chrom() != methy.back().chrom())
    return std::cerr << "Error: multiple chromosomes in a BED file is not yet supported\n", EXIT_FAILURE;

  vcf.reset_bounds(reg);
  if (!vcf)
    return std::cerr << "Error: could not seek to region " << reg.chromosome() << ":" << reg.from() << "-"  << reg.to() << " in VCF\n", EXIT_FAILURE;

  
  savvy::compressed_vector<std::int8_t> gts;
  savvy::variant rec;
  while (methy.size() && vcf >> rec)
  {
    if (rec.alts().size() != 1)
      return std::cerr << "Error: all variant records must be biallelic\n", EXIT_FAILURE;

    while (methy.size() && methy.front().distance(rec, 0) > 0)
    {
      write_methy(methy.front(), output_file, args);
      methy.pop_front();
    }

    for (auto it = methy.begin(); it != methy.end() && it->distance(rec, 0) == 0; ++it)
    {
      rec.get_format("GT", gts);
      if (args.method() == method_t::mask || args.method() == method_t::filter)
      {
        it->update_mask(gts);
      }
      else
      {
        auto stride = gts.size() / vcf.samples().size();
        if (stride > 2)
          return std::cerr << "Error: ploidy greater than 2 is currently not supported with regress method.\n", EXIT_FAILURE;
        savvy::stride_reduce(gts, stride);
        it->add_predictor(gts);
      }
    }
  }

  while (methy.size())  
  {
    write_methy(methy.front(), output_file, args);
    methy.pop_front();
  }

  return output_file.good() && !vcf.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}
