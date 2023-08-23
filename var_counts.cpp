#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <limits>
#include <cmath>

#include <savvy/reader.hpp>
#include <shrinkwrap/istream.hpp>

#include "src/inv_norm.hpp"

typedef float signal_t;

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

class methy_t
{
private:
  std::string chrom_;
  std::int64_t pos_;
  std::int64_t pos_end_;
  std::string id_;
  std::vector<signal_t> signals_;
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

  void mask_signal(std::size_t idx) { signals_[idx] = -std::numeric_limits<signal_t>::infinity(); }
  
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
    ofs.put('\n');
  }
};

int main(int argc, char** argv)
{

  std::string mask_string = "NA";
  bool inv_norm = true;
  std::int64_t probe_length = 50;

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
 
  //auto vcf_sample_ids = vcf.subset_samples({bed_sample_ids.begin(), bed_sample_ids.end()});
  //if (vcf_sample_ids != bed_sample_ids)
  //  return std::cerr << "Error: all sample IDs in BED file must exist in VCF file and in the same order\n", EXIT_FAILURE;

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

  std::vector<std::size_t> variant_counts_per_probe(methy.size());

  
  savvy::compressed_vector<std::int8_t> gts;
  savvy::variant rec;
  while (methy.size() && vcf >> rec)
  {
    while (methy.size() && methy.front().distance(rec, probe_length) > 0)
    {
      if (inv_norm)
        methy.front().inv_norm();
      //methy_t::serialize(methy.front(), output_file, mask_string);
      methy.pop_front();
    }

    std::size_t count_idx = 0;
    for (auto it = methy.begin(); it != methy.end() && it->distance(rec, probe_length) == 0; ++it,++count_idx)
    {
      ++variant_counts_per_probe[(variant_counts_per_probe.size() - methy.size()) + count_idx];
      /*rec.get_format("GT", gts);
      std::size_t stride = gts.size() / it->signals().size();
      for (auto a = gts.begin(); a != gts.end(); ++a)
      {
        it->mask_signal(a.offset() / stride);
      }*/
    }
  }

  while (methy.size())  
  {
    if (inv_norm)
      methy.front().inv_norm();
    //methy_t::serialize(methy.front(), output_file, mask_string);
    methy.pop_front();
  }

  for (std::size_t i = 0; i < variant_counts_per_probe.size(); ++i)
    output_file << variant_counts_per_probe[i] << "\n";

  return EXIT_SUCCESS;
}
