/*

   calin/io/resequenced_data_source.hpp -- Stephen Fegan -- 2019-06-06

   Resequenced(RandomAccess)DataSource - resequence data from a number of
   (RandomAccess)DataSource

   Copyright 2019, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#pragma once

#include <util/log.hpp>
#include <io/data_source.hpp>
#include <io/buffered_data_source.hpp>

namespace calin { namespace io { namespace data_source {

template<typename T> class Sequencer
{
public:
  virtual ~Sequencer() { /* nothing to see here */ }
  virtual bool is_ordered(const T& a, const T& b) = 0;
};

template<typename DST> class BasicResequencedDataSource: public DST
{
public:
  CALIN_TYPEALIAS(data_type, typename DST::data_type);

  BasicResequencedDataSource(
      Sequencer<data_type>* sequencer, std::vector<DST*> sources,
      bool adopt_sequencer = false, bool adopt_sources = false):
    DST(), sequencer_(sequencer), sources_(sources),
    adopt_sequencer_(adopt_sequencer), adopt_sources_(adopt_sources)
  {
    // nothing to see here
  }

  ~BasicResequencedDataSource()
  {
    // nothing to see here
  }

  data_type* get_next(uint64_t& seq_index_out,
    google::protobuf::Arena** arena = nullptr) override
  {
    Payload<data_type> the_next = next_[0];
    for(unsigned inext=1; index<next_.size(); inext++) {
      if(next_[inext].ptr != nullptr) {
        if(the_next.ptr == nullptr or
            sequencer_->is_ordered(*next_[inext].ptr, *the_next.ptr))
          the_next = next_[inext];
      }
    }
  }

protected:
  Sequencer<data_type>* sequencer_ = nullptr;
  std::vector<DST*> sources_;
  std::vector<Payload<data_type> > next_;
};

#ifndef SWIG
template<typename T> using ResequencedDataSource =
  BasicResequencedDataSource<DataSource<T> >;
#endif

#if 0
template<typename RADST> class BasicChainedRandomAccessDataSource:
  public BasicChainedDataSource<RADST>
{
public:
  CALIN_TYPEALIAS(data_type, typename RADST::data_type);

  BasicChainedRandomAccessDataSource(DataSourceOpener<RADST>* opener,
    bool adopt_opener = false):
    BasicChainedDataSource<RADST>(opener, adopt_opener)
  {
    // Base class can't call virtual open_file function so we complete it
    if(source_)add_source_index(source_);
  }

  virtual ~BasicChainedRandomAccessDataSource() {
    /* nothing to see here */ }

  virtual uint64_t size() override
  {
    for(unsigned i = chained_file_index_.size(); i<opener_->num_sources(); i++)
    {
      std::unique_ptr<RADST> source(opener_->open(i));
      // Throw exception in case opener doesn't do so itself
      if(!source)throw std::runtime_error("Could not open file " +
        std::to_string(i));
      add_source_index(source.get());
    }

    if(opener_->num_sources())return chained_file_index_.back();
    else return 0;
  }

  void set_next_index(uint64_t next_index) override
  {
    if(opener_->num_sources() == 0)return;

    if(next_index < chained_file_index_.back())
    {
      // next_index belongs to one of the files we have already indexed
      auto file_index = std::upper_bound(chained_file_index_.begin(),
        chained_file_index_.end(), next_index);
      unsigned new_isource = int(file_index-chained_file_index_.begin());
      if(new_isource != isource_)
      {
        isource_ = new_isource;
        BasicChainedDataSource<RADST>::open_file();
      }
    }
    else
    {
      while(isource_<opener_->num_sources() and
        chained_file_index_.back()<=next_index)
      {
        ++isource_;
        open_file();
      }
    }

    if(source_)
    {
      seq_index_ = next_index;
      if(isource_ == 0)
        source_->set_next_index(next_index);
      else
        source_->set_next_index(next_index - chained_file_index_[isource_-1]);
    }
  }

protected:
  using BasicChainedDataSource<RADST>::source_;
  using BasicChainedDataSource<RADST>::isource_;
  using BasicChainedDataSource<RADST>::opener_;
  using BasicChainedDataSource<RADST>::seq_index_;

  void add_source_index(RADST* source)
  {
    if(chained_file_index_.empty())
      chained_file_index_.push_back(source->size());
    else
      chained_file_index_.push_back(source->size()+chained_file_index_.back());
  }

  void open_file() override
  {
    BasicChainedDataSource<RADST>::open_file();
    if(source_ and isource_ >= chained_file_index_.size())
    {
      add_source_index(source_);
    }
  }

  std::vector<uint64_t> chained_file_index_;
};
#endif
} } } // namespace calin::io::data_source
