#pragma once
#include <any>
#include <memory>
#include <optional>
#include <stdexcept>

#include "types.hpp"

template <typename FT> struct SplitStrategy;

template <typename FT>
std::unique_ptr<SplitStrategy<FT>> CreateSplitStrategy(json config);

template <typename FT> struct SplitContext {
  std::optional<FT> max_error_place;
};

template <typename FT> struct SplitStrategy {
  virtual FT split(FT begin, FT end, SplitContext<FT> params) = 0;

  virtual ~SplitStrategy() = default;

protected:
  SplitStrategy() {}
  friend std::unique_ptr<SplitStrategy<FT>>
  CreateSplitStrategy<FT>(json config);
};

template <typename FT> struct SplitAtHalf : SplitStrategy<FT> {
  FT split(FT begin, FT end, SplitContext<FT> params) override {
    return (begin + end) / 2;
  }

protected:
  SplitAtHalf() = default;
  friend std::unique_ptr<SplitStrategy<FT>> CreateSplitStrategy<FT>(json config);
};

template <typename FT> struct SplitAtMaxError : SplitStrategy<FT> {
  FT split(FT begin, FT end, SplitContext<FT> params) override {
    return params.max_error_place.value();
  }

protected:
  SplitAtMaxError() = default;
  friend std::unique_ptr<SplitStrategy<FT>>
  CreateSplitStrategy<FT>(json config);
};

// Factory
template <typename FT>
std::unique_ptr<SplitStrategy<FT>> CreateSplitStrategy(json config) {
  std::string strat_name = config.value("split_strategy", "split_at_half");
  if (strat_name == "split_at_half") {
    return std::unique_ptr<SplitAtHalf<FT>>(new SplitAtHalf<FT>());
  } else if (strat_name == "split_at_max_error") {
    return std::unique_ptr<SplitAtMaxError<FT>>(new SplitAtMaxError<FT>());
  } else {
    throw std::logic_error("Unknown split strategy: " + strat_name);
  }
}
