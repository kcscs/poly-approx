#pragma once
#include <memory>
#include "types.hpp"

template <typename FT>
struct StopStrategy;
template <typename FT>
std::unique_ptr<StopStrategy<FT>> CreateStopStrategy(json config);


template <typename FT>
struct StopStrategy {
    virtual bool stop(FT cur_error) = 0;
    virtual void reset() {}

    virtual ~StopStrategy() = default;

protected:
    StopStrategy() = default;
    friend std::unique_ptr<StopStrategy<FT>> CreateStopStrategy<>(json config);
};




template <typename FT>
struct StopIfErrorChangeIsBelowThreshold : public StopStrategy<FT> {
    virtual bool stop(FT cur_error) override {
        if(last_error && abs(last_error.value() - cur_error) < stop_threshold)
            return true;
        last_error = cur_error;
        return false;
    }

    virtual void reset() override {
        last_error = std::nullopt;
    }

    virtual ~StopIfErrorChangeIsBelowThreshold() = default;

protected:
    FT stop_threshold;
    std::optional<FT> last_error;

    StopIfErrorChangeIsBelowThreshold(FT stop_threshold) : stop_threshold(stop_threshold) {}
    friend std::unique_ptr<StopStrategy<FT>> CreateStopStrategy<>(json config);
};




template <typename FT>
struct StopIfErrorIsBelowThreshold : public StopStrategy<FT> {
    virtual bool stop(FT cur_error) override {
        return abs(cur_error) < stop_threshold;
    }

    virtual ~StopIfErrorIsBelowThreshold() = default; 

protected:
    FT stop_threshold;

    StopIfErrorIsBelowThreshold(FT stop_threshold) : stop_threshold(stop_threshold) {}
    friend std::unique_ptr<StopStrategy<FT>> CreateStopStrategy<>(json config);
};






// ---------- Factory ------------

template <typename FT>
std::unique_ptr<StopStrategy<FT>> CreateStopStrategy(json config) {
    std::string strat_name = config.value("stop_strategy", "threshold");
    FT threshold = config.value("stop_threshold", config.value("target_precision", 1e-5f));
    if (strat_name == "threshold") {
        return std::unique_ptr<StopIfErrorIsBelowThreshold<FT>>(new StopIfErrorIsBelowThreshold<FT>(threshold));
    } else if(strat_name == "error_change_threshold") {
        return std::unique_ptr<StopIfErrorChangeIsBelowThreshold<FT>>(new StopIfErrorChangeIsBelowThreshold<FT>(threshold));
    } else {
        throw std::logic_error("Unknown stop strategy: " + strat_name);
    }
}

