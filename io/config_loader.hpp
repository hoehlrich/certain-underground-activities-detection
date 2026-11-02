#pragma once

#include <filesystem>

#include "core/config.hpp"

namespace cfar::io {

// Reads the processing sections ([band], [acq], [false_alarm], [cfar]) of a TOML
// config into the core structs; other sections are ignored. Calls
// cfar::validate(). Throws std::runtime_error on a parse error or a bad value.
ProcessingConfig load_processing_config(const std::filesystem::path& path);

}  // namespace cfar::io
