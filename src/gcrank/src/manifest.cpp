#include "gcrank/manifest.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

#include "gcrank/detail/crc32.h"
#include "gcrank/error.h"

namespace gcrank {

namespace {

// ---- minimal JSON value --------------------------------------------------

struct JsonValue {
  enum class Type { Null, Bool, Number, String, Array, Object };

  Type type = Type::Null;
  bool boolean = false;
  bool is_integer = false;  // Number that fits in uint64 without fraction
  uint64_t integer = 0;
  double number = 0.0;
  std::string text;
  std::vector<JsonValue> array;
  // Object members are kept as parallel key/value arrays, in insertion order.
  // A std::vector<std::pair<std::string, JsonValue>> member would need the
  // pair to be complete while JsonValue is still being defined, which is not
  // portable (Clang rejects it; GCC accepts it by accident). The writer sorts
  // the keys for canonical output.
  std::vector<std::string> keys;
  std::vector<JsonValue> values;

  const JsonValue* find(std::string_view key) const {
    for (std::size_t i = 0; i < values.size(); ++i)
      if (keys[i] == key) return &values[i];
    return nullptr;
  }
  static JsonValue make_string(std::string s) {
    JsonValue v;
    v.type = Type::String;
    v.text = std::move(s);
    return v;
  }
  static JsonValue make_uint(uint64_t u) {
    JsonValue v;
    v.type = Type::Number;
    v.is_integer = true;
    v.integer = u;
    return v;
  }
  static JsonValue make_bool(bool b) {
    JsonValue v;
    v.type = Type::Bool;
    v.boolean = b;
    return v;
  }
  static JsonValue make_object() {
    JsonValue v;
    v.type = Type::Object;
    return v;
  }
  void set(const std::string& key, JsonValue value) {
    keys.push_back(key);
    values.push_back(std::move(value));
  }
};

// ---- parser --------------------------------------------------------------

class JsonParser {
public:
  JsonParser(std::string_view text, std::string origin)
      : text_(text), origin_(std::move(origin)) {}

  JsonValue parse_document() {
    skip_whitespace();
    JsonValue value = parse_value();
    skip_whitespace();
    if (pos_ != text_.size()) fail("trailing content after the top-level value");
    return value;
  }

private:
  [[noreturn]] void fail(const std::string& message) const {
    throw FormatError("malformed JSON in " + origin_ + " at byte " +
                      std::to_string(pos_) + ": " + message);
  }

  // Depth guard: an unbounded recursion could exhaust the stack on hostile
  // input. The manifest is at most three levels deep.
  struct DepthGuard {
    std::size_t& depth;
    ~DepthGuard() { --depth; }
  };

  void skip_whitespace() {
    while (pos_ < text_.size()) {
      const char c = text_[pos_];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        ++pos_;
      else
        break;
    }
  }

  char peek() {
    if (pos_ >= text_.size()) fail("unexpected end of input");
    return text_[pos_];
  }

  void expect(char c) {
    if (peek() != c) fail(std::string("expected '") + c + "'");
    ++pos_;
  }

  void expect_literal(std::string_view literal) {
    if (text_.substr(pos_, literal.size()) != literal)
      fail("expected '" + std::string(literal) + "'");
    pos_ += literal.size();
  }

  JsonValue parse_value() {
    switch (peek()) {
      case '{':
        return parse_object();
      case '[':
        return parse_array();
      case '"': {
        JsonValue v;
        v.type = JsonValue::Type::String;
        v.text = parse_string();
        return v;
      }
      case 't':
        expect_literal("true");
        return JsonValue::make_bool(true);
      case 'f':
        expect_literal("false");
        return JsonValue::make_bool(false);
      case 'n':
        expect_literal("null");
        return JsonValue();
      default:
        return parse_number();
    }
  }

  JsonValue parse_object() {
    if (depth_ >= 64) fail("nesting is too deep");
    DepthGuard guard{depth_};
    ++depth_;
    JsonValue v = JsonValue::make_object();
    expect('{');
    skip_whitespace();
    if (peek() == '}') {
      ++pos_;
      return v;
    }
    for (;;) {
      skip_whitespace();
      std::string key = parse_string();
      skip_whitespace();
      expect(':');
      skip_whitespace();
      JsonValue member = parse_value();
      v.keys.push_back(std::move(key));
      v.values.push_back(std::move(member));
      skip_whitespace();
      const char c = peek();
      if (c == ',') {
        ++pos_;
        continue;
      }
      if (c == '}') {
        ++pos_;
        return v;
      }
      fail("expected ',' or '}' in object");
    }
  }

  JsonValue parse_array() {
    if (depth_ >= 64) fail("nesting is too deep");
    DepthGuard guard{depth_};
    ++depth_;
    JsonValue v;
    v.type = JsonValue::Type::Array;
    expect('[');
    skip_whitespace();
    if (peek() == ']') {
      ++pos_;
      return v;
    }
    for (;;) {
      v.array.push_back(parse_value());
      skip_whitespace();
      const char c = peek();
      if (c == ',') {
        ++pos_;
        skip_whitespace();
        continue;
      }
      if (c == ']') {
        ++pos_;
        return v;
      }
      fail("expected ',' or ']' in array");
    }
  }

  unsigned parse_hex4() {
    if (pos_ + 4 > text_.size()) fail("truncated \\u escape");
    unsigned value = 0;
    for (int i = 0; i < 4; ++i) {
      const char c = text_[pos_++];
      unsigned digit = 0;
      if (c >= '0' && c <= '9')
        digit = static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f')
        digit = static_cast<unsigned>(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F')
        digit = static_cast<unsigned>(c - 'A' + 10);
      else
        fail("invalid hex digit in \\u escape");
      value = (value << 4) | digit;
    }
    return value;
  }

  static void append_utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
      out.push_back(static_cast<char>(cp));
    } else if (cp < 0x800) {
      out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp < 0x10000) {
      out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else {
      out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
  }

  std::string parse_string() {
    expect('"');
    std::string out;
    for (;;) {
      if (pos_ >= text_.size()) fail("unterminated string");
      const unsigned char c = static_cast<unsigned char>(text_[pos_++]);
      if (c == '"') return out;
      if (c < 0x20) fail("unescaped control character in string");
      if (c != '\\') {
        out.push_back(static_cast<char>(c));
        continue;
      }
      if (pos_ >= text_.size()) fail("unterminated escape");
      const char e = text_[pos_++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/': out.push_back('/'); break;
        case 'b': out.push_back('\b'); break;
        case 'f': out.push_back('\f'); break;
        case 'n': out.push_back('\n'); break;
        case 'r': out.push_back('\r'); break;
        case 't': out.push_back('\t'); break;
        case 'u': {
          unsigned cp = parse_hex4();
          if (cp >= 0xDC00 && cp <= 0xDFFF)
            fail("lone low surrogate in \\u escape");
          if (cp >= 0xD800 && cp <= 0xDBFF) {
            // Surrogate pair.
            if (pos_ + 1 < text_.size() && text_[pos_] == '\\' && text_[pos_ + 1] == 'u') {
              pos_ += 2;
              const unsigned low = parse_hex4();
              if (low >= 0xDC00 && low <= 0xDFFF)
                cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
              else
                fail("invalid low surrogate in \\u escape");
            } else {
              fail("missing low surrogate in \\u escape");
            }
          }
          append_utf8(out, cp);
          break;
        }
        default:
          fail("unknown escape sequence");
      }
    }
  }

  JsonValue parse_number() {
    const std::size_t start = pos_;
    if (pos_ < text_.size() && text_[pos_] == '-') ++pos_;
    if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
      fail("invalid value");
    const std::size_t first_digit = pos_;
    while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9') ++pos_;
    if (pos_ - first_digit > 1 && text_[first_digit] == '0')
      fail("leading zeros are not allowed in numbers");
    bool is_integer = true;
    if (pos_ < text_.size() && text_[pos_] == '.') {
      is_integer = false;
      ++pos_;
      if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
        fail("expected a digit after '.'");
      while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9') ++pos_;
    }
    if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E')) {
      is_integer = false;
      ++pos_;
      if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-')) ++pos_;
      if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
        fail("expected a digit in the exponent");
      while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9') ++pos_;
    }
    const std::string_view token = text_.substr(start, pos_ - start);
    JsonValue v;
    v.type = JsonValue::Type::Number;
    const char* first = token.data();
    const char* last = first + token.size();
    if (is_integer && token[0] != '-') {
      uint64_t value = 0;
      const std::from_chars_result parsed = std::from_chars(first, last, value);
      if (parsed.ec == std::errc{} && parsed.ptr == last) {
        v.is_integer = true;
        v.integer = value;
        return v;
      }
      is_integer = false;
    }
    double number = 0.0;
    const std::from_chars_result parsed =
        std::from_chars(first, last, number, std::chars_format::general);
    if (parsed.ec == std::errc::result_out_of_range)
      number = token.starts_with('-') ? -std::numeric_limits<double>::infinity()
                                      : std::numeric_limits<double>::infinity();
    v.is_integer = false;
    v.number = number;
    return v;
  }

  const std::string_view text_;
  std::string origin_;
  std::size_t pos_ = 0;
  std::size_t depth_ = 0;
};

// ---- writer --------------------------------------------------------------

void write_json_string(std::string_view s, std::string& out) {
  out.push_back('"');
  for (const char c : s) {
    const unsigned char u = static_cast<unsigned char>(c);
    switch (u) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\b': out += "\\b"; break;
      case '\f': out += "\\f"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (u < 0x20) {
          std::format_to(std::back_inserter(out), "\\u{:04x}", u);
        } else {
          out.push_back(c);  // UTF-8 passes through
        }
    }
  }
  out.push_back('"');
}

void write_json(const JsonValue& value, std::string& out, int depth) {
  const std::string pad(static_cast<std::size_t>(depth) * 2, ' ');
  const std::string pad_inner(static_cast<std::size_t>(depth + 1) * 2, ' ');
  switch (value.type) {
    case JsonValue::Type::Null:
      out += "null";
      break;
    case JsonValue::Type::Bool:
      out += value.boolean ? "true" : "false";
      break;
    case JsonValue::Type::Number:
      if (value.is_integer) {
        out += std::to_string(value.integer);
      } else {
        std::format_to(std::back_inserter(out), "{:.17g}", value.number);
      }
      break;
    case JsonValue::Type::String:
      write_json_string(value.text, out);
      break;
    case JsonValue::Type::Array: {
      if (value.array.empty()) {
        out += "[]";
        break;
      }
      out += "[\n";
      for (std::size_t i = 0; i < value.array.size(); ++i) {
        out += pad_inner;
        write_json(value.array[i], out, depth + 1);
        out += (i + 1 < value.array.size()) ? ",\n" : "\n";
      }
      out += pad;
      out += "]";
      break;
    }
    case JsonValue::Type::Object: {
      if (value.values.empty()) {
        out += "{}";
        break;
      }
      std::vector<std::size_t> order(value.values.size());
      for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
      std::ranges::sort(order, std::less<>{}, [&value](std::size_t i) -> const std::string& {
        return value.keys[i];
      });
      out += "{\n";
      for (std::size_t i = 0; i < order.size(); ++i) {
        const std::size_t at = order[i];
        out += pad_inner;
        write_json_string(value.keys[at], out);
        out += ": ";
        write_json(value.values[at], out, depth + 1);
        out += (i + 1 < order.size()) ? ",\n" : "\n";
      }
      out += pad;
      out += "}";
      break;
    }
  }
}

// ---- schema helpers ------------------------------------------------------

const JsonValue& require_field(const JsonValue& object, std::string_view key,
                               const std::string& ctx) {
  if (object.type != JsonValue::Type::Object)
    throw FormatError(ctx + ": the top-level value must be a JSON object");
  const JsonValue* value = object.find(key);
  if (value == nullptr)
    throw FormatError(ctx + ": missing required key '" + std::string(key) + "'");
  return *value;
}

std::string require_string(const JsonValue& object, std::string_view key,
                           const std::string& ctx) {
  const JsonValue& value = require_field(object, key, ctx);
  if (value.type != JsonValue::Type::String)
    throw FormatError(ctx + ": key '" + std::string(key) + "' must be a string");
  return value.text;
}

uint64_t require_uint(const JsonValue& object, std::string_view key, const std::string& ctx) {
  const JsonValue& value = require_field(object, key, ctx);
  if (value.type != JsonValue::Type::Number || !value.is_integer)
    throw FormatError(ctx + ": key '" + std::string(key) +
                      "' must be a non-negative integer");
  return value.integer;
}

const JsonValue& require_object(const JsonValue& object, std::string_view key,
                                const std::string& ctx) {
  const JsonValue& value = require_field(object, key, ctx);
  if (value.type != JsonValue::Type::Object)
    throw FormatError(ctx + ": key '" + std::string(key) + "' must be an object");
  return value;
}

bool require_bool(const JsonValue& object, std::string_view key, const std::string& ctx) {
  const JsonValue& value = require_field(object, key, ctx);
  if (value.type != JsonValue::Type::Bool)
    throw FormatError(ctx + ": key '" + std::string(key) + "' must be a boolean");
  return value.boolean;
}

uint32_t parse_crc_hex(std::string_view text, const std::string& ctx) {
  if (text.size() != 8) throw FormatError(ctx + ": crc32 must be 8 hex digits");
  uint32_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value, 16);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw FormatError(ctx + ": crc32 contains a non-hex digit");
  return value;
}

}  // namespace

std::string created_from_env() {
  std::time_t when = std::time(nullptr);
  const char* env = std::getenv("SOURCE_DATE_EPOCH");
  if (env != nullptr && *env != '\0') {
    const std::string_view text(env);
    long long value = 0;
    const std::from_chars_result parsed =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || value < 0)
      throw UsageError(std::string("SOURCE_DATE_EPOCH is not a valid epoch: '") + env + "'");
    when = static_cast<std::time_t>(value);
  }
  const auto seconds =
      std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::from_time_t(when));
  return std::format("{:%Y-%m-%dT%H:%M:%SZ}", seconds);
}

std::map<std::string, std::string> default_semantics() {
  return {
      {"ambiguous_bases",
       "N, R, Y, S, W, K, M, B, D, H, V: never GC and never valid; removed from the "
       "denominator by gc_fraction(remove_ambiguous = true)"},
      {"base_encoding", "genome.bv: 1 = unambiguous G/C; valid.bv: 1 = unambiguous A/C/G/T"},
      {"gc_bases", "GC"},
      {"remove_ambiguous",
       "when true the denominator is valid.bv (A+C+G+T), otherwise the window length"},
      {"uracil", "U is the RNA spelling of T and is folded onto T: valid, never GC"},
      {"valid_bases", "ACGT"},
      {"window", "half-open [start, end), 0-based, over the concatenated contigs"},
  };
}

// ---- Manifest ------------------------------------------------------------

std::string Manifest::to_json() const {
  JsonValue root = JsonValue::make_object();
  root.set("format_version", JsonValue::make_string(format_version));
  root.set("builder_version", JsonValue::make_string(builder_version));
  root.set("created", JsonValue::make_string(created));

  JsonValue source_json = JsonValue::make_object();
  source_json.set("path", JsonValue::make_string(source.path));
  source_json.set("sha256", JsonValue::make_string(source.sha256));
  source_json.set("file_bytes", JsonValue::make_uint(source.file_bytes));
  source_json.set("total_bases", JsonValue::make_uint(source.total_bases));
  source_json.set("num_contigs", JsonValue::make_uint(source.num_contigs));
  root.set("source", std::move(source_json));

  JsonValue index_json = JsonValue::make_object();
  index_json.set("bit_order", JsonValue::make_string(bit_order));
  index_json.set("endianness", JsonValue::make_string(endianness));
  index_json.set("superblock_bits", JsonValue::make_uint(superblock_bits));
  index_json.set("basic_block_bits", JsonValue::make_uint(basic_block_bits));
  index_json.set("remove_ambiguous_supported", JsonValue::make_bool(remove_ambiguous_supported));
  if (has_totals) {
    index_json.set("total_gc", JsonValue::make_uint(total_gc));
    index_json.set("total_valid", JsonValue::make_uint(total_valid));
  }
  root.set("index", std::move(index_json));

  JsonValue files_json = JsonValue::make_object();
  for (const auto& [name, info] : files) {
    JsonValue entry = JsonValue::make_object();
    entry.set("size", JsonValue::make_uint(info.size));
    entry.set("crc32", JsonValue::make_string(crc32_hex(info.crc32)));
    files_json.set(name, std::move(entry));
  }
  root.set("files", std::move(files_json));

  JsonValue semantics_json = JsonValue::make_object();
  for (const auto& [key, text] : semantics)
    semantics_json.set(key, JsonValue::make_string(text));
  root.set("semantics", std::move(semantics_json));

  std::string out;
  write_json(root, out, 0);
  out.push_back('\n');
  return out;
}

Manifest Manifest::from_json(std::string_view text, const std::string& origin) {
  const JsonValue root = JsonParser(text, origin).parse_document();
  const std::string ctx = origin;

  Manifest manifest;
  manifest.format_version = require_string(root, "format_version", ctx);
  const int major = parse_major_version(manifest.format_version);
  if (major != parse_major_version(kFormatVersion))
    throw FormatError(ctx + ": unsupported format_version '" + manifest.format_version +
                      "' (this build reads major version " +
                      std::to_string(parse_major_version(kFormatVersion)) + ")");
  manifest.builder_version = require_string(root, "builder_version", ctx);
  manifest.created = require_string(root, "created", ctx);

  const JsonValue& source = require_object(root, "source", ctx);
  manifest.source.path = require_string(source, "path", ctx);
  manifest.source.sha256 = require_string(source, "sha256", ctx);
  manifest.source.total_bases = require_uint(source, "total_bases", ctx);
  manifest.source.num_contigs = require_uint(source, "num_contigs", ctx);
  if (source.find("file_bytes") != nullptr)
    manifest.source.file_bytes = require_uint(source, "file_bytes", ctx);

  const JsonValue& index = require_object(root, "index", ctx);
  manifest.bit_order = require_string(index, "bit_order", ctx);
  if (manifest.bit_order != kBitOrderMsb0)
    throw FormatError(ctx + ": unsupported bit_order '" + manifest.bit_order + "'");
  manifest.endianness = require_string(index, "endianness", ctx);
  if (manifest.endianness != kEndiannessLittle && manifest.endianness != kEndiannessBig)
    throw FormatError(ctx + ": unsupported endianness '" + manifest.endianness + "'");
  manifest.superblock_bits = static_cast<uint32_t>(require_uint(index, "superblock_bits", ctx));
  manifest.basic_block_bits = static_cast<uint32_t>(require_uint(index, "basic_block_bits", ctx));
  if (manifest.superblock_bits != kSuperblockBits || manifest.basic_block_bits != kBasicBlockBits)
    throw FormatError(ctx + ": unsupported rank block sizes (" +
                      std::to_string(manifest.superblock_bits) + "/" +
                      std::to_string(manifest.basic_block_bits) + ")");
  manifest.remove_ambiguous_supported = require_bool(index, "remove_ambiguous_supported", ctx);
  if (!manifest.remove_ambiguous_supported)
    throw FormatError(ctx +
                      ": remove_ambiguous_supported = false is not supported by this build");
  if (index.find("total_gc") != nullptr) {
    manifest.total_gc = require_uint(index, "total_gc", ctx);
    manifest.total_valid = require_uint(index, "total_valid", ctx);
    manifest.has_totals = true;
    if (manifest.total_gc > manifest.source.total_bases ||
        manifest.total_valid > manifest.source.total_bases)
      throw FormatError(ctx + ": total_gc/total_valid exceed total_bases");
  } else if (index.find("total_valid") != nullptr) {
    throw FormatError(ctx + ": missing required key 'index.total_gc'");
  }

  const JsonValue& files = require_object(root, "files", ctx);
  for (const char* name : kDataFiles) {
    const JsonValue& info = require_object(files, name, ctx);
    FileInfo entry;
    entry.size = require_uint(info, "size", ctx);
    entry.crc32 = parse_crc_hex(require_string(info, "crc32", ctx), ctx);
    manifest.files[name] = entry;
  }

  const JsonValue& semantics = require_object(root, "semantics", ctx);
  for (std::size_t i = 0; i < semantics.values.size(); ++i) {
    const std::string& key = semantics.keys[i];
    const JsonValue& value = semantics.values[i];
    if (value.type != JsonValue::Type::String)
      throw FormatError(ctx + ": semantics." + key + " must be a string");
    manifest.semantics[key] = value.text;
  }
  for (const char* key : {"gc_bases", "valid_bases"})
    if (manifest.semantics.find(key) == manifest.semantics.end())
      throw FormatError(ctx + ": missing required key 'semantics." + key + "'");

  return manifest;
}

void Manifest::save(const std::string& path) const {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path + "' for writing");
  const std::string text = to_json();
  out.write(text.data(), static_cast<std::streamsize>(text.size()));
  out.flush();
  if (!out) throw OutputError("write error on '" + path + "'");
}

Manifest Manifest::load(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  if (in.bad()) throw IoError("read error on '" + path + "'");
  return from_json(text, path);
}

}  // namespace gcrank
