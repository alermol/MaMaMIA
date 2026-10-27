#include "gcrank/detail/sha256.h"

#include <array>
#include <bit>
#include <cstring>
#include <fstream>
#include <span>

// GCRANK_DISABLE_SHA_NI exists so the scalar fallback can be built and tested
// on a machine that supports the extension.
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(GCRANK_DISABLE_SHA_NI)
#include <cpuid.h>
#include <immintrin.h>
#define GCRANK_HAVE_SHA_NI 1
#endif

#include "gcrank/error.h"

namespace gcrank {

namespace {

const uint32_t kK[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u,
};

constexpr char kHexDigits[] = "0123456789abcdef";
constexpr std::size_t kBufferBytes = 1u << 18;

// Lowercase hex rendering of a digest.
std::string hex_encode(std::span<const uint8_t> digest) {
  std::string out(digest.size() * 2, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[2 * i] = kHexDigits[digest[i] >> 4];
    out[2 * i + 1] = kHexDigits[digest[i] & 0x0F];
  }
  return out;
}

}  // namespace

void SHA256::reset() {
  h_[0] = 0x6a09e667u;
  h_[1] = 0xbb67ae85u;
  h_[2] = 0x3c6ef372u;
  h_[3] = 0xa54ff53au;
  h_[4] = 0x510e527fu;
  h_[5] = 0x9b05688cu;
  h_[6] = 0x1f83d9abu;
  h_[7] = 0x5be0cd19u;
  buf_len_ = 0;
  total_bits_ = 0;
}

namespace {

void transform_scalar(uint32_t state[8], const uint8_t block[64]) {
  uint32_t w[64];
  for (int i = 0; i < 16; ++i) {
    w[i] = (static_cast<uint32_t>(block[4 * i]) << 24) |
           (static_cast<uint32_t>(block[4 * i + 1]) << 16) |
           (static_cast<uint32_t>(block[4 * i + 2]) << 8) |
           static_cast<uint32_t>(block[4 * i + 3]);
  }
  for (int i = 16; i < 64; ++i) {
    const uint32_t s0 = std::rotr(w[i - 15], 7) ^ std::rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
    const uint32_t s1 = std::rotr(w[i - 2], 17) ^ std::rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }

  uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
  uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
  for (int i = 0; i < 64; ++i) {
    const uint32_t s1 = std::rotr(e, 6) ^ std::rotr(e, 11) ^ std::rotr(e, 25);
    const uint32_t ch = (e & f) ^ (~e & g);
    const uint32_t t1 = h + s1 + ch + kK[i] + w[i];
    const uint32_t s0 = std::rotr(a, 2) ^ std::rotr(a, 13) ^ std::rotr(a, 22);
    const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
    const uint32_t t2 = s0 + maj;
    h = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }
  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;
  state[4] += e;
  state[5] += f;
  state[6] += g;
  state[7] += h;
}

#if defined(GCRANK_HAVE_SHA_NI)
// SHA-256 using the x86 SHA extensions (SHA-NI). The transformation between
// the A-H state and the ABEF/CDGH layout expected by sha256rnds2 is the
// standard one from Intel's SHA extension guide; the message schedule uses
// sha256msg1/sha256msg2, unrolled because the schedule registers must stay in
// the register file.
__attribute__((target("sha,sse4.1,ssse3"))) void transform_sha_ni(uint32_t state[8],
                                                                  const uint8_t block[64]) {
  __m128i state0, state1, msg, tmp;
  __m128i m0, m1, m2, m3;
  const __m128i mask =
      _mm_set_epi64x(0x0c0d0e0f08090a0bULL, 0x0405060700010203ULL);

  // State layout: [A B C D] and [E F G H] -> ABEF and CDGH.
  tmp = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&state[0]));
  state1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&state[4]));
  tmp = _mm_shuffle_epi32(tmp, 0xB1);
  state1 = _mm_shuffle_epi32(state1, 0x1B);
  state0 = _mm_alignr_epi8(tmp, state1, 8);
  state1 = _mm_blend_epi16(state1, tmp, 0xF0);
  const __m128i save0 = state0;
  const __m128i save1 = state1;

  m0 = _mm_shuffle_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i*>(block + 0)), mask);
  m1 = _mm_shuffle_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i*>(block + 16)), mask);
  m2 = _mm_shuffle_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i*>(block + 32)), mask);
  m3 = _mm_shuffle_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i*>(block + 48)), mask);

#if defined(__clang__)
#pragma clang loop unroll(full)
#elif defined(__GNUC__)
#pragma GCC unroll 16
#endif
  for (int g = 0; g < 16; ++g) {
    // The four message registers hold groups g%4, (g+1)%4, ...; a runtime
    // index is fine here because the loop is unrolled and they become
    // compile-time constants.
    __m128i* const regs[4] = {&m0, &m1, &m2, &m3};
    __m128i& cur = *regs[g & 3];
    msg = _mm_add_epi32(
        cur, _mm_loadu_si128(reinterpret_cast<const __m128i*>(kK + 4 * g)));
    state1 = _mm_sha256rnds2_epu32(state1, state0, msg);
    msg = _mm_shuffle_epi32(msg, 0x0E);
    state0 = _mm_sha256rnds2_epu32(state0, state1, msg);

    if (g == 1) {
      m0 = _mm_sha256msg1_epu32(m0, m1);
    } else if (g == 2) {
      m1 = _mm_sha256msg1_epu32(m1, m2);
    } else if (g >= 3) {
      __m128i& next = *regs[(g + 1) & 3];
      __m128i& prev = *regs[(g + 3) & 3];
      tmp = _mm_alignr_epi8(cur, prev, 4);
      next = _mm_add_epi32(next, tmp);
      next = _mm_sha256msg2_epu32(next, cur);
      prev = _mm_sha256msg1_epu32(prev, cur);
    }
  }

  state0 = _mm_add_epi32(state0, save0);
  state1 = _mm_add_epi32(state1, save1);

  // ABEF and CDGH -> [A B C D] and [E F G H].
  tmp = _mm_shuffle_epi32(state0, 0x1B);
  state1 = _mm_shuffle_epi32(state1, 0xB1);
  state0 = _mm_blend_epi16(tmp, state1, 0xF0);
  state1 = _mm_alignr_epi8(state1, tmp, 8);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&state[0]), state0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&state[4]), state1);
}

using TransformFn = void (*)(uint32_t*, const uint8_t*);

// CPUID, rather than __builtin_cpu_supports: Clang 18 and earlier do not accept
// any feature string for SHA, while the CPUID bits are the same everywhere.
// SHA-NI implies SSE4.1, but check the prerequisites explicitly before using
// the target("sha,sse4.1,ssse3") function.
bool sha_ni_available() {
  unsigned eax = 0;
  unsigned ebx = 0;
  unsigned ecx = 0;
  unsigned edx = 0;
  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx) == 0) return false;
  constexpr unsigned kSsse3 = 1u << 9;
  constexpr unsigned kSse41 = 1u << 19;
  if ((ecx & (kSsse3 | kSse41)) != (kSsse3 | kSse41)) return false;
  if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx) == 0) return false;
  constexpr unsigned kSha = 1u << 29;
  return (ebx & kSha) != 0;
}

TransformFn select_transform() {
  return sha_ni_available() ? transform_sha_ni : transform_scalar;
}
#else
using TransformFn = void (*)(uint32_t*, const uint8_t*);
TransformFn select_transform() { return transform_scalar; }
#endif

}  // namespace

void SHA256::transform(const uint8_t block[64]) {
  // Resolved once on first use; static local initialization is thread-safe.
  static const TransformFn transform_fn = select_transform();
  transform_fn(h_, block);
}

void SHA256::update(const void* data, std::size_t n) {
  const uint8_t* p = static_cast<const uint8_t*>(data);
  total_bits_ += static_cast<uint64_t>(n) * 8u;
  if (buf_len_ > 0) {
    const std::size_t want = 64 - buf_len_;
    const std::size_t take = n < want ? n : want;
    std::memcpy(buf_ + buf_len_, p, take);
    buf_len_ += take;
    p += take;
    n -= take;
    if (buf_len_ == 64) {
      transform(buf_);
      buf_len_ = 0;
    }
  }
  while (n >= 64) {
    transform(p);
    p += 64;
    n -= 64;
  }
  if (n > 0) {
    std::memcpy(buf_, p, n);
    buf_len_ = n;
  }
}

void SHA256::finish(uint8_t out[32]) {
  const uint64_t bits = total_bits_;
  // Padding: 0x80, zeros, 8-byte big-endian length.
  uint8_t pad[72];
  std::memset(pad, 0, sizeof(pad));
  pad[0] = 0x80;
  const std::size_t rem = buf_len_ % 64;
  const std::size_t pad_len = (rem < 56) ? (56 - rem) : (120 - rem);
  update(pad, pad_len);
  uint8_t lenbuf[8];
  for (int i = 0; i < 8; ++i) lenbuf[7 - i] = static_cast<uint8_t>(bits >> (8 * i));
  update(lenbuf, 8);
  for (int i = 0; i < 8; ++i) {
    out[4 * i] = static_cast<uint8_t>(h_[i] >> 24);
    out[4 * i + 1] = static_cast<uint8_t>(h_[i] >> 16);
    out[4 * i + 2] = static_cast<uint8_t>(h_[i] >> 8);
    out[4 * i + 3] = static_cast<uint8_t>(h_[i]);
  }
}

std::string sha256_hex(const void* data, std::size_t n) {
  SHA256 h;
  h.update(data, n);
  uint8_t digest[32];
  h.finish(digest);
  return hex_encode(digest);
}

std::string sha256_hex_from_digest(const uint8_t digest[32]) {
  return hex_encode(std::span<const uint8_t>(digest, 32));
}

std::string sha256_file_hex(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  SHA256 h;
  std::array<char, kBufferBytes> buf;
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got > 0) h.update(buf.data(), static_cast<std::size_t>(got));
  }
  if (in.bad()) throw IoError("read error on '" + path + "'");
  uint8_t digest[32];
  h.finish(digest);
  return hex_encode(digest);
}

}  // namespace gcrank
