// bam-slim: remove BAM fields that MaMaMIA never reads.
//
// Kept:
//   - @HD, @SQ and @RG header lines (SO:coordinate is required by MaMaMIA)
//   - FLAG, RNAME, POS, MAPQ and CIGAR on every record
//   - the RG:Z tag unless --drop-rg is used (needed for multi-sample BAMs)
//
// Removed by default:
//   - QUAL
//   - read names (written as "*")
//   - every optional tag except RG
//   - @PG and @CO header lines
//
// SEQ is kept by default so the result remains useful to other BAM tools.
// --drop-seq also removes SEQ; MaMaMIA infers read length from CIGAR, and the
// resulting records still pass samtools quickcheck.
//
// Build:
//   g++ -O2 -std=c++23 -o bam-slim tools/bam-slim.cpp
//       $(pkg-config --cflags --libs htslib)
//
// Usage:
//   bam-slim [options] <in.bam> <out.bam>

#include <htslib/sam.h>

#include <getopt.h>
#include <sys/stat.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Options {
    bool drop_seq = false;
    bool keep_qual = false;
    bool keep_name = false;
    bool drop_rg = false;
    bool index = false;
    int threads = 2;
};

void usage(FILE *out) {
    std::fprintf(
        out,
        "Usage: bam-slim [options] <in.bam> <out.bam>\n"
        "\n"
        "Options:\n"
        "  --drop-seq       Also remove SEQ (read length is inferred from CIGAR)\n"
        "  --keep-qual      Keep base qualities (requires SEQ)\n"
        "  --keep-name      Keep read names (default writes '*')\n"
        "  --drop-rg        Drop RG tags (only safe for single-sample BAMs)\n"
        "  -@, --threads N  htslib I/O threads (default 2)\n"
        "  --index          Write a CSI index next to the output BAM\n"
        "  -h, --help       Show this help\n"
        "\n"
        "By default QUAL, read names and all optional tags except RG are\n"
        "removed; SEQ is kept. @PG and @CO header lines are always removed.\n");
}

bool header_tag(const char *line, size_t len, const char *tag) {
    const size_t n = std::strlen(tag);
    if (len < n || std::strncmp(line, tag, n) != 0) {
        return false;
    }
    return len == n || line[n] == '\t';
}

std::string slim_header(sam_hdr_t *header) {
    const char *text = sam_hdr_str(header);
    std::string out;
    if (text == nullptr) {
        return out;
    }
    const char *line = text;
    while (*line != '\0') {
        const char *end = std::strchr(line, '\n');
        const size_t len =
            end == nullptr ? std::strlen(line) : static_cast<size_t>(end - line);
        if (header_tag(line, len, "@HD") || header_tag(line, len, "@SQ") ||
            header_tag(line, len, "@RG")) {
            out.append(line, len);
            out.push_back('\n');
        }
        if (end == nullptr) {
            break;
        }
        line = end + 1;
    }
    return out;
}

long long file_size(const char *path) {
    struct stat info;
    if (::stat(path, &info) != 0) {
        return -1;
    }
    return static_cast<long long>(info.st_size);
}

struct SamCloser {
    void operator()(samFile *file) const {
        if (file != nullptr) {
            sam_close(file);
        }
    }
};

struct HeaderDeleter {
    void operator()(sam_hdr_t *header) const {
        if (header != nullptr) {
            sam_hdr_destroy(header);
        }
    }
};

struct BamDeleter {
    void operator()(bam1_t *record) const {
        if (record != nullptr) {
            bam_destroy1(record);
        }
    }
};

using SamPtr = std::unique_ptr<samFile, SamCloser>;
using HeaderPtr = std::unique_ptr<sam_hdr_t, HeaderDeleter>;
using BamPtr = std::unique_ptr<bam1_t, BamDeleter>;

bool set_aux_rg(bam1_t *dst, const bam1_t *src) {
    const uint8_t *rg = bam_aux_get(src, "RG");
    if (rg == nullptr) {
        return true;
    }
    const char *value = bam_aux2Z(rg);
    if (value == nullptr) {
        std::fprintf(stderr, "bam-slim: read has a malformed RG tag\n");
        return false;
    }
    return bam_aux_append(
               dst, "RG", 'Z', static_cast<uint32_t>(std::strlen(value) + 1),
               reinterpret_cast<const uint8_t *>(value)
           ) >= 0;
}

int run(const Options &options, const char *in_path, const char *out_path) {
    SamPtr input(sam_open(in_path, "r"));
    if (!input) {
        std::fprintf(stderr, "bam-slim: cannot open '%s'\n", in_path);
        return 1;
    }
    if (options.threads > 1 && hts_set_threads(input.get(), options.threads) != 0) {
        std::fprintf(stderr, "bam-slim: cannot set input threads\n");
        return 1;
    }

    HeaderPtr header(sam_hdr_read(input.get()));
    if (!header) {
        std::fprintf(stderr, "bam-slim: cannot read the header of '%s'\n", in_path);
        return 1;
    }

    const std::string header_text = slim_header(header.get());
    if (header_text.empty()) {
        std::fprintf(stderr, "bam-slim: '%s' has no usable header lines\n", in_path);
        return 1;
    }
    HeaderPtr out_header(
        sam_hdr_parse(header_text.size(), header_text.c_str())
    );
    if (!out_header) {
        std::fprintf(stderr, "bam-slim: cannot rebuild the output header\n");
        return 1;
    }

    SamPtr output(sam_open(out_path, "wb"));
    if (!output) {
        std::fprintf(stderr, "bam-slim: cannot create '%s'\n", out_path);
        return 1;
    }
    if (options.threads > 1 && hts_set_threads(output.get(), options.threads) != 0) {
        std::fprintf(stderr, "bam-slim: cannot set output threads\n");
        return 1;
    }
    if (sam_hdr_write(output.get(), out_header.get()) < 0) {
        std::fprintf(stderr, "bam-slim: cannot write the output header\n");
        return 1;
    }

    BamPtr record(bam_init1());
    BamPtr slim(bam_init1());
    if (!record || !slim) {
        std::fprintf(stderr, "bam-slim: cannot allocate a BAM record\n");
        return 1;
    }

    long long reads = 0;
    while (true) {
        const int read_status = sam_read1(input.get(), header.get(), record.get());
        if (read_status < 0) {
            if (read_status < -1) {
                std::fprintf(stderr, "bam-slim: error reading '%s'\n", in_path);
                return 1;
            }
            break;
        }

        const bam1_t *src = record.get();
        const char *name = options.keep_name ? bam_get_qname(src) : nullptr;
        const size_t l_qname =
            options.keep_name ? std::strlen(name) : 0;

        std::string sequence;
        std::string qualities;
        const char *seq_ptr = nullptr;
        const char *qual_ptr = nullptr;
        size_t l_seq = 0;
        if (!options.drop_seq && src->core.l_qseq > 0) {
            l_seq = static_cast<size_t>(src->core.l_qseq);
            sequence.resize(l_seq);
            const uint8_t *packed = bam_get_seq(src);
            for (size_t i = 0; i < l_seq; ++i) {
                sequence[i] = seq_nt16_str[bam_seqi(packed, i)];
            }
            seq_ptr = sequence.data();

            if (options.keep_qual) {
                const uint8_t *raw = bam_get_qual(src);
                // htslib uses 0xff for a missing quality string.
                if (raw[0] != 0xff) {
                    qualities.resize(l_seq);
                    for (size_t i = 0; i < l_seq; ++i) {
                        qualities[i] = static_cast<char>(raw[i] + 33);
                    }
                    qual_ptr = qualities.data();
                }
            }
        }

        const uint32_t *cigar =
            src->core.n_cigar > 0 ? bam_get_cigar(src) : nullptr;
        if (bam_set1(
                slim.get(), l_qname, name, src->core.flag, src->core.tid,
                src->core.pos, src->core.qual, static_cast<size_t>(src->core.n_cigar),
                cigar, src->core.mtid, src->core.mpos, src->core.isize, l_seq,
                seq_ptr, qual_ptr, 0
            ) < 0) {
            std::fprintf(stderr, "bam-slim: cannot rebuild a BAM record\n");
            return 1;
        }
        if (!options.drop_rg && !set_aux_rg(slim.get(), src)) {
            return 1;
        }
        if (sam_write1(output.get(), out_header.get(), slim.get()) < 0) {
            std::fprintf(stderr, "bam-slim: cannot write a BAM record\n");
            return 1;
        }
        ++reads;
    }

    if (sam_close(output.release()) < 0) {
        std::fprintf(stderr, "bam-slim: cannot finish '%s'\n", out_path);
        return 1;
    }
    if (options.index &&
        sam_index_build3(out_path, nullptr, 14, options.threads) < 0) {
        std::fprintf(stderr, "bam-slim: cannot index '%s'\n", out_path);
        return 1;
    }

    const long long in_bytes = file_size(in_path);
    const long long out_bytes = file_size(out_path);
    if (in_bytes > 0 && out_bytes >= 0) {
        std::fprintf(
            stderr,
            "bam-slim: %lld reads; %lld -> %lld bytes (%.1f%% smaller)\n",
            reads, in_bytes, out_bytes,
            100.0 * static_cast<double>(in_bytes - out_bytes) /
                static_cast<double>(in_bytes)
        );
    } else {
        std::fprintf(stderr, "bam-slim: %lld reads\n", reads);
    }
    return 0;
}

}  // namespace

int main(int argc, char **argv) {
    Options options;
    static const struct option long_options[] = {
        {"drop-seq", no_argument, nullptr, 's'},
        {"keep-qual", no_argument, nullptr, 'q'},
        {"keep-name", no_argument, nullptr, 'n'},
        {"drop-rg", no_argument, nullptr, 'r'},
        {"threads", required_argument, nullptr, '@'},
        {"index", no_argument, nullptr, 'i'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0},
    };

    int option = 0;
    while ((option = getopt_long(
                argc, argv, "h@:", long_options, nullptr
            )) != -1) {
        switch (option) {
        case 's':
            options.drop_seq = true;
            break;
        case 'q':
            options.keep_qual = true;
            break;
        case 'n':
            options.keep_name = true;
            break;
        case 'r':
            options.drop_rg = true;
            break;
        case '@':
            options.threads = std::atoi(optarg);
            if (options.threads < 1) {
                std::fprintf(stderr, "bam-slim: --threads must be at least 1\n");
                return 2;
            }
            break;
        case 'i':
            options.index = true;
            break;
        case 'h':
            usage(stdout);
            return 0;
        default:
            usage(stderr);
            return 2;
        }
    }

    if (optind + 2 != argc) {
        usage(stderr);
        return 2;
    }
    if (options.drop_seq && options.keep_qual) {
        std::fprintf(stderr, "bam-slim: --keep-qual requires SEQ; drop --drop-seq\n");
        return 2;
    }

    return run(options, argv[optind], argv[optind + 1]);
}
