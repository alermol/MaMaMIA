# Changelog

All notable changes to gcrank are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/). The on-disk format carries its own
`format_version`, which is versioned independently of the package.

## [0.2.0] - 2026-10-25

### Removed

- Unused internal helpers: `detail::println()`, `CRC32::reset()`,
  `ThreadPool::size()`, `PlainInputStream(FILE*, std::string)`, and the
  `close_enough()` test helper.
- Includes that had no direct users.

No public API or on-disk format changes.

## [0.1.0] - 2026-10-16

First release.
