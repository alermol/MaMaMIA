# Bundled Inter subsets

Both files are unmodified subsets from the Google Fonts distribution of the
Inter project, licensed under the SIL Open Font License 1.1 (`OFL.txt`).

| File | Inter version | Upstream URL | Size | SHA-256 |
|---|---|---|---|---|
| `inter-latin.woff2` | v20 (Google Fonts) | `https://fonts.gstatic.com/s/inter/v20/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa1ZL7.woff2` | 48,256 bytes | `3100e775e8616cd2611beecfa23a4263d7037586789b43f035236a2e6fbd4c62` |
| `inter-cyrillic.woff2` | v20 (Google Fonts) | `https://fonts.gstatic.com/s/inter/v20/UcC73FwrK3iLTeHuS_nVMrMxCp50SjIa0ZL7SUc.woff2` | 18,748 bytes | `71d5ee93cc1e9f1d520a3a8b66456de18c7879d8df09d57fcd2eaff75fef0075` |

The Latin URL above is the value published by the Google Fonts CSS API for the
current `inter-latin` file; the file itself predates this record, so its
checksum is recorded here. The Cyrillic subset was downloaded on 2026-11-12
from the same Inter v20 release so the two subsets stay metric-compatible.
The `@font-face` blocks in `../mamamia-font.css` use the unicode ranges
published alongside those URLs.
