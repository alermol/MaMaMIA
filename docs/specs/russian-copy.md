# MaMaMIA — Russian UI Copy (EN → RU)

- **Status:** final, ready-to-use copy for the P0–P2 i18n dictionary. Backs `russian-glossary.md`.
- **Key convention:** lowercase dot-path leaf keys grouped by UI region — `header.*`, `tabs.*`, `sidebar.*`, `params.*`, `results.*`, `errors.*`, `tooltip.*`, `karyotype.*`, `export.*`. Keys are identical in both languages.
- **Placeholders:** `{name}`/`%s`/`%d` are preserved exactly; `%s`/`%d` keep the same order as the English `sprintf` call sites.
- **Length flag:** a value marked **>1.6×** is an overflow risk for the fixed header/narrow sidebar (see `russian-localization.md` §6).
- **Source** column = `file:line` where the English literal lives. `p` marks strings built with `paste(...)` (adjacent literals).
- Rules applied: infinitive buttons, sentence case, ё used, Latin acronyms kept (see glossary §0).

---

## header.*

| key | EN | RU | source | flag |
|---|---|---|---|---|
| header.subtitle | Max Mapping for Meticulous Introgression Analysis | Max Mapping for Meticulous Introgression Analysis *(keep English — acronym expansion, glossary §2.1)* | R/app-bam.R:3528-3531 | — |
| header.note | Breadth-aware introgression mapping | Картирование интрогрессий по широте покрытия | R/app-bam.R:3537 | >1.6× (nowrap+uppercase; needs CSS `white-space: normal`) |
| header.language | (no EN literal — switcher uses EN/RU buttons; static aria-label) | Язык | — | — |

## tabs.*

| key | EN | RU | source |
|---|---|---|---|
| tabs.karyotype | Karyotype | Кариотип | R/app-bam.R:2946, R/dashboard.R:1041 |
| tabs.calls | Calls | Находки | R/app-bam.R:2950, R/dashboard.R:1045 |
| tabs.export | Export | Экспорт | R/app-bam.R:2970 |

## sidebar.*

| key | EN | RU | source | flag |
|---|---|---|---|---|
| sidebar.input | Input | Входные данные | R/app-bam.R:3551 | >1.6× (uppercase label) |
| sidebar.example_data | Example data | Пример данных | R/app-bam.R:3555 | — |
| sidebar.example_data_title | Load the bundled example dataset | Загрузить встроенный пример данных | R/app-bam.R:3558 | — |
| sidebar.clear_example_data | Clear example data | Очистить пример данных | R/app-bam.R:3746 | — |
| sidebar.analysis | Analysis | Анализ | R/app-bam.R:2724 | — |
| sidebar.assembly | Assembly | Сборка | R/app-bam.R:2728 | — |
| sidebar.advanced | Advanced | Дополнительно | R/app-bam.R:2732 | — |
| sidebar.bam_file | BAM file | BAM-файл | R/app-bam.R:2236,2259 | — |
| sidebar.max_size | Max %s | Макс. %s | R/app-bam.R:2240 | — |
| sidebar.about_upload_limit | About the BAM upload limit | О лимите загрузки BAM | R/app-bam.R:2245-2246 | — |
| sidebar.drag_browse | Drag & drop or browse | Перетащите файл или выберите | R/app-bam.R:2253 | — |
| sidebar.choose_bam | Choose BAM file | Выбрать BAM-файл | R/app-bam.R:2263 | — |
| sidebar.no_bam | No BAM selected | BAM не выбран | R/app-bam.R:2264 | — |
| sidebar.time_to_finish | Time to finish | Время до завершения | R/app-bam.R:2270 | — |
| sidebar.recipient_species | Recipient species | Вид-реципиент | R/app-bam.R:2699,2682 | — |
| sidebar.donor_species | Donor species | Вид-донор | R/app-bam.R:2704,2685 | — |
| sidebar.swap_species | Swap donor and recipient | Поменять донора и реципиента | R/app-bam.R:2709 | — |
| sidebar.no_species_registry | This registry has no species section, so the assembly is detected from the BAM header. | В реестре нет раздела видов, поэтому сборка определяется из заголовка BAM. | R/app-bam.R:2647-2650 p | — |
| sidebar.fixed_example_species | Fixed for the bundled example. Upload a BAM to change the assembly. | Зафиксировано для встроенного примера. Загрузите BAM, чтобы изменить сборку. | R/app-bam.R:2690-2691 p | — |
| sidebar.fixed_example_params | Analysis parameters are fixed for the bundled example. Upload a BAM to edit every profile parameter. | Параметры анализа зафиксированы для встроенного примера. Загрузите BAM, чтобы изменить параметры. | R/app-bam.R:2796-2797 p | — |
| sidebar.mode_example | Example data: events on eight chromosomes. | Пример данных: события на восьми хромосомах. | R/app-bam.R:3686 | — |
| sidebar.mode_intro | One BAM per run. The BAM header is matched to the configured assembly. | Один BAM на запуск. Заголовок BAM сопоставляется со сборкой. | R/app-bam.R:3693-3696 p | — |
| sidebar.mode_ready | BAM ready. The assembly was matched from the header. | BAM готов. Сборка определена по заголовку. | R/app-bam.R:3704 | — |

## sidebar.recovery.* (preflight failure panel)

| key | EN | RU | source |
|---|---|---|---|
| sidebar.recovery.title_too_large | BAM is too large | BAM слишком большой | R/app-bam.R:2384 |
| sidebar.recovery.title_assembly | No matching assembly | Сборка не найдена | R/app-bam.R:2386 |
| sidebar.recovery.title_attention | Setup needs attention | Требуется настройка | R/app-bam.R:2388 |
| sidebar.recovery.body_too_large | This BAM is larger than this server's upload limit, so it was not analysed. | Этот BAM превышает лимит загрузки сервера, поэтому не был проанализирован. | R/app-bam.R:2310-2311 p |
| sidebar.recovery.body_missing_contigs | This BAM does not contain the contigs required by the selected assembly. | Этот BAM не содержит контигов, требуемых выбранной сборкой. | R/app-bam.R:2319-2320 p |
| sidebar.recovery.note_missing_targets | The calculation log lists the missing targets. | Недостающие цели перечислены в журнале расчёта. | R/app-bam.R:2322 |
| sidebar.recovery.body_no_assembly | This BAM's header does not match any assembly in the registry. | Заголовок BAM не соответствует ни одной сборке в реестре. | R/app-bam.R:2330-2331 p |
| sidebar.recovery.note_choose_pair | Choose a matching species pair or add the assembly to the registry. | Выберите подходящую пару видов или добавьте сборку в реестр. | R/app-bam.R:2334-2335 p |
| sidebar.recovery.body_unknown_assembly | The BAM header names an assembly that is not in the registry. | Заголовок BAM ссылается на сборку, отсутствующую в реестре. | R/app-bam.R:2341-2342 p |
| sidebar.recovery.note_add_assembly | Add the assembly to the registry or choose a matching species pair. | Добавьте сборку в реестр или выберите подходящую пару видов. | R/app-bam.R:2345-2346 p |
| sidebar.recovery.body_generic | This BAM does not match a configured assembly. | Этот BAM не соответствует ни одной сборке. | R/app-bam.R:2351 |
| sidebar.recovery.note_check_header | Check the BAM header and the assembly setup. | Проверьте заголовок BAM и настройки сборки. | R/app-bam.R:2352 |
| sidebar.recovery.size_line | This file is %s; the limit is %s. | Размер файла — %s; лимит — %s. | R/app-bam.R:2399 |
| sidebar.recovery.next_step | Next step: | Следующий шаг: | R/app-bam.R:2430 |
| sidebar.recovery.slim_hint | Slim it locally with | Уменьшите его локально с помощью | R/app-bam.R:2412 |
| sidebar.recovery.slim_tail | , then upload the smaller BAM. | , затем загрузите уменьшенный BAM. | R/app-bam.R:2419 |
| sidebar.recovery.open_bamhusk | Open bamhusk | Открыть bamhusk | R/app-bam.R:2441,2629 |
| sidebar.recovery.choose_assembly | Choose assembly | Выбрать сборку | R/app-bam.R:2446 |
| sidebar.recovery.review_setup | Review setup | Проверить настройки | R/app-bam.R:2446 |

## sidebar.upload_help.* (upload-limit help modal)

| key | EN | RU | source |
|---|---|---|---|
| sidebar.upload_help.modal_title | BAM upload limit | Лимит загрузки BAM | R/app-bam.R:2618 |
| sidebar.upload_help.lead | This server accepts BAM files up to | Сервер принимает BAM-файлы размером до | R/app-bam.R:2526 |
| sidebar.upload_help.lead_tail | . The size is checked in your browser before upload. | . Размер проверяется в браузере до загрузки. | R/app-bam.R:2528 |
| sidebar.upload_help.limit | Limit | Лимит | R/app-bam.R:2533 |
| sidebar.upload_help.check | Check | Проверка | R/app-bam.R:2537 |
| sidebar.upload_help.before_upload | Before upload | До загрузки | R/app-bam.R:2538 |
| sidebar.upload_help.typical_saving | Typical saving | Типичная экономия | R/app-bam.R:2541 |
| sidebar.upload_help.slim_callout | If the BAM is larger, slim it locally with bamhusk. The command writes a new file and leaves the original untouched. | Если BAM больше, уменьшите его локально с помощью bamhusk. Команда создаёт новый файл и не изменяет исходный. | R/app-bam.R:2548-2550 p |
| sidebar.upload_help.slim_section | Slim it locally | Уменьшение локально | R/app-bam.R:2554 |
| sidebar.upload_help.get_bamhusk | Get bamhusk | Установить bamhusk | R/app-bam.R:2558 |
| sidebar.upload_help.run | Run | Запустить | R/app-bam.R:2569 |
| sidebar.upload_help.upload | Upload | Загрузить | R/app-bam.R:2573 |
| sidebar.upload_help.what_changes | What changes | Что меняется | R/app-bam.R:2585 |
| sidebar.upload_help.kept | Kept | Сохраняется | R/app-bam.R:2588 |
| sidebar.upload_help.removed | Removed | Удаляется | R/app-bam.R:2595 |
| sidebar.upload_help.optional | Optional | Необязательно | R/app-bam.R:2602 |
| sidebar.upload_help.command_aria | Command to slim the BAM | Команда для уменьшения BAM | R/app-bam.R:2580 |

## params.*

| key | EN | RU | source | flag |
|---|---|---|---|---|
| params.analysis_profile | Analysis profile | Профиль анализа | R/app-bam.R:4274,4246 | — |
| params.profile_custom | Custom analysis settings. | Пользовательские настройки анализа. | R/app-bam.R:4268 | — |
| params.profile_custom_based | Custom settings based on | Пользовательские настройки на основе | R/app-bam.R:4262 | — |
| params.profile_custom_tail | . Highlighted Advanced values differ from the profile defaults. | . Выделенные значения в «Дополнительно» отличаются от значений профиля по умолчанию. | R/app-bam.R:4264-4265 p | — |
| params.reset_all | Reset all to defaults | Сбросить всё к значениям по умолчанию | R/app-bam.R:2825 | — |
| params.reset_title | Restore every profile parameter and control sample to the selected profile's defaults | Вернуть все параметры профиля и контрольные образцы к значениям выбранного профиля по умолчанию | R/app-bam.R:2829-2830 p | — |
| params.heading_resolution | Resolution | Разрешение | R/app-bam.R:2834 | — |
| params.heading_scoring | Scoring | Оценка | R/app-bam.R:2836 | — |
| params.heading_masking | Masking | Маскирование | R/app-bam.R:2877 | — |
| params.heading_calibration | Calibration | Калибровка | R/app-bam.R:2885 | — |
| params.masking_note | Cohort recurrence cutoff is derived automatically from the fitted cohort mixture. Single-sample runs are unmasked. | Порог рекуррентности когорты выводится автоматически из подобранной смеси. Запуски с одним образцом не маскируются. | R/app-bam.R:2881-2883 p | — |
| params.legend_bold | Bold | Жирным | R/app-bam.R:2919 | — |
| params.legend_values |  values are the selected profile's defaults; |  выделены значения профиля по умолчанию; | R/app-bam.R:2920 | — |
| params.legend_highlighted | highlighted | выделенные | R/app-bam.R:2921 | — |
| params.legend_custom |  values were changed and make the run Custom. |  значения были изменены и делают запуск «Пользовательским». | R/app-bam.R:2922 | — |
| params.label_resolution | Resolution | Разрешение | R/app-bam.R:4369 | — |
| params.label_metric | Metric | Метрика | R/app-bam.R:4375 | — |
| params.label_score | Score | Оценка | R/app-bam.R:4376 | — |
| params.label_donor_aggregation | Donor aggregation | Агрегация доноров | R/app-bam.R:4377 | — |
| params.label_mapq | Minimum MAPQ | Мин. MAPQ | R/app-bam.R:4378 | — |
| params.label_calibration_fdr | Calibration FDR | FDR калибровки | R/app-bam.R:4379 | — |
| params.label_calibration_min_mb | Calibration minimum length (Mb) | Мин. длина калибровки (Мб) | R/app-bam.R:4380 | — |
| params.label_calibration_merge_gap | Calibration merge gap (windows) | Зазор слияния калибровки (окна) | R/app-bam.R:4381 | — |
| params.floor_label | Floor %s Mb | Порог %s Мб | R/app-bam.R:2760 | — |
| params.metric_breadth | Breadth | Широта покрытия | R/app-bam.R:2842 | — |
| params.metric_count | Count | Число прочтений | R/app-bam.R:2843 | — |
| params.metric_combined | Combined | Комбинированный | R/app-bam.R:2844 | — |
| params.score_paired | Paired log-ratio | Парный лог-коэффициент | R/app-bam.R:2854 | — |
| params.score_donor_only | Donor only | Только донор | R/app-bam.R:2855 | — |
| params.agg_max | Max | Макс. | R/app-bam.R:2864 | — |
| params.agg_mean | Mean | Среднее | R/app-bam.R:2864 | — |
| params.agg_sum | Sum | Сумма | R/app-bam.R:2864 | — |
| params.control_samples | Calibration controls | Контрольные образцы калибровки | R/app-bam.R:4581 | — |
| params.control_placeholder | All samples (default) | Все образцы (по умолчанию) | R/app-bam.R:4595 | — |
| params.control_note | Used only to calibrate call thresholds. Leave empty to pool all samples. | Используются только для калибровки порогов. Оставьте пустым, чтобы объединить все образцы. | R/app-bam.R:4601-4603 p | — |
| params.control_multi_note | Control samples apply only to multi-sample cohorts. | Контрольные образцы применимы только к когортам из нескольких образцов. | R/app-bam.R:4572 | — |
| params.cohort_warning_title | Cohort recurrence unavailable. | Рекуррентность когорты недоступна. | R/app-bam.R:2780 | — |
| params.cohort_warning_body | This run has one sample. Cohort recurrence masking requires at least two samples, so no windows will be masked. | В этом запуске один образец. Маскирование по рекуррентности когорты требует минимум два образца, поэтому окна не будут замаскированы. | R/app-bam.R:2782-2783 p | — |
| params.example_profile_note | Fixed for the bundled example. | Зафиксировано для встроенного примера. | R/app-bam.R:4251 | — |
| params.resolution_aria | Exact resolution in windows per chromosome | Точное разрешение в окнах на хромосому | R/app-bam.R:443-444 | — |
| params.resolution_locked_title | Resolution is fixed for the bundled example | Разрешение зафиксировано для встроенного примера | R/app-bam.R:482-483 | — |
| params.resolution_fixed_note | Fixed at %s windows per chromosome. Upload a BAM to change it. | Зафиксировано: %s окон на хромосому. Загрузите BAM, чтобы изменить. | R/app-bam.R:489-492 p | — |

## params.help.* (parameter tooltips)

| key | EN | RU | source |
|---|---|---|---|
| params.help.resolution | Windows per chromosome. More windows give finer breakpoints but make the run slower. | Окон на хромосому. Больше окон даёт более точные границы, но замедляет расчёт. | R/app-bam.R:1670-1673 p |
| params.help.metric | Window metric used by the score. Breadth uses covered bases, count uses read depth, and combined uses both. | Метрика окна для оценки. «Широта покрытия» считает покрытые основания, «Число прочтений» — глубину, «Комбинированный» — обе. | R/app-bam.R:1674-1677 p |
| params.help.score | Paired log-ratio compares donor and recipient coverage. Donor only ignores the recipient and tolerates copy-number changes in the recipient. | Парный лог-коэффициент сравнивает покрытие донора и реципиента. «Только донор» игнорирует реципиента и допускает изменения числа копий у реципиента. | R/app-bam.R:1678-1682 p |
| params.help.donor_aggregation | How several donor components are combined in one window: max, mean, or sum. | Как несколько компонентов донора объединяются в одном окне: максимум, среднее или сумма. | R/app-bam.R:1683-1686 p |
| params.help.mapq | Minimum mapping quality for a read to count. Higher values exclude more ambiguous alignments. | Минимальное качество картирования для учёта прочтения. Более высокие значения исключают больше неоднозначных выравниваний. | R/app-bam.R:1687-1690 p |
| params.help.calibration_fdr | Target false discovery rate for segment thresholds. Lower values make calls more conservative. | Целевой уровень ложных находок (FDR) для порогов сегментов. Меньшие значения делают находки консервативнее. | R/app-bam.R:1691-1694 p |
| params.help.calibration_min_mb | Minimum call width in megabases. Shorter segments are excluded. | Минимальная ширина находки в мегабазах. Более короткие сегменты исключаются. | R/app-bam.R:1695-1698 p |
| params.help.calibration_merge_gap | Largest gap, in windows, between positive segments that may be merged into one call. | Максимальный зазор в окнах между положительными сегментами, которые могут быть объединены в одну находку. | R/app-bam.R:1699-1702 p |
| params.help.calibration_floor | Minimum mean score for a call in the %s length bin. Raise it to reject weak calls in this size range. | Минимальная средняя оценка для находки в интервале длины %s. Повысьте, чтобы отсеять слабые находки в этом диапазоне. | R/app-bam.R:1661-1665 p |
| params.help.control_samples | Samples used only to calibrate thresholds. Leave empty to pool every sample in the cohort. | Образцы, используемые только для калибровки порогов. Оставьте пустым, чтобы объединить все образцы когорты. | R/app-bam.R:1703-1706 p |

## run.* / progress.* (run button, loader, stage labels)

| key | EN | RU | source |
|---|---|---|---|
| run.run_analysis | Run analysis | Запустить анализ | R/app-bam.R:4812 |
| run.preparing | Preparing analysis | Подготовка анализа | R/app-bam.R:5066, R/dashboard.R:852 |
| run.analysis_eyebrow | Analysis | Анализ | R/dashboard.R:850 |
| run.processing_time | Processing time | Время обработки | R/dashboard.R:858 |
| progress.preflight | Checking BAM | Проверка BAM | R/app-bam.R:1937 |
| progress.index | Indexing BAM | Индексирование BAM | R/app-bam.R:1938 |
| progress.windows | Building windows | Построение окон | R/app-bam.R:1939 |
| progress.gc | Reading GC content | Чтение GC-состава | R/app-bam.R:1940 |
| progress.metrics | Counting reads | Подсчёт прочтений | R/app-bam.R:1941 |
| progress.score | Scoring windows | Оценка окон | R/app-bam.R:1942 |
| progress.mask | Masking windows | Маскирование окон | R/app-bam.R:1943 |
| progress.segmentation | Segmenting profiles | Сегментация профилей | R/app-bam.R:1944 |
| progress.calibration | Calibrating calls | Калибровка находок | R/app-bam.R:1945 |

## footer.*

| key | EN | RU | source |
|---|---|---|---|
| footer.show_log | Show log | Показать журнал | R/app-bam.R:3594,3250 |
| footer.hide_log | Hide log | Скрыть журнал | R/app-bam.R:3250 |
| footer.expand_log | Expand the calculation log | Развернуть журнал расчёта | R/app-bam.R:3593 |
| footer.collapse_log | Collapse the calculation log | Свернуть журнал расчёта | R/app-bam.R:3255 |
| footer.preflight_ok | Preflight OK | Проверка пройдена | R/app-bam.R:4625 |
| footer.waiting | Waiting | Ожидание | R/app-bam.R:4633 |
| footer.upload_hint | Upload a BAM or use the example data. | Загрузите BAM или используйте пример данных. | R/app-bam.R:4635 |
| footer.bam_too_large | BAM too large | BAM слишком большой | R/app-bam.R:4641 |
| footer.needs_assembly | Needs assembly | Нужна сборка | R/app-bam.R:4643 |
| footer.preflight_failed | Preflight failed | Проверка не пройдена | R/app-bam.R:4645 |
| footer.see_log | See the calculation log. | См. журнал расчёта. | R/app-bam.R:4653 |
| footer.slim_smaller | Slim it locally with bamhusk, then upload the smaller BAM. | Уменьшите локально с помощью bamhusk и загрузите уменьшенный BAM. | R/app-bam.R:4655-4658 p |
| footer.choose_assembly | Choose a matching assembly in the sidebar. | Выберите подходящую сборку на боковой панели. | R/app-bam.R:4660 |
| footer.example_prefix | Example data \| | Пример данных \| | R/app-bam.R:4669 |
| footer.bam_file_count | %d BAM file%s | BAM-файлов: %d *(plural rework — see §plural)* | R/app-bam.R:4628-4630 |

## errors.* (notifications, modals, validation messages)

| key | EN | RU | source |
|---|---|---|---|
| errors.modal_too_large_title | BAM is too large to upload | BAM слишком большой для загрузки | R/app-bam.R:2470 |
| errors.file_label | File | Файл | R/app-bam.R:202-220 |
| errors.file_label_plural | Files | Файлы | R/app-bam.R:202-220 |
| errors.file_more | (+%d more) | (ещё %d) | R/app-bam.R:202-220 |
| errors.modal_choose_another | Choose another file | Выбрать другой файл | R/app-bam.R:2504 |
| errors.modal_close | Close | Закрыть | R/app-bam.R:2501,2631; R/dashboard.R:1217 |
| errors.fix_preflight | Fix the preflight problems before running the analysis. | Устраните проблемы проверки перед запуском анализа. | R/app-bam.R:4923 |
| errors.run_refused_log | Run refused: preflight is not ready. Fix the problems reported in this log first. | Запуск отклонён: проверка не готова. Сначала устраните проблемы из журнала. | R/app-bam.R:4918-4920 p |
| errors.bam_count | The dashboard analyses one BAM per run; %d files were supplied. A multi-sample BAM is one file. | Приложение анализирует один BAM за запуск; передано файлов: %d. Многообразцовый BAM — это один файл. *(plural rework)* | R/app-bam.R:585-588 p |
| errors.resolution_range | Resolution must be between %s and %s in steps of %s | Разрешение должно быть от %s до %s с шагом %s | R/app-bam.R:531 |
| errors.choose_resolution | Choose an analysis resolution | Выберите разрешение анализа | R/app-bam.R:523 |
| errors.profile_resolution | Analysis profile '%s' uses %s windows per chromosome; this configuration allows %s-%s in steps of %s | Профиль анализа '%s' использует %s окон на хромосому; эта конфигурация допускает %s–%s с шагом %s | R/app-bam.R:1783-1786 |
| errors.whole_number | %s must be a whole number between %s and %s | %s должно быть целым числом от %s до %s | R/app-bam.R:547,1506 |
| errors.number | %s must be a number between %s and %s | %s должно быть числом от %s до %s | R/app-bam.R:1743 |
| errors.at_least_sample | At least one sample is required | Требуется хотя бы один образец | R/app-bam.R:1768 |
| errors.calibration_floor | Each calibration floor must be a non-negative number | Каждый порог калибровки должен быть неотрицательным числом | R/app-bam.R:1855 |
| errors.unknown_control | Unknown control sample(s): %s | Неизвестный контрольный образец (образцы): %s | R/app-bam.R:1884 |
| errors.choose_pair | Choose a recipient species and a donor species | Выберите вид-реципиент и вид-донор | R/app-bam.R:4024 |
| errors.unknown_species | Unknown species selection ('%s' recipient, '%s' donor) | Неизвестный выбор видов ('%s' — реципиент, '%s' — донор) | R/app-bam.R:4035 |
| errors.short_contigs | Contig%s %s %s shorter than the requested %s windows; choose a lower resolution. | *(needs two-variant template — see §plural)* | R/app-bam.R:1307-1312 |

## log.* (footer calculation log — high-visibility entries; full tail below in notes)

| key | EN | RU | source |
|---|---|---|---|
| log.started | MaMaMIA session started. | Сеанс MaMaMIA запущен. | R/app-bam.R:916 |
| log.config_validated | Configuration validated. | Конфигурация проверена. | R/app-bam.R:917 |
| log.header_parsed | BAM header parsed. | Заголовок BAM разобран. | R/app-bam.R:918 |
| log.example_selected | Example data selected. | Выбран пример данных. | R/app-bam.R:3749 |
| log.example_cleared | Example data cleared. Upload a BAM to start a new run. | Пример данных очищен. Загрузите BAM для нового запуска. | R/app-bam.R:3736 |
| log.example_note | Example dataset: %s recipient and %s donor; events on eight chromosomes. | Пример данных: %s — реципиент, %s — донор; события на восьми хромосомах. | R/example-data.R:277 |
| log.detected_species | Detected species from the BAM header: %s. | Вид определён из заголовка BAM: %s. | R/app-bam.R:3961 |
| log.detected_pair | %s recipient, %s donor | %s — реципиент, %s — донор | R/app-bam.R:3925 |
| log.analysis_start | Starting analysis: %s windows/chromosome, metric '%s', score '%s', donor aggregation '%s'. | Запуск анализа: %s окон/хромосома, метрика '%s', оценка '%s', агрегация доноров '%s'. | R/app-bam.R:5011-5014 p |
| log.analysis_finished | Analysis finished: %d %s. | Анализ завершён — откалиброванных находок: %d. *(plural rework — drops 2nd %s)* | R/app-bam.R:5112 |
| log.coverage_warning | Coverage warning: %s | Предупреждение о покрытии: %s | R/app-bam.R:5123 |
| log.calibration_controls | Calibration controls: | Контрольные образцы калибровки: | R/app-bam.R:5019-5021 |
| log.suggested_fix | Suggested fix: %s | Рекомендуемое решение: %s | R/app-bam.R:845,4946,5089 |
| log.reset_params | Reset every Advanced parameter to the %s defaults. | Сброшены все параметры «Дополнительно» к значениям %s по умолчанию. | R/app-bam.R:4331 |
| log.reset_params_notify | Advanced parameters reset to %s defaults. | Параметры «Дополнительно» сброшены к значениям %s по умолчанию. | R/app-bam.R:4334 |
| log.restored_cached | Restored a cached analysis result for these inputs; no recalculation was needed. | Восстановлен кэшированный результат для этих входных данных; пересчёт не потребовался. | R/app-bam.R:4845-4846 p |
| log.cleared_result | Analysis inputs changed; the previous result was cleared. Run the analysis for the new inputs. | Входные данные изменились; предыдущий результат очищен. Запустите анализ для новых данных. | R/app-bam.R:4855-4857 p |
| log.using_cached | Using the cached analysis result; no recalculation was needed. | Используется кэшированный результат; пересчёт не потребовался. | R/app-bam.R:4996-4998 p |
| log.using_cached_notify | Using the cached analysis result for these inputs. | Используется кэшированный результат для этих входных данных. | R/app-bam.R:5445 |
| log.setup_panels | Setup panels opened after a preflight failure. | Панели настройки открыты после ошибки проверки. | R/app-bam.R:3769 |

## results.*

| key | EN | RU | source | flag |
|---|---|---|---|---|
| results.sample | Sample | Образец | R/app-bam.R:5136 | — |
| results.call_note | Brief view of passing calls. Hover a column heading for its description; use Additional columns to add diagnostics. Shift-click headings to sort by several columns. Turn off Passing only to include calls that failed calibration. The full table, BEDPE files and images are in the results bundle. | Краткий вид прошедших находок. Наведите курсор на заголовок столбца для описания; используйте «Дополнительные столбцы» для диагностики. Shift+клик по заголовкам сортирует по нескольким столбцам. Отключите «Только прошедшие», чтобы включить находки, не прошедшие калибровку. Полная таблица, файлы BEDPE и изображения — в архиве результатов. | R/app-bam.R:2959-2964 p | >1.6× |
| results.filter_class | Call class | Класс находки | R/dashboard.R:71 | — |
| results.filter_class_title | Class | Класс | R/dashboard.R:73 | — |
| results.filter_status | Class status | Статус класса | R/dashboard.R:81 | — |
| results.filter_status_title | Status | Статус | R/dashboard.R:82 | — |
| results.filter_copy | Copy state | Копийность | R/dashboard.R:89 | — |
| results.filter_copy_title | Copy | Копии | R/dashboard.R:90 | — |
| results.filter_confidence | Confidence | Достоверность | R/dashboard.R:99 | — |
| results.filter_passing | Passing only | Только прошедшие | R/dashboard.R:106 | — |
| results.picker_selected | {0} selected | {0} выбрано | R/dashboard.R:63 | — |
| results.min_width | Minimum call width (Mb) | Мин. ширина находки (Мб) | R/dashboard.R:116 | — |
| results.additional_columns | Additional columns | Дополнительные столбцы | R/dashboard.R:125 | — |
| results.add_columns | Add columns | Добавить столбцы | R/dashboard.R:133-134 | — |
| results.picker_added | {0} added | {0} добавлено | R/dashboard.R:136 | — |
| results.select_visible | Select visible | Выбрать видимые | R/dashboard.R:175 | — |
| results.clear | Clear | Очистить | R/dashboard.R:178 | — |
| results.no_description | No description available. | Описание недоступно. | R/export-bundle.R:141 | — |
| results.run_to_list | Run the analysis to list calls | Запустите анализ, чтобы увидеть список находок | R/dashboard.R:1350 | — |
| results.run_to_plot | Run the analysis to plot | Запустите анализ, чтобы построить график | R/dashboard.R:1334 | — |
| results.run_hint | Load a BAM or the bundled example data, then run the analysis to map introgressions. | Загрузите BAM или встроенный пример данных, затем запустите анализ для картирования интрогрессий. | R/dashboard.R:1335-1337 p | — |
| results.run_cohort | Run the analysis to see the cohort summary | Запустите анализ, чтобы увидеть сводку когорты | R/app-bam.R:5161 | — |
| results.run_before_download | Run the analysis before downloading | Запустите анализ перед скачиванием | R/dashboard.R:1469,1521 | — |
| results.select_before_download | Select at least one call before downloading | Выберите хотя бы одну находку перед скачиванием | R/dashboard.R:1558 | — |
| results.selection_status | call(s) selected for BEDPE export | Выбрано для BEDPE-экспорта: %d *(plural rework)* | R/dashboard.R:1192 | — |
| results.no_analysis | No analysis loaded yet | Анализ ещё не загружен | R/dashboard.R:1190 | — |
| results.unplaced_none | No unplaced calls | Нет неразмещённых находок | R/dashboard.R:1168 | — |
| results.unplaced_hidden | %d unplaced hidden | Скрыто неразмещённых: %d *(plural rework)* | R/dashboard.R:1171 | — |
| results.unplaced_shown | %d unplaced shown | Показано неразмещённых: %d *(plural rework)* | R/dashboard.R:1173 | — |
| results.deletions_shown | Deletions shown | Делеции показаны | R/dashboard.R:1306 | — |
| results.deletions_none | No deletions omitted | Делеции не опущены | R/dashboard.R:1311 | — |
| results.deletions_omitted | %d segment%s omitted | Опущено сегментов: %d *(plural rework)* | R/dashboard.R:1313-1317 | — |

## karyotype.* (plot controls + widget JS labels)

| key | EN | RU | source | flag |
|---|---|---|---|---|
| karyotype.plot | Plot | Вид | R/dashboard.R:883 | — |
| karyotype.heatmap | Heatmap | Тепловая карта | R/dashboard.R:884 | — |
| karyotype.calls | Calls | Находки | R/dashboard.R:884 | — |
| karyotype.hide_deletions | Hide deletions | Скрыть делеции | R/dashboard.R:891 | — |
| karyotype.hide_unplaced | Hide unplaced | Скрыть неразмещённые | R/dashboard.R:898 | — |
| karyotype.download_current_title | Download the karyotype exactly as filtered and displayed (SVG). | Скачать кариотип в точности как отфильтровано и показано (SVG). | R/dashboard.R:914-915 p | — |
| karyotype.current_view | Current view | Текущий вид | R/dashboard.R:918 | — |
| karyotype.download_full_title | Download every call for this sample, ignoring the active filters (SVG). | Скачать все находки для этого образца, игнорируя активные фильтры (SVG). | R/dashboard.R:926-927 p | — |
| karyotype.full_view | Full view | Полный вид | R/dashboard.R:930 | — |
| karyotype.tooltip_help | Tooltip help | Справка по подсказкам | R/dashboard.R:933 | — |
| karyotype.aria_label | Karyotype in assembly order | Кариотип в порядке сборки | www/htmlwidgets/mamamiaKaryotypeCard.js:751 | — |
| karyotype.empty | No chromosomes in this result. | В этом результате нет хромосом. | mamamiaKaryotypeCard.js:661 | — |
| karyotype.legend_substitution | Substitution | Замещение | R/karyotype-scene.R:86 | — |
| karyotype.legend_addition | Addition | Добавление | R/karyotype-scene.R:87 | — |
| karyotype.legend_deletion | Deletion | Делеция | R/karyotype-scene.R:88 | — |
| karyotype.legend_complex | Complex | Сложный | R/karyotype-scene.R:89 | — |
| karyotype.legend_other | Other | Прочее | R/karyotype-scene.R:90 | — |
| karyotype.legend_deleted_recipient | Deleted recipient | Реципиент с делецией | R/karyotype-scene.R:91 | — |
| karyotype.legend_added_donor | Added donor | Добавленный донор | R/karyotype-scene.R:92 | — |
| karyotype.legend_unplaced_donor | Unplaced donor | Неразмещённый донор | R/karyotype-scene.R:93 | — |
| karyotype.legend_score | Score | Оценка | mamamiaKaryotypeCard.js:110 | — |
| karyotype.caption_deleted | Deleted | Делетирован | R/karyotype-scene.R:373 | — |
| karyotype.caption_unplaced | Unplaced | Не размещён | R/karyotype-scene.R:375 | — |
| karyotype.caption_added | Added | Добавлен | R/karyotype-scene.R:377 | — |

## tooltip.* (no-call card + matched-call card, JS + R)

| key | EN | RU | source |
|---|---|---|---|
| tooltip.call | call | находка | R/call-tooltip.R:129,234; mamamiaKaryotypeCard.js:290,385 |
| tooltip.no_call | no call | нет находки | mamamiaKaryotypeCard.js:417 |
| tooltip.unplaced_donor_title | unplaced donor | неразмещённый донор | mamamiaKaryotypeCard.js:407 |
| tooltip.deleted_recipient_title | deleted recipient | реципиент с делецией | mamamiaKaryotypeCard.js:409 |
| tooltip.added_donor_title | added donor | добавленный донор | mamamiaKaryotypeCard.js:411 |
| tooltip.recipient_chromosome_title | recipient chromosome | хромосома реципиента | mamamiaKaryotypeCard.js:412 |
| tooltip.component | Component | Компонент | mamamiaKaryotypeCard.js:424 |
| tooltip.state | State | Состояние | mamamiaKaryotypeCard.js:427,432,435,437 |
| tooltip.state_unplaced | unplaced donor chromosome | неразмещённая донорная хромосома | mamamiaKaryotypeCard.js:428 |
| tooltip.state_deleted | deleted recipient (no coverage) | реципиент с делецией (нет покрытия) | mamamiaKaryotypeCard.js:432 |
| tooltip.state_added | added donor chromosome | добавленная донорная хромосома | mamamiaKaryotypeCard.js:435 |
| tooltip.state_recipient | recipient chromosome | хромосома реципиента | mamamiaKaryotypeCard.js:437 |
| tooltip.length | Length | Длина | mamamiaKaryotypeCard.js:440 |
| tooltip.chromosome | Chromosome | Хромосома | mamamiaKaryotypeCard.js:442 |
| tooltip.position | Position | Позиция | mamamiaKaryotypeCard.js:461,469 |
| tooltip.score | Score | Оценка | mamamiaKaryotypeCard.js:465,110 |
| tooltip.click_to_filter | Click to filter this group | Нажмите, чтобы отфильтровать эту группу | mamamiaKaryotypeCard.js:472-473 |
| tooltip.group | group | группа | mamamiaKaryotypeCard.js:400; R/call-tooltip.R:251 |
| tooltip.no_segment | No assigned donor segment at this position. | В этой позиции нет назначенного донорного сегмента. | mamamiaKaryotypeCard.js:447 |
| tooltip.placement | Placement | Размещение | R/call-tooltip.R:447 |
| tooltip.evidence | Evidence | Свидетельства | R/call-tooltip.R:456 |
| tooltip.geometry | Geometry | Геометрия | R/call-tooltip.R:465 |
| tooltip.relation | Relation | Связь | R/call-tooltip.R:326 |
| tooltip.assignment | Assignment | Размещение | R/call-tooltip.R:337 |
| tooltip.group_class | Group class | Класс группы | R/call-tooltip.R:343 |
| tooltip.alternative | Alternative | Альтернатива | R/call-tooltip.R:351 |
| tooltip.donor_width | Donor width | Ширина донора | R/call-tooltip.R:415 |
| tooltip.call_width | Call width | Ширина находки | R/call-tooltip.R:415 |
| tooltip.recipient_run | Recipient run | Делеция реципиента | R/call-tooltip.R:425 |
| tooltip.overlap | Overlap | Перекрытие | R/call-tooltip.R:432 |
| tooltip.width_ratio | Width ratio | Отношение ширин | R/call-tooltip.R:439 |
| tooltip.no_segment_assigned | No donor or recipient segment assigned. | Не назначен ни донорный, ни реципиентный сегмент. | R/call-tooltip.R:452 |
| tooltip.same_group | same group | та же группа | R/call-tooltip.R:324 |
| tooltip.cross_group | cross-group | межгрупповой | R/call-tooltip.R:322 |
| tooltip.not_placed | not placed | не размещено | R/call-tooltip.R:320 |
| tooltip.unassigned | unassigned | не назначено | R/call-tooltip.R:277 |
| tooltip.no_host | no host | нет реципиента | R/call-tooltip.R:275 |
| tooltip.unknown | unknown | неизвестно | R/call-tooltip.R:302 |
| tooltip.confident | confident | уверенный | R/call-tooltip.R:123 |
| tooltip.ambiguous | ambiguous | неоднозначный | R/call-tooltip.R:125 |
| tooltip.needs_review | needs_review | требует проверки | R/call-tooltip.R:127 |
| tooltip.chip_class_confidence | class confidence | достоверность класса | R/call-tooltip.R:385 |
| tooltip.chip_evidence | evidence | свидетельства | R/call-tooltip.R:390 |
| tooltip.chip_copy | copy | копии | R/call-tooltip.R:393 |
| tooltip.chip_total | total | суммарно | R/call-tooltip.R:398 |
| tooltip.chip_components_heterogeneous | components heterogeneous | компоненты неоднородны | R/call-tooltip.R:401 |
| tooltip.chip_low_mappability | low mappability | низкая картированность | R/call-tooltip.R:404 |
| tooltip.chip_not_passing | not passing | не проходит | R/call-tooltip.R:407 |
| tooltip.legacy_unavailable | Call details unavailable in this legacy payload. | Детали находки недоступны в этом устаревшем наборе данных. | mamamiaKaryotypeCard.js:392 |

## export.* (Export tab)

| key | EN | RU | source |
|---|---|---|---|
| export.results_bundle | Results bundle | Архив результатов | R/dashboard.R:211 |
| export.results_lead | Download full and confident call tables, all-call, passing and confident BEDPE files, and matching karyotype SVGs in one ZIP archive. The archive has one folder per sample plus run-wide summaries at the top level. | Скачайте полные и уверенные таблицы находок, BEDPE-файлы (все, прошедшие и уверенные) и соответствующие SVG-кариотипы одним ZIP-архивом. В архиве по одной папке на образец плюс сводки по запуску на верхнем уровне. | R/dashboard.R:215-219 p |
| export.download_bundle | Download results bundle (.zip) | Скачать архив результатов (.zip) | R/dashboard.R:224 |
| export.preparing_bundle | Preparing bundle... | Подготовка архива... | R/dashboard.R:226 |
| export.all_calls | All calls | Все находки | R/dashboard.R:231 |
| export.all_calls_lead | Download every call in the current run as one BEDPE file, without the rest of the bundle. | Скачать все находки текущего запуска одним BEDPE-файлом, без остального архива. | R/dashboard.R:235-236 p |
| export.download_all_bedpe | Download all BEDPE | Скачать весь BEDPE | R/dashboard.R:241 |
| export.selected_calls | Selected calls | Выбранные находки | R/dashboard.R:248 |
| export.download_selected_bedpe | Download selected BEDPE | Скачать выбранный BEDPE | R/dashboard.R:252 |
| export.preparing_bedpe | Preparing BEDPE... | Подготовка BEDPE... | R/dashboard.R:243,254 |
| export.preparing_download | Preparing download... | Подготовка скачивания... | R/app-bam.R:3405; R/dashboard.R:184 |
| export.download_failed | Download failed. Try again. | Скачивание не удалось. Попробуйте снова. | R/app-bam.R:3424 |
| export.tooltip_guide | Tooltip guide | Справка по подсказкам | R/dashboard.R:1213 |

---

## Plural-avoidance rewrite list (Phase 1: no plural machinery)

These templates join a bare number to a noun in English. The Russian value is rephrased to "Label: N" (number last) so no 1/2-4/5+ agreement is needed. **The `%s`/`%d` arguments change at these call sites** — the engineer must adapt the `sprintf` call, not just the template. Each is listed with the required code-side change.

| key | EN template | RU template | code change |
|---|---|---|---|
| errors.bam_count | The dashboard analyses one BAM per run; %d files were supplied. A multi-sample BAM is one file. | Приложение анализирует один BAM за запуск; передано файлов: %d. Многообразцовый BAM — это один файл. | keep `%d` (count) only |
| footer.bam_file_count | %d BAM file%s | BAM-файлов: %d | drop the trailing `%s` suffix argument |
| log.analysis_finished | Analysis finished: %d %s. | Анализ завершён — откалиброванных находок: %d. | pass `n_calls` only; drop `call_word` |
| results.selection_status | call(s) selected for BEDPE export | Выбрано для BEDPE-экспорта: %d | wrap `paste()` → `sprintf` with one `%d` |
| results.unplaced_hidden | %d unplaced hidden | Скрыто неразмещённых: %d | — |
| results.unplaced_shown | %d unplaced shown | Показано неразмещённых: %d | — |
| results.deletions_omitted | %d segment%s omitted | Опущено сегментов: %d | drop the `%s` suffix argument |
| errors.short_contigs | Contig%s %s %s shorter than the requested %s windows; choose a lower resolution. | *(two variants, no single template)* | use a 2-key split (singular/plural) or `mamamia_tr_plural`; RU: "Контиг %s короче запрошенных %s окон; выберите меньшее разрешение." / "Контиги %s короче запрошенных %s окон; выберите меньшее разрешение." |

The notification `run_summary` (`R/app-bam.R:5115`) mirrors `log.analysis_finished`; apply the same single-`%d` rephrasing.

## Key-group coverage confirmation (P0–P2)

- **Fully covered:** `header.*`, `tabs.*`, `sidebar.*` (incl. `sidebar.recovery.*`, `sidebar.upload_help.*`), `params.*` (incl. `params.help.*`), `run.*`, `progress.*`, `footer.*`, `errors.*` (notifications/modals/validation), `log.*` (high-visibility entries), `results.*` (filters, table, counters, empty states), `karyotype.*` (controls + legend + captions + widget JS labels), `tooltip.*` (no-call card + matched-call card), `export.*` (Export tab cards + buttons).
- **Partially covered / flagged:** the full "Tooltip guide" modal body (`explorer_tooltip_help_ui()`, `R/dashboard.R:260-820`, ~100 strings: "Quick interpretation", "Header", "Placement", "Evidence", "Status chips", "Geometry", "Position and footer", "No-call card", "Hide deletions", "Unplaced calls", "Caveats" sections) is **not** fully enumerated here — its section headings and row terms map 1:1 to the glossary (class/copy-state/status/placement values in `russian-glossary.md` §6) and the `tooltip.*` labels above; the engineer should extract its remaining prose as a `tooltip_guide.*` block in a follow-up commit, reusing the glossary values.
- **No literal in source (do not invent keys):** `header.language` has no English literal (the switcher is EN/RU buttons with a static bilingual aria-label); `params.help.calibration_floor` is generated from `Floor %s Mb` by `sub()` — localize `params.floor_label` instead.
- **Deferred (P3, out of scope — listed for parity):** `export-bundle*`, `cohort-export` (`"Calls that pass calibration filters."` R/cohort-export.R:368, `"Calls by class"` R/export-bundle-plots.R:385, `"Calls by class across the run"` R/export-bundle-cohort.R:101), demo-data content, Russian plural machinery, persistence.

## English-string correctness flags (report only — no silent fix)

1. **`unplaced donor`** (legend + no-call title) is a *donor contig state*, not a call — do not translate as a donor-call concept (glossary §3.1).
2. **`Deleted recipient`** (legend) is a *recipient chromosome with no coverage*; the spec's "deleted recipient segment" does not exist as a literal (glossary §3.2).
3. **"component evidence"** is not a real field: the app has `class_evidence` (which side) and `component_status` (consistency across subgenomes) — glossary §3.3. If the spec's `payload$labels` emits a `tooltip.component` key, it must map to the *Component* no-call row, not to evidence.
4. **`Breadth-aware introgression mapping`** — "breadth" = fraction of window covered (not read depth); "с учётом покрытия" overstates it (glossary §3.4).
5. **`Copy state`** filter label means relative read-depth direction (gain/loss/balanced…), i.e. **копийность**, not "состояние копии" (glossary §6).
6. **`Run analysis` progress label "Preparing analysis"** also appears as the loader eyebrow/stage; keep them in sync.
7. Minor: `call_tooltip_status()` falls back to the literal `"call"` when neither class status nor confidence resolves (R/call-tooltip.R:129) — in Russian this must be `находка` (already `tooltip.call`).

## Number-formatting leak (not copy — R-side follow-up)

`format(big.mark = ",", …)` (e.g. R/app-bam.R:424,532,1304,1363,1789) and `formatC(…, format="f", digits=…)` / `call_tooltip_signed` / `call_tooltip_bp` (R/call-tooltip.R:46-68) render English "1,000" / "1.23". Russian UI needs a locale-aware formatter ("1 000", "1,23"). No Russian string in this dictionary hardcodes a formatted number.

---

## tooltip_guide.* (Tooltip help modal, `explorer_tooltip_help_ui()` R/dashboard.R:260-820)

Convention: `row(term, description)` → `<dt>`/`<dd>`; a `description` built from `shiny::tagList(...)` with inline `tags$strong(...)` value names is given ONE key with the assembled sentence and the `strong()` tokens listed in the flag column. Value tokens follow `russian-glossary.md` §6; "call/calls" = "находка/находки"; "breadth" = "широта покрытия". Modal title/button (R/dashboard.R:1213,1217) reuse `export.tooltip_guide` / `errors.modal_close`.

| key | EN | RU | source | flag |
|---|---|---|---|---|
| tooltip_guide.modal_title | Tooltip guide | Справка по подсказкам | R/dashboard.R:1213 | = export.tooltip_guide |
| tooltip_guide.modal_close | Close | Закрыть | R/dashboard.R:1217 | = errors.modal_close |
| tooltip_guide.lead | Hover a chromosome bar to inspect the nearest introgression call. The card is a call-level summary, not a per-window trace. It follows the current class, class-status, copy-state, confidence and width filters, so a filtered-out call cannot appear. | Наведите курсор на полосу хромосомы, чтобы изучить ближайшую находку интрогрессии. Карточка — это сводка на уровне находки, а не покадровый трейс. Она учитывает текущие фильтры класса, статуса класса, копийности, достоверности и ширины, поэтому отфильтрованная находка не может появиться. | R/dashboard.R:314-319 p | >1.6× |
| tooltip_guide.callout | Read top to bottom: the header identifies the call and its class confidence; Placement says where it was matched; Evidence says why; Geometry describes what changed; Position (heatmap only) adds local context. | Читайте сверху вниз: заголовок определяет находку и достоверность её класса; «Размещение» сообщает, где она была сопоставлена; «Свидетельства» — почему; «Геометрия» описывает, что изменилось; «Позиция» (только тепловая карта) добавляет локальный контекст. | R/dashboard.R:324-329 p | — |
| tooltip_guide.quick_interpretation.title | Quick interpretation | Быстрая интерпретация | R/dashboard.R:332 | — |
| tooltip_guide.quick_interpretation.intro | Use these patterns as a starting point, then check placement and geometry before treating a call as a confirmed event. | Используйте эти паттерны как отправную точку, затем проверьте размещение и геометрию, прежде чем считать находку подтверждённым событием. | R/dashboard.R:333-336 p | — |
| tooltip_guide.quick_interpretation.balanced_substitution | Balanced substitution | Сбалансированное замещение | R/dashboard.R:338 | — |
| tooltip_guide.quick_interpretation.balanced_substitution.desc | Donor excess, matching recipient loss, and a near-neutral total. Usually class substitution with copy state balanced. | Избыток донора, соответствующая потеря у реципиента и почти нейтральная сумма. Обычно класс «замещение» с копийностью «сбалансированное». | R/dashboard.R:340-346 | strong: near-neutral total / substitution / balanced |
| tooltip_guide.quick_interpretation.addition | Addition (gain) | Добавление (прирост) | R/dashboard.R:349 | — |
| tooltip_guide.quick_interpretation.addition.desc | Donor excess without reciprocal recipient loss, so total depth rises. Usually copy state gain. | Избыток донора без ответной потери у реципиента, поэтому суммарная глубина растёт. Обычно копийность «прирост». | R/dashboard.R:351-354 | strong: gain |
| tooltip_guide.quick_interpretation.deletion | Deletion (loss) | Делеция (потеря) | R/dashboard.R:357 | — |
| tooltip_guide.quick_interpretation.deletion.desc | Recipient loss without donor support, so total depth falls. Usually copy state loss. | Потеря у реципиента без поддержки донора, поэтому суммарная глубина падает. Обычно копийность «потеря». | R/dashboard.R:359-362 | strong: loss |
| tooltip_guide.quick_interpretation.complex | Complex | Сложный | R/dashboard.R:365 | — |
| tooltip_guide.quick_interpretation.complex.desc | Donor and recipient evidence disagree, or the donor and recipient spans are geometrically incompatible. Check Evidence and Geometry together. | Свидетельства донора и реципиента расходятся, либо диапазоны донора и реципиента геометрически несовместимы. Проверьте «Свидетельства» и «Геометрию» вместе. | R/dashboard.R:367-370 p | — |
| tooltip_guide.quick_interpretation.other | Other | Прочее | R/dashboard.R:373 | — |
| tooltip_guide.quick_interpretation.other.desc | A scored segment with enough signal to report but not enough support for a more specific class. Treat it as a review lead, not a confirmed event. | Оценённый сегмент с достаточным сигналом для отчёта, но без достаточной поддержки для более конкретного класса. Рассматривайте его как повод для проверки, а не как подтверждённое событие. | R/dashboard.R:375-378 p | — |
| tooltip_guide.quick_interpretation.no_host | No host | Нет реципиента | R/dashboard.R:381 | — |
| tooltip_guide.quick_interpretation.no_host.desc | The donor segment is assigned by chromosome group but there is no recipient deletion to anchor it. The class can still be confident; the host is not. | Донорный сегмент привязан по хромосомной группе, но нет делеции реципиента, чтобы закрепить его. Класс может быть уверенным; реципиент-хозяин — нет. | R/dashboard.R:383-386 p | — |
| tooltip_guide.header.title | Header | Заголовок | R/dashboard.R:390 | — |
| tooltip_guide.header.intro | The header always answers three questions: what is the call, how confident is the class, and which group and donor span does it describe? | Заголовок всегда отвечает на три вопроса: что это за находка, насколько достоверен класс и какую группу и донорный диапазон он описывает? | R/dashboard.R:392-394 p | — |
| tooltip_guide.header.class | Class | Класс | R/dashboard.R:397 | — |
| tooltip_guide.header.class.desc | substitution, addition, deletion, complex or other. The left border and class title use the same colour as the legend. | замещение, добавление, делеция, сложный или прочее. Левая граница и заголовок класса используют тот же цвет, что и легенда. | R/dashboard.R:399-409 | strong: substitution/addition/deletion/complex/other |
| tooltip_guide.header.status | Status | Статус | R/dashboard.R:413 | — |
| tooltip_guide.header.status.desc | Class confidence: confident, ambiguous or needs_review. This describes the class, not the recipient placement. | Достоверность класса: уверенный, неоднозначный или требует проверки. Это описывает класс, а не размещение у реципиента. | R/dashboard.R:415-421 | strong: confident/ambiguous/needs_review |
| tooltip_guide.header.context | Context | Контекст | R/dashboard.R:425 | — |
| tooltip_guide.header.context.desc | Sample name, chromosome group and donor span. Compare the donor span with the recipient span under Placement. | Название образца, хромосомная группа и донорный диапазон. Сравните донорный диапазон с диапазоном реципиента в разделе «Размещение». | R/dashboard.R:426-429 p | — |
| tooltip_guide.header.accent_colour | Accent colour | Акцентный цвет | R/dashboard.R:432 | — |
| tooltip_guide.header.accent_colour.desc | The card border and class title inherit the call colour from the legend; the status pill uses the class-status palette. | Рамка карточки и заголовок класса наследуют цвет находки из легенды; статусная плашка использует палитру статуса класса. | R/dashboard.R:433-437 p | — |
| tooltip_guide.placement.title | Placement | Размещение | R/dashboard.R:441 | — |
| tooltip_guide.placement.intro | Placement answers where the donor segment was matched. Placement confidence and class confidence are separate: a confident substitution can still have an ambiguous host. | Размещение отвечает, где был сопоставлен донорный сегмент. Достоверность размещения и достоверность класса независимы: уверенное замещение может иметь неоднозначного реципиента-хозяина. | R/dashboard.R:443-445 p | — |
| tooltip_guide.placement.recipient | Recipient | Реципиент | R/dashboard.R:448 | — |
| tooltip_guide.placement.recipient.desc | Assigned recipient chromosome and interval. The interval is the deletion run used for placement, not necessarily the full donor span. | Назначенная хромосома реципиента и интервал. Интервал — это делеционный участок, использованный для размещения, а не обязательно весь донорный диапазон. | R/dashboard.R:450-452 p | — |
| tooltip_guide.placement.donor | Donor | Донор | R/dashboard.R:456 | — |
| tooltip_guide.placement.donor.desc | Donor chromosome, group and interval. | Хромосома донора, группа и интервал. | R/dashboard.R:457 | — |
| tooltip_guide.placement.relation | Relation | Связь | R/dashboard.R:460 | — |
| tooltip_guide.placement.relation.desc | same group, cross-group candidate or not placed. A cross-group call is a translocation candidate, not a confirmed host. | та же группа, межгрупповой кандидат или не размещено. Межгрупповая находка — кандидат на транслокацию, а не подтверждённый реципиент-хозяин. | R/dashboard.R:462-468 | strong: same group / cross-group candidate / not placed |
| tooltip_guide.placement.assignment | Assignment | Размещение | R/dashboard.R:472 | — |
| tooltip_guide.placement.assignment.desc | assigned, cross-group candidate, ambiguous, no host or group only. This is the placement state used in the Calls table. | размещён, межгрупповой кандидат, неоднозначный, нет реципиента или только группа. Это состояние размещения, используемое в таблице «Находки». | R/dashboard.R:474-483 | strong: assigned / cross-group candidate / ambiguous / no host / group only |
| tooltip_guide.placement.margin | Margin | Отрыв | R/dashboard.R:487 | — |
| tooltip_guide.placement.margin.desc | Lead of the chosen host over the next candidate. A small margin is a reason to inspect the alternative placement. | Отрыв выбранного реципиента-хозяина над следующим кандидатом. Малый отрыв — повод проверить альтернативное размещение. | R/dashboard.R:489-491 p | — |
| tooltip_guide.placement.group_class | Group class | Класс группы | R/dashboard.R:494 | — |
| tooltip_guide.placement.group_class.desc | Shown only when recipient placement changes the call class, for example addition to substitution. | Показывается только когда размещение у реципиента меняет класс находки, например «добавление» на «замещение». | R/dashboard.R:496-498 p | — |
| tooltip_guide.placement.alternative | Alternative | Альтернатива | R/dashboard.R:501 | — |
| tooltip_guide.placement.alternative.desc | Another possible host shown when the assignment is ambiguous. Compare its group, interval and margin before choosing a host. | Другой возможный реципиент-хозяин, показываемый когда размещение неоднозначно. Сравните его группу, интервал и отрыв, прежде чем выбрать реципиента-хозяина. | R/dashboard.R:502-505 p | — |
| tooltip_guide.placement.reading | Reading placement | Как читать размещение | R/dashboard.R:509 | — |
| tooltip_guide.placement.reading.desc | Strongest: same group and assigned. Cross-group is a candidate only. No host means there is no recipient loss to anchor the call. Group only is a fallback when no contig-level host is available. | Сильнее всего: та же группа и «размещён». Межгрупповой — только кандидат. «Нет реципиента» означает отсутствие потери у реципиента для закрепления находки. «Только группа» — запасной вариант, когда нет реципиента-хозяина на уровне контига. | R/dashboard.R:510-514 p | — |
| tooltip_guide.evidence.title | Evidence | Свидетельства | R/dashboard.R:519 | — |
| tooltip_guide.evidence.intro | Evidence compares donor and recipient breadth in the same windows. This is the core reason for the call; the numbers are signed residuals, not raw read counts. | Свидетельства сравнивают широту покрытия донора и реципиента в одних и тех же окнах. Это основная причина находки; числа — это знаковые остатки, а не сырые счётчики прочтений. | R/dashboard.R:521-523 p | — |
| tooltip_guide.evidence.score | Score | Оценка | R/dashboard.R:526 | — |
| tooltip_guide.evidence.score.desc | Paired donor/recipient statistic for the segment. Larger positive values mean stronger donor excess. Thresholds are calibrated for the window count and segment width. | Парная статистика донор/реципиент для сегмента. Большие положительные значения означают более сильный избыток донора. Пороги откалиброваны под число окон и ширину сегмента. | R/dashboard.R:528-530 p | — |
| tooltip_guide.evidence.donor | Donor | Донор | R/dashboard.R:534 | — |
| tooltip_guide.evidence.donor.desc | Mean donor-side residual. Positive means donor gain; negative means donor loss. | Средний остаток со стороны донора. Положительный — прирост донора; отрицательный — потеря донора. | R/dashboard.R:535 | — |
| tooltip_guide.evidence.recipient | Recipient | Реципиент | R/dashboard.R:538 | — |
| tooltip_guide.evidence.recipient.desc | Mean recipient-side residual. Positive means recipient gain; negative means recipient loss. | Средний остаток со стороны реципиента. Положительный — прирост реципиента; отрицательный — потеря реципиента. | R/dashboard.R:539 | — |
| tooltip_guide.evidence.total | Total | Суммарно | R/dashboard.R:542 | — |
| tooltip_guide.evidence.total.desc | Relative total-depth residual. Near zero means balanced; positive means total gain; negative means total loss. This is what separates a replacement from an addition or deletion. | Относительный остаток суммарной глубины. Около нуля — сбалансированное; положительный — суммарный прирост; отрицательный — суммарная потеря. Это отличает замещение от добавления или делеции. | R/dashboard.R:544-547 p | — |
| tooltip_guide.evidence.colour | Colour | Цвет | R/dashboard.R:551 | — |
| tooltip_guide.evidence.colour.desc | Metrics use two decimals: green is positive, red is negative and muted is missing. The sign matters more than the absolute size for interpretation. | Метрики используют два знака после запятой: зелёный — положительное, красный — отрицательное, приглушённый — отсутствует. Для интерпретации знак важнее абсолютного значения. | R/dashboard.R:552-555 p | — |
| tooltip_guide.evidence.copy_state | Copy state | Копийность | R/dashboard.R:559 | — |
| tooltip_guide.evidence.copy_state.desc | balanced, gain, loss, compound or neutral. This is a relative read-depth direction, not a quantitative measure. | сбалансированное, прирост, потеря, составное или нейтральное. Это относительное направление изменения глубины прочтений, а не количественная мера. | R/dashboard.R:561-571 | strong: balanced/gain/loss/compound/neutral |
| tooltip_guide.evidence.evidence_source | Evidence source | Источник свидетельств | R/dashboard.R:575 | — |
| tooltip_guide.evidence.evidence_source.desc | donor, recipient, donor+recipient or none. It says which side drove the class. | донор, реципиент, донор+реципиент или нет. Указывает, какая сторона определила класс. | R/dashboard.R:577-584 | strong: donor/recipient/donor+recipient/none |
| tooltip_guide.evidence.typical_pattern | Typical pattern | Типичный паттерн | R/dashboard.R:588 | — |
| tooltip_guide.evidence.typical_pattern.desc | Substitution: donor positive, recipient negative, total near zero. Addition: donor positive, recipient near zero, total positive. Deletion: donor near zero, recipient negative, total negative. Complex: conflicting signs or a width mismatch. | Замещение: донор положительный, реципиент отрицательный, сумма около нуля. Добавление: донор положительный, реципиент около нуля, сумма положительная. Делеция: донор около нуля, реципиент отрицательный, сумма отрицательная. Сложный: противоречивые знаки или несовпадение ширин. | R/dashboard.R:590-594 p | — |
| tooltip_guide.status.title | Status chips | Статусные плашки | R/dashboard.R:599 | — |
| tooltip_guide.status.intro | Chips appear only when their condition holds. They are shortcuts into the call table and the full result columns. | Плашки появляются только когда выполняется их условие. Это ссылки в таблицу находок и полные столбцы результата. | R/dashboard.R:601-602 p | — |
| tooltip_guide.status.class_confidence | class high / moderate / low | класс: высокая / средняя / низкая | R/dashboard.R:605 | §6 + confidence enum |
| tooltip_guide.status.class_confidence.desc | Class confidence reported by the classifier. | Достоверность класса, сообщаемая классификатором. | R/dashboard.R:606 | — |
| tooltip_guide.status.evidence_sides | donor / recipient / donor+recipient / none | донор / реципиент / донор+реципиент / нет | R/dashboard.R:609 | — |
| tooltip_guide.status.evidence_sides.desc | Which evidence sides supported the class. | Какие стороны свидетельств поддержали класс. | R/dashboard.R:610 | — |
| tooltip_guide.status.copy_states | balanced / gain / loss / compound / neutral | сбалансированное / прирост / потеря / составное / нейтральное | R/dashboard.R:613 | — |
| tooltip_guide.status.copy_states.desc | The copy-state hypothesis for the segment. | Гипотеза копийности для сегмента. | R/dashboard.R:614 | — |
| tooltip_guide.status.total | total gain / total loss | суммарный прирост / суммарная потеря | R/dashboard.R:617 | — |
| tooltip_guide.status.total.desc | Total-depth evidence is not neutral. | Свидетельство суммарной глубины не нейтрально. | R/dashboard.R:618 | — |
| tooltip_guide.status.components | components heterogeneous | компоненты неоднородны | R/dashboard.R:621 | = tooltip.chip_components_heterogeneous |
| tooltip_guide.status.components.desc | Component spread or donor-component support is uneven; the group-level score may hide a component-specific event. | Разброс компонентов или поддержка донорных компонентов неравномерны; оценка на уровне группы может скрывать событие, специфичное для компонента. | R/dashboard.R:622-625 p | — |
| tooltip_guide.status.low_mappability | low mappability | низкая картированность | R/dashboard.R:628 | = tooltip.chip_low_mappability |
| tooltip_guide.status.low_mappability.desc | BAM-derived high-MAPQ fraction is below 0.5; placement is less reliable. | Доля прочтений с высоким MAPQ (из BAM) ниже 0,5; размещение менее надёжно. | R/dashboard.R:629-632 p | — |
| tooltip_guide.status.not_passing | not passing | не проходит | R/dashboard.R:635 | = tooltip.chip_not_passing |
| tooltip_guide.status.not_passing.desc | The call is shown for review but did not pass the configured filters. | Находка показана для проверки, но не прошла настроенные фильтры. | R/dashboard.R:636 | — |
| tooltip_guide.geometry.title | Geometry | Геометрия | R/dashboard.R:640 | — |
| tooltip_guide.geometry.intro | Geometry checks whether the donor and recipient spans are compatible. It can turn a substitution into a complex call even when the depth evidence looks clean. | Геометрия проверяет, совместимы ли диапазоны донора и реципиента. Она может превратить замещение в сложную находку, даже когда свидетельства глубины выглядят чистыми. | R/dashboard.R:642-644 p | — |
| tooltip_guide.geometry.donor_width | Donor width | Ширина донора | R/dashboard.R:646 | — |
| tooltip_guide.geometry.donor_width.desc | Width of the donor segment. | Ширина донорного сегмента. | R/dashboard.R:646 | — |
| tooltip_guide.geometry.recipient_run | Recipient run | Делеция реципиента | R/dashboard.R:648 | = tooltip.recipient_run |
| tooltip_guide.geometry.recipient_run.desc | Width of the recipient deletion run used for placement. | Ширина делеционного участка реципиента, использованного для размещения. | R/dashboard.R:649 | — |
| tooltip_guide.geometry.overlap | Overlap | Перекрытие | R/dashboard.R:652 | = tooltip.overlap |
| tooltip_guide.geometry.overlap.desc | Reciprocal overlap of donor and recipient intervals in contig-fraction space. Low overlap weakens the pairing. | Взаимное перекрытие интервалов донора и реципиента в пространстве долей контига. Низкое перекрытие ослабляет парность. | R/dashboard.R:653-656 p | — |
| tooltip_guide.geometry.width_ratio | Width ratio | Отношение ширин | R/dashboard.R:659 | = tooltip.width_ratio |
| tooltip_guide.geometry.width_ratio.desc | Donor fraction divided by recipient-run fraction. A value near 1 supports a balanced replacement; values far from 1 suggest an addition, deletion or complex event. | Доля донора, делённая на долю делеционного участка реципиента. Значение около 1 поддерживает сбалансированное замещение; значения, далёкие от 1, указывают на добавление, делецию или сложное событие. | R/dashboard.R:660-664 p | — |
| tooltip_guide.geometry.interpretation | Interpretation | Интерпретация | R/dashboard.R:667 | — |
| tooltip_guide.geometry.interpretation.desc | Overlap below 0.5 or a width ratio far from 1 is a reason to prefer complex or review the recipient assignment. | Перекрытие ниже 0,5 или отношение ширин, далёкое от 1, — повод предпочесть «сложный» или пересмотреть размещение у реципиента. | R/dashboard.R:668-671 p | — |
| tooltip_guide.position.title | Position and footer | Позиция и нижний колонтитул | R/dashboard.R:675 | — |
| tooltip_guide.position.intro | Position is available in Heatmap mode; Calls mode ends with the footer instead. | Позиция доступна в режиме «Тепловая карта»; режим «Находки» заканчивается нижним колонтитулом. | R/dashboard.R:677-678 p | — |
| tooltip_guide.position.position | Position | Позиция | R/dashboard.R:681 | = tooltip.position |
| tooltip_guide.position.position.desc | Hovered position along the chromosome, in bp, kb or Mb. | Позиция под курсором вдоль хромосомы, в bp, kb или Mb. | R/dashboard.R:682 | — |
| tooltip_guide.position.score | Score | Оценка | R/dashboard.R:684 | — |
| tooltip_guide.position.score.desc | Heatmap score at that relative position. | Оценка тепловой карты в этой относительной позиции. | R/dashboard.R:684 | — |
| tooltip_guide.position.mode | Mode | Режим | R/dashboard.R:686 | — |
| tooltip_guide.position.mode.desc | Heatmap shows the score profile; Calls shows coloured call segments. | «Тепловая карта» показывает профиль оценок; «Находки» показывает цветные сегменты находок. | R/dashboard.R:687 | — |
| tooltip_guide.position.click | Click | Клик | R/dashboard.R:690 | — |
| tooltip_guide.position.click.desc | Click a chromosome to filter its group; click the background to clear the filter. | Кликните хромосому, чтобы отфильтровать её группу; кликните фон, чтобы снять фильтр. | R/dashboard.R:691-694 p | — |
| tooltip_guide.no_call.title | No-call card | Карточка без находки | R/dashboard.R:698 | — |
| tooltip_guide.no_call.intro | When the pointer is not over a call, the tooltip describes the chromosome instead of a call. | Когда курсор не над находкой, подсказка описывает хромосому вместо находки. | R/dashboard.R:700-701 p | — |
| tooltip_guide.no_call.header | Header | Заголовок | R/dashboard.R:704 | — |
| tooltip_guide.no_call.header.desc | Chromosome state: recipient chromosome, deleted recipient or added donor. | Состояние хромосомы: хромосома реципиента, реципиент с делецией или добавленный донор. | R/dashboard.R:705 | — |
| tooltip_guide.no_call.chromosome_block | Chromosome block | Блок хромосомы | R/dashboard.R:708 | — |
| tooltip_guide.no_call.chromosome_block.desc | Component, state and configured length. | Компонент, состояние и настроенная длина. | R/dashboard.R:709 | — |
| tooltip_guide.no_call.partial_addition | Partial addition | Частичное добавление | R/dashboard.R:712 | — |
| tooltip_guide.no_call.partial_addition.desc | A donor chromosome bar appears when a partial addition covers at least 20% of that donor. Only the added part is blue; the remainder is an empty dashed donor deletion. | Полоса донорной хромосомы появляется, когда частичное добавление покрывает не менее 20% этого донора. Только добавленная часть синяя; остальное — пустая пунктирная делеция донора. | R/dashboard.R:713-716 p | — |
| tooltip_guide.no_call.added_donor_fill | Added donor fill | Заливка добавленного донора | R/dashboard.R:720 | — |
| tooltip_guide.no_call.added_donor_fill.desc | Calls mode uses the light-blue added-donor state colour. Heatmap mode uses the heatmap zero-coverage colour outside the score bins. Calls drawn on the bar keep their class colour. | Режим «Находки» использует светло-синий цвет состояния «добавленный донор». Режим «Тепловая карта» использует цвет нулевого покрытия тепловой карты вне бинов оценок. Находки, нарисованные на полосе, сохраняют цвет своего класса. | R/dashboard.R:721-725 p | — |
| tooltip_guide.no_call.call_block | Call block | Блок находки | R/dashboard.R:729 | — |
| tooltip_guide.no_call.call_block.desc | Calls mode only: says that no assigned donor segment is present at the hovered position. | Только режим «Находки»: сообщает, что в позиции под курсором нет назначенного донорного сегмента. | R/dashboard.R:730-732 p | — |
| tooltip_guide.no_call.position_block | Position block | Блок позиции | R/dashboard.R:736 | — |
| tooltip_guide.no_call.position_block.desc | Heatmap mode only: the hovered position and score. | Только режим «Тепловая карта»: позиция под курсором и оценка. | R/dashboard.R:737 | — |
| tooltip_guide.hide_deletions.title | Hide deletions | Скрыть делеции | R/dashboard.R:741 | = karyotype.hide_deletions |
| tooltip_guide.hide_deletions.intro | The switch beside the Plot buttons removes deletions from the karyotype instead of drawing them as dashed or red segments. The counter beside it reports how many segments were omitted in the current view. | Переключатель рядом с кнопками «Вид» убирает делеции из кариотипа, вместо того чтобы рисовать их пунктирными или красными сегментами. Счётчик рядом сообщает, сколько сегментов опущено в текущем виде. | R/dashboard.R:743-746 p | — |
| tooltip_guide.hide_deletions.partial_deletion | Partial deletion | Частичная делеция | R/dashboard.R:749 | — |
| tooltip_guide.hide_deletions.partial_deletion.desc | The deleted interval is spliced out and the remaining recipient chromosome closes the gap and becomes shorter. | Делетированный интервал вырезается, и оставшаяся хромосома реципиента закрывает разрыв и становится короче. | R/dashboard.R:750-753 p | — |
| tooltip_guide.hide_deletions.complete_deletion | Complete deletion | Полная делеция | R/dashboard.R:756 | — |
| tooltip_guide.hide_deletions.complete_deletion.desc | The recipient chromosome is not drawn at all. | Хромосома реципиента не рисуется вовсе. | R/dashboard.R:757 | — |
| tooltip_guide.hide_deletions.substitution | Substitution | Замещение | R/dashboard.R:760 | — |
| tooltip_guide.hide_deletions.substitution.desc | The added donor chromosome takes the place of the deleted recipient chromosome. | Добавленная донорная хромосома занимает место делетированной хромосомы реципиента. | R/dashboard.R:761-764 p | — |
| tooltip_guide.unplaced.title | Unplaced calls | Неразмещённые находки | R/dashboard.R:768 | — |
| tooltip_guide.unplaced.intro | A call is unplaced when it cannot be drawn on a recipient chromosome or an added donor bar. With Hide unplaced off, its donor contig is drawn as a neutral unplaced bar and its calls remain in the table. With Hide unplaced on, the bar and table rows are removed. Donor-only additions are placed on added donor bars and stay visible. The counter reports unplaced calls in the current view; full exports keep every call. | Находка неразмещена, когда её нельзя нарисовать на хромосоме реципиента или на полосе добавленного донора. При выключенном «Скрыть неразмещённые» её донорный контиг рисуется как нейтральная неразмещённая полоса, а находки остаются в таблице. При включённом «Скрыть неразмещённые» полоса и строки таблицы удаляются. Чисто донорные добавления размещаются на полосах добавленного донора и остаются видимыми. Счётчик сообщает о неразмещённых находках в текущем виде; полный экспорт сохраняет каждую находку. | R/dashboard.R:770-776 p | >1.6× |
| tooltip_guide.caveats.title | Caveats | Оговорки | R/dashboard.R:780 | — |
| tooltip_guide.caveats.relative | Relative, not absolute | Относительное, а не абсолютное | R/dashboard.R:783 | — |
| tooltip_guide.caveats.relative.desc | Copy state describes a relative read-depth direction. It is not a quantitative measure. | Копийность описывает относительное направление изменения глубины прочтений. Это не количественная мера. | R/dashboard.R:784-786 p | — |
| tooltip_guide.caveats.class_vs_placement | Class vs placement | Класс против размещения | R/dashboard.R:790 | — |
| tooltip_guide.caveats.class_vs_placement.desc | Class status describes confidence in the label; the Assignment row describes confidence in the host. | Статус класса описывает достоверность метки; строка «Размещение» описывает достоверность реципиента-хозяина. | R/dashboard.R:792-793 p | — |
| tooltip_guide.caveats.cross_group | Cross-group | Межгрупповой | R/dashboard.R:797 | — |
| tooltip_guide.caveats.cross_group.desc | A cross-group result is a translocation candidate, not a confirmed donor-recipient pair. | Межгрупповой результат — кандидат на транслокацию, а не подтверждённая пара донор-реципиент. | R/dashboard.R:798-800 p | — |
| tooltip_guide.caveats.filters | Filters | Фильтры | R/dashboard.R:804 | — |
| tooltip_guide.caveats.filters.desc | The card follows the visible calls. Clear a filter to inspect calls that are currently hidden. | Карточка следует за видимыми находками. Снимите фильтр, чтобы изучить находки, которые сейчас скрыты. | R/dashboard.R:805-807 p | — |
| tooltip_guide.caveats.heatmap_score | Heatmap score | Оценка тепловой карты | R/dashboard.R:811 | — |
| tooltip_guide.caveats.heatmap_score.desc | A heatmap value is local context, not call-level evidence. Use the Score, Donor, Recipient and Total metrics for the call itself. | Значение тепловой карты — это локальный контекст, а не свидетельство на уровне находки. Используйте метрики «Оценка», «Донор», «Реципиент» и «Суммарно» для самой находки. | R/dashboard.R:813-815 p | — |

---

## Coverage note — `tooltip_guide.*`

The `tooltip_guide.*` block is now **complete**: every section heading, intro, `<dt>` row label and `<dd>` description in `explorer_tooltip_help_ui()` (R/dashboard.R:260-820) is enumerated, plus the modal title/button (R/dashboard.R:1213,1217). Row descriptions built from `shiny::tagList(...)` with inline `tags$strong(...)` value names are given one key per sentence, with the `strong()` value tokens listed in the flag column so the engineer can re-wrap the translated tokens in `<strong>` mechanically.

Resolved term conflicts / decisions made here (no source edit):
1. **Confidence enum `high / moderate / low`** — not in glossary §6. Rendered as feminine "высокая / средняя / низкая" to agree with "достоверность". This also implies the *confidence filter* choices "High"/"Moderate"/"Low" (`R/dashboard.R:99-100`) are currently **missing** from the `results.*` block and should be added with the same values — flagged for the engineer, not added here (out of this append's scope).
2. **"host"** in placement rows is rendered "реципиент-хозяин" (the recipient chromosome anchoring the donor segment), consistent with `tooltip.no_host` = "нет реципиента".
3. **"total"** (metric/chip) rendered "суммарно"/"суммарный" to avoid collision with "Total" meaning "итог"; kept consistent with `tooltip.chip_total`.

Untranslated / kept-as-is by rule: none in this block except `bp`/`kb`/`Mb`, `MAPQ`, `Heatmap`/`Calls` (as mode names inside prose), and `donor+recipient` (kept with the "+" as in the app's enum). No exported-bundle strings were touched.

---

## Remaining UI strings (completeness pass)

Verified against the working tree (line numbers are current). Grouped under `errors.*`, `sidebar.upload_help.*`, `log.*`, `call_table.*`. Rules as in the rest of the file: infinitive buttons, sentence case, ё, glossary terms (находка, широта покрытия, неразмещённый донорный контиг, копийность), `%s`/`%d` order preserved.

| key | EN | RU | source | flag |
|---|---|---|---|---|
| errors.modal_default_message | The BAM is too large for this server. | BAM слишком большой для этого сервера. | R/app-bam.R:2536 | — |
| errors.modal_detail | This BAM is %s; this server accepts at most %s. | Этот BAM имеет размер %s; сервер принимает не более %s. | R/app-bam.R:2540 | — |
| errors.modal_rejected_note | The upload was rejected before staging or analysis; no BAM data was processed. | Загрузка была отклонена до размещения или анализа; данные BAM не обрабатывались. | R/app-bam.R:2575-2576 p | — |
| errors.example_resolution | The bundled example supports resolution(s): %s. | Встроенный пример поддерживает разрешения: %s. | R/app-bam.R:5227 | — |
| sidebar.upload_help.kept_desc | Everything MaMaMIA reads: FLAG, RNAME, POS, MAPQ, CIGAR, @HD/@SQ/@RG and RG tags. | Всё, что читает MaMaMIA: FLAG, RNAME, POS, MAPQ, CIGAR, @HD/@SQ/@RG и теги RG. | R/app-bam.R:2682-2683 p | — |
| sidebar.upload_help.removed_desc | QUAL, read names, optional tags except RG, and @PG/@CO header lines. | QUAL, имена прочтений, необязательные теги кроме RG, а также строки заголовка @PG/@CO. | R/app-bam.R:2691-2692 p | — |
| sidebar.upload_help.optional_desc | SEQ with --drop-seq. Read length still comes from CIGAR, so the analysis is unchanged. | SEQ с --drop-seq. Длина прочтений по-прежнему берётся из CIGAR, поэтому анализ не меняется. | R/app-bam.R:2700-2701 p | — |
| sidebar.upload_help.from | " from " | " с " | R/app-bam.R:2642 | leading+trailing space preserved |
| sidebar.upload_help.run_tail | " the command below on the large BAM." | " команду ниже для большого BAM." | R/app-bam.R:2655 | leading space preserved |
| sidebar.upload_help.upload_tail | " the new, smaller BAM here." | " новый, уменьшенный BAM сюда." | R/app-bam.R:2661 | leading space preserved |
| log.coverage_no_reads | Sample '%s' has no reads on the %s targets. | У образца '%s' нет прочтений на целях %s. | R/app-bam.R:2062 | — |
| log.coverage_near_zero | Sample '%s' has near-zero coverage on the %s side (%s reads over %d windows). | У образца '%s' почти нулевое покрытие на стороне %s (%s прочтений в %d окнах). | R/app-bam.R:2067 | — |
| log.preflight_warning_no_rg | No @RG SM value; the sample name is derived from the file name. | Нет значения @RG SM; имя образца выведено из имени файла. | R/preflight.R:68 | — |
| log.preflight_warning_targets | The BAM header contains %d target(s) not in assembly '%s'; they will be ignored. | Заголовок BAM содержит целей вне сборки '%s': %d; они будут проигнорированы. | R/preflight.R:77 | plural rework (number last) |
| log.preflight_warning_bai | A BAI index is present but ignored; MaMaMIA requires a CSI index. | Обнаружен BAI-индекс; он игнорируется — MaMaMIA требуется CSI-индекс. | R/preflight.R:86 | — |
| log.preflight_warning_csi_create | A CSI index will be created next to the BAM. | Рядом с BAM будет создан CSI-индекс. | R/preflight.R:91 | — |
| log.preflight_warning_csi_rebuild | The existing CSI index is stale or corrupt and will be rebuilt. | Существующий CSI-индекс устарел или повреждён и будет пересобран. | R/preflight.R:92 | — |
| log.preflight_warning_csi_stage | The BAM directory is not writable; the BAM will be staged in session temp and indexed there. | Каталог BAM недоступен для записи; BAM будет размещён во временном каталоге сеанса и проиндексирован там. | R/preflight.R:93 | — |

### call_table.* — Calls-tab column labels (user-visible DT headers, `R/call-schema.R:111-129`)

The Calls tab renders these via `bam_call_table_header()` → `bam_call_column_metadata()` → `call_schema_labels()`. They are in scope (not export-only). Value-enum names inside labels stay Latin per rule (bp, MAPQ, ID, Z, SE).

| key | EN | RU | source | flag |
|---|---|---|---|---|
| call_table.label.sample | Sample | Образец | R/call-schema.R:112 | — |
| call_table.label.group | Group | Группа | R/call-schema.R:112 | — |
| call_table.label.start_iid | Start window | Начальное окно | R/call-schema.R:112 | — |
| call_table.label.end_iid | End window | Конечное окно | R/call-schema.R:112 | — |
| call_table.label.start_bp | Group start (bp) | Начало группы (п.н.) | R/call-schema.R:113 | — |
| call_table.label.end_bp | Group end (bp) | Конец группы (п.н.) | R/call-schema.R:113 | — |
| call_table.label.width_bp | Width (bp) | Ширина (п.н.) | R/call-schema.R:113 | — |
| call_table.label.score_mean | Score | Оценка | R/call-schema.R:113 | — |
| call_table.label.donor_mean | Donor mean | Среднее донора | R/call-schema.R:113 | — |
| call_table.label.recipient_mean | Recipient mean | Среднее реципиента | R/call-schema.R:113 | — |
| call_table.label.total_mean | Total mean | Суммарное среднее | R/call-schema.R:114 | — |
| call_table.label.donor_component_spread | Donor spread | Разброс донора | R/call-schema.R:114 | — |
| call_table.label.recipient_component_spread | Recipient spread | Разброс реципиента | R/call-schema.R:114 | — |
| call_table.label.n_donor_supported | Donor sets | Наборы донора | R/call-schema.R:114 | — |
| call_table.label.mapq_fraction | MAPQ fraction | Доля MAPQ | R/call-schema.R:114 | — |
| call_table.label.n_windows | Windows | Окна | R/call-schema.R:115 | — |
| call_table.label.method | Method | Метод | R/call-schema.R:115 | — |
| call_table.label.mode | Mode | Режим | R/call-schema.R:115 | — |
| call_table.label.metric | Metric | Метрика | R/call-schema.R:115 | — |
| call_table.label.length_bin | Length bin | Бин длины | R/call-schema.R:115 | — |
| call_table.label.threshold | Threshold | Порог | R/call-schema.R:115 | — |
| call_table.label.pass | Pass | Проходит | R/call-schema.R:115 | — |
| call_table.label.class | Class | Класс | R/call-schema.R:115 | — |
| call_table.label.class_evidence | Class evidence | Свидетельства класса | R/call-schema.R:115 | — |
| call_table.label.copy_state | Copy state | Копийность | R/call-schema.R:115 | — |
| call_table.label.total_evidence | Total evidence | Суммарные свидетельства | R/call-schema.R:115 | — |
| call_table.label.class_confidence | Class confidence | Достоверность класса | R/call-schema.R:115 | — |
| call_table.label.class_status | Status | Статус | R/call-schema.R:115 | — |
| call_table.label.component_status | Component status | Статус компонентов | R/call-schema.R:115 | — |
| call_table.label.mappability_status | Mappability | Картируемость | R/call-schema.R:116 | — |
| call_table.label.confidence | Confidence | Достоверность | R/call-schema.R:116 | — |
| call_table.label.class_group | Group class | Класс группы | R/call-schema.R:116 | — |
| call_table.label.recipient_component | Recipient set | Набор реципиента | R/call-schema.R:116 | — |
| call_table.label.recipient_chr_id | Recipient ID | ID реципиента | R/call-schema.R:116 | — |
| call_table.label.recipient_chr_name | Recipient | Реципиент | R/call-schema.R:116 | — |
| call_table.label.recipient_group | Recipient group | Группа реципиента | R/call-schema.R:116 | — |
| call_table.label.recipient_start_bp | Recipient start (bp) | Начало реципиента (п.н.) | R/call-schema.R:117 | — |
| call_table.label.recipient_end_bp | Recipient end (bp) | Конец реципиента (п.н.) | R/call-schema.R:117 | — |
| call_table.label.recipient_run_start_bp | Depletion start (bp) | Начало делеции (п.н.) | R/call-schema.R:118 | — |
| call_table.label.recipient_run_end_bp | Depletion end (bp) | Конец делеции (п.н.) | R/call-schema.R:118 | — |
| call_table.label.recipient_deletion_width_bp | Depletion width (bp) | Ширина делеции (п.н.) | R/call-schema.R:118 | — |
| call_table.label.recipient_expected_breadth | Expected breadth | Ожидаемая широта покрытия | R/call-schema.R:119 | >1.6× |
| call_table.label.recipient_observed_breadth | Observed breadth | Наблюдаемая широта покрытия | R/call-schema.R:119 | >1.6× |
| call_table.label.recipient_mapq_fraction | Recipient MAPQ | MAPQ реципиента | R/call-schema.R:119 | — |
| call_table.label.recipient_mean_depletion | Mean depletion | Средняя делеция | R/call-schema.R:119 | — |
| call_table.label.recipient_depletion_se | Depletion SE | СО делеции | R/call-schema.R:119 | — |
| call_table.label.recipient_depletion_lower | Depletion lower | Нижняя граница делеции | R/call-schema.R:119 | — |
| call_table.label.recipient_depletion_z | Depletion Z | Z делеции | R/call-schema.R:119 | — |
| call_table.label.recipient_width_overlap | Width overlap | Перекрытие ширин | R/call-schema.R:119 | — |
| call_table.label.recipient_fraction_ratio | Width ratio | Отношение ширин | R/call-schema.R:119 | — |
| call_table.label.copy_state_group | Group copy state | Копийность группы | R/call-schema.R:120 | — |
| call_table.label.recipient_weight | Recipient weight | Вес реципиента | R/call-schema.R:120 | — |
| call_table.label.recipient_margin | Recipient margin | Отрыв реципиента | R/call-schema.R:120 | — |
| call_table.label.recipient_assignment | Placement | Размещение | R/call-schema.R:120 | — |
| call_table.label.recipient_cross_group | Cross-group | Межгрупповой | R/call-schema.R:120 | — |
| call_table.label.recipient_reliable | Recipient reliable | Реципиент надёжен | R/call-schema.R:120 | — |
| call_table.label.recipient_strong | Recipient strong | Реципиент сильный | R/call-schema.R:120 | — |
| call_table.label.recipient_alternative | Alternative recipient | Альтернативный реципиент | R/call-schema.R:121 | — |
| call_table.label.recipient_alternative_score | Alternative score | Оценка альтернативы | R/call-schema.R:121 | — |
| call_table.label.donor_component | Donor set | Набор донора | R/call-schema.R:121 | — |
| call_table.label.donor_chr_id | Donor ID | ID донора | R/call-schema.R:121 | — |
| call_table.label.donor_chr_name | Donor | Донор | R/call-schema.R:121 | — |
| call_table.label.donor_start_bp | Donor start (bp) | Начало донора (п.н.) | R/call-schema.R:121 | — |
| call_table.label.donor_end_bp | Donor end (bp) | Конец донора (п.н.) | R/call-schema.R:121 | — |

### call_table.* — column descriptions (hover tooltips, `R/call-schema.R:168-232`)

| key | EN | RU | source | flag |
|---|---|---|---|---|
| call_table.desc.sample | Sample name read from the BAM @RG SM tag. | Имя образца из тега BAM @RG SM. | R/call-schema.R:169 | — |
| call_table.desc.group | Chromosome homology group used for scoring and segmentation, for example 7. | Хромосомная группа гомологии, используемая для оценки и сегментации, например 7. | R/call-schema.R:170 | — |
| call_table.desc.start_iid | First analysis window index in the call, zero-based within the group. | Индекс первого окна анализа в находке, отсчёт с нуля внутри группы. | R/call-schema.R:171 | — |
| call_table.desc.end_iid | Last analysis window index in the call, zero-based and inclusive. | Индекс последнего окна анализа в находке, отсчёт с нуля, включительно. | R/call-schema.R:172 | — |
| call_table.desc.start_bp | Call start in group-level coordinates, in base pairs, scaled to the donor group length. | Начало находки в координатах группы, в парах оснований, с масштабированием к длине группы донора. | R/call-schema.R:173 | — |
| call_table.desc.end_bp | Call end in group-level coordinates, in base pairs. | Конец находки в координатах группы, в парах оснований. | R/call-schema.R:174 | — |
| call_table.desc.width_bp | Call width in base pairs, calculated as end_bp minus start_bp. | Ширина находки в парах оснований (end_bp − start_bp). | R/call-schema.R:175 | — |
| call_table.desc.score_mean | Mean introgression score across the call; higher values mean stronger donor evidence. | Средняя оценка интрогрессии по находке; больше значения — сильнее свидетельства донора. | R/call-schema.R:176 | — |
| call_table.desc.donor_mean | Mean donor residual across the call; positive values mean donor enrichment. | Средний остаток донора по находке; положительные значения — обогащение донора. | R/call-schema.R:177 | — |
| call_table.desc.recipient_mean | Mean recipient residual across the call; negative values mean recipient loss. | Средний остаток реципиента по находке; отрицательные значения — потеря реципиента. | R/call-schema.R:178 | — |
| call_table.desc.total_mean | Mean total-depth residual across the call; positive values mean a total gain and negative values a total loss. | Средний остаток суммарной глубины по находке; положительные значения — суммарный прирост, отрицательные — суммарная потеря. | R/call-schema.R:179 | — |
| call_table.desc.donor_component_spread | Spread of donor residuals across donor copy sets; large values flag component heterogeneity. | Разброс остатков донора по наборам копий донора; большие значения указывают на неоднородность компонентов. | R/call-schema.R:180 | — |
| call_table.desc.recipient_component_spread | Spread of recipient residuals across recipient copy sets; large values flag component heterogeneity. | Разброс остатков реципиента по наборам копий реципиента; большие значения указывают на неоднородность компонентов. | R/call-schema.R:181 | — |
| call_table.desc.n_donor_supported | Number of donor copy sets that support the call. | Число наборов копий донора, поддерживающих находку. | R/call-schema.R:182 | — |
| call_table.desc.mapq_fraction | Mean fraction of high-MAPQ reads across the call; values below 0.5 flag low mappability. | Средняя доля прочтений с высоким MAPQ по находке; значения ниже 0,5 указывают на низкую картированность. | R/call-schema.R:183 | — |
| call_table.desc.n_windows | Number of analysis windows covered by the call. | Число окон анализа, покрытых находкой. | R/call-schema.R:184 | — |
| call_table.desc.method | Scoring method: paired_lr for a paired donor/recipient log ratio, or donor_only. | Метод оценки: paired_lr — парный лог-коэффициент донор/реципиент, или donor_only. | R/call-schema.R:185 | — |
| call_table.desc.mode | Analysis mode: single for one sample, or multi for a cohort. | Режим анализа: single — один образец, multi — когорта. | R/call-schema.R:186 | — |
| call_table.desc.metric | Coverage metric used for scoring: breadth, count or combined. | Метрика покрытия для оценки: breadth, count или combined. | R/call-schema.R:187 | — |
| call_table.desc.length_bin | Segment-length bin used to select the calibration threshold. | Бин длины сегмента, используемый для выбора порога калибровки. | R/call-schema.R:188 | — |
| call_table.desc.threshold | Calibrated score threshold for this call's mode and length bin. | Откалиброванный порог оценки для режима и бина длины этой находки. | R/call-schema.R:189 | — |
| call_table.desc.pass | TRUE when score_mean exceeds the calibrated threshold. | TRUE, если score_mean превышает откалиброванный порог. | R/call-schema.R:190 | — |
| call_table.desc.class | Final call class after recipient assignment: substitution, addition, deletion, complex or other. | Итоговый класс находки после размещения у реципиента: substitution, addition, deletion, complex или other. | R/call-schema.R:191 | — |
| call_table.desc.class_evidence | Evidence that supported the call: donor, recipient, donor+recipient, total imbalance or none. | Свидетельства, поддержавшие находку: donor, recipient, donor+recipient, дисбаланс суммы или none. | R/call-schema.R:192 | — |
| call_table.desc.copy_state | Copy-number interpretation: balanced, gain, loss, compound or neutral. | Интерпретация числа копий: balanced, gain, loss, compound или neutral. | R/call-schema.R:193 | — |
| call_table.desc.total_evidence | Total-depth evidence: gain, loss, neutral, or blank when unavailable. | Свидетельство суммарной глубины: gain, loss, neutral или пусто, если недоступно. | R/call-schema.R:194 | — |
| call_table.desc.class_confidence | Class-level confidence: high, moderate or low. | Достоверность на уровне класса: high, moderate или low. | R/call-schema.R:195 | — |
| call_table.desc.class_status | Review status: confident, ambiguous or needs_review. | Статус проверки: confident, ambiguous или needs_review. | R/call-schema.R:196 | — |
| call_table.desc.component_status | Component consistency: consistent or heterogeneous. | Согласованность компонентов: consistent или heterogeneous. | R/call-schema.R:197 | — |
| call_table.desc.mappability_status | Mappability flag: ok or low; low means the high-MAPQ read fraction is below 0.5. | Флаг картированности: ok или low; low означает, что доля прочтений с высоким MAPQ ниже 0,5. | R/call-schema.R:198 | — |
| call_table.desc.confidence | Legacy donor-effect confidence: high, moderate or low. | Устаревшая достоверность донорного эффекта: high, moderate или low. | R/call-schema.R:199 | — |
| call_table.desc.class_group | Group-level call class before recipient assignment. | Класс находки на уровне группы до размещения у реципиента. | R/call-schema.R:200 | — |
| call_table.desc.recipient_component | Copy-set label of the assigned recipient contig. | Метка набора копий назначенного контига реципиента. | R/call-schema.R:201 | — |
| call_table.desc.recipient_chr_id | Assembly contig identifier of the assigned recipient. | Идентификатор контига сборки назначенного реципиента. | R/call-schema.R:202 | — |
| call_table.desc.recipient_chr_name | Display name of the assigned recipient contig. | Отображаемое имя назначенного контига реципиента. | R/call-schema.R:203 | — |
| call_table.desc.recipient_group | Homology group of the assigned recipient contig; differs from chromosome_group for cross-group placements. | Группа гомологии назначенного контига реципиента; отличается от chromosome_group при межгрупповых размещениях. | R/call-schema.R:204 | — |
| call_table.desc.recipient_start_bp | Start of the corresponding interval on the assigned recipient contig, in base pairs. | Начало соответствующего интервала на назначенном контиге реципиента, в парах оснований. | R/call-schema.R:205 | — |
| call_table.desc.recipient_end_bp | End of that interval on the assigned recipient contig, in base pairs. | Конец этого интервала на назначенном контиге реципиента, в парах оснований. | R/call-schema.R:206 | — |
| call_table.desc.recipient_run_start_bp | Start of the contiguous recipient depletion run used for assignment, in base pairs. | Начало непрерывного участка делеции реципиента, использованного для размещения, в парах оснований. | R/call-schema.R:207 | — |
| call_table.desc.recipient_run_end_bp | End of that depletion run, in base pairs. | Конец этого участка делеции, в парах оснований. | R/call-schema.R:208 | — |
| call_table.desc.recipient_deletion_width_bp | Width of the recipient depletion run, in base pairs; blank when no run was found. | Ширина участка делеции реципиента, в парах оснований; пусто, если участок не найден. | R/call-schema.R:209 | — |
| call_table.desc.recipient_expected_breadth | Baseline breadth expected for the recipient contig. | Базовая ожидаемая широта покрытия для контига реципиента. | R/call-schema.R:210 | — |
| call_table.desc.recipient_observed_breadth | Mean breadth observed inside the call interval on the recipient contig. | Средняя широта покрытия, наблюдаемая внутри интервала находки на контиге реципиента. | R/call-schema.R:211 | — |
| call_table.desc.recipient_mapq_fraction | Mean high-MAPQ read fraction in the recipient interval. | Средняя доля прочтений с высоким MAPQ в интервале реципиента. | R/call-schema.R:212 | — |
| call_table.desc.recipient_mean_depletion | Mean corrected log-depth depletion in the recipient interval; positive values mean loss. | Средняя скорректированная лог-глубина делеции в интервале реципиента; положительные значения означают потерю. | R/call-schema.R:213 | — |
| call_table.desc.recipient_depletion_se | Standard error of the mean recipient depletion. | Стандартная ошибка средней делеции реципиента. | R/call-schema.R:214 | — |
| call_table.desc.recipient_depletion_lower | Lower confidence bound of the recipient depletion estimate. | Нижняя доверительная граница оценки делеции реципиента. | R/call-schema.R:215 | — |
| call_table.desc.recipient_depletion_z | Z score of the recipient depletion estimate. | Z-оценка оценки делеции реципиента. | R/call-schema.R:216 | — |
| call_table.desc.recipient_width_overlap | Reciprocal overlap between the donor segment and the recipient depletion run, from 0 to 1. | Взаимное перекрытие донорного сегмента и участка делеции реципиента, от 0 до 1. | R/call-schema.R:217 | — |
| call_table.desc.recipient_fraction_ratio | Donor width fraction divided by recipient run width fraction; values near 1 mean matched widths. | Доля ширины донора, делённая на долю ширины участка реципиента; значения около 1 — совпадающие ширины. | R/call-schema.R:218 | — |
| call_table.desc.copy_state_group | Copy state before recipient assignment. | Копийность до размещения у реципиента. | R/call-schema.R:219 | — |
| call_table.desc.recipient_weight | Weight of the chosen recipient among plausible candidates, from 0 to 1. | Вес выбранного реципиента среди правдоподобных кандидатов, от 0 до 1. | R/call-schema.R:220 | — |
| call_table.desc.recipient_margin | Score difference between the chosen recipient and the best alternative; larger values are more decisive. | Разница оценок между выбранным реципиентом и лучшей альтернативой; большие значения более решающие. | R/call-schema.R:221 | — |
| call_table.desc.recipient_assignment | Recipient placement outcome: assigned, ambiguous, cross_group or no_recipient_loss. | Итог размещения у реципиента: assigned, ambiguous, cross_group или no_recipient_loss. | R/call-schema.R:222 | — |
| call_table.desc.recipient_cross_group | TRUE when the assigned recipient contig is outside the call's homology group. | TRUE, когда назначенный контиг реципиента вне группы гомологии находки. | R/call-schema.R:223 | — |
| call_table.desc.recipient_reliable | TRUE when the recipient candidate passed the calibrated reliability screen. | TRUE, когда кандидат-реципиент прошёл откалиброванную проверку надёжности. | R/call-schema.R:224 | — |
| call_table.desc.recipient_strong | TRUE when the recipient candidate passed the stronger depletion screen. | TRUE, когда кандидат-реципиент прошёл более строгую проверку делеции. | R/call-schema.R:225 | — |
| call_table.desc.recipient_alternative | Display name of the next-best recipient contig, when available. | Отображаемое имя следующего по качеству контига реципиента, если доступно. | R/call-schema.R:226 | — |
| call_table.desc.recipient_alternative_score | Assignment score of the next-best recipient contig. | Оценка размещения следующего по качеству контига реципиента. | R/call-schema.R:227 | — |
| call_table.desc.donor_component | Copy-set label of the donor contig that contributed the call. | Метка набора копий донорного контига, давшего находку. | R/call-schema.R:228 | — |
| call_table.desc.donor_chr_id | Assembly contig identifier of the donor contig. | Идентификатор контига сборки донорного контига. | R/call-schema.R:229 | — |
| call_table.desc.donor_chr_name | Display name of the donor contig. | Отображаемое имя донорного контига. | R/call-schema.R:230 | — |
| call_table.desc.donor_start_bp | Donor segment start in base pairs. | Начало донорного сегмента в парах оснований. | R/call-schema.R:231 | — |
| call_table.desc.donor_end_bp | Donor segment end in base pairs. | Конец донорного сегмента в парах оснований. | R/call-schema.R:232 | — |

---

## Coverage note — completeness pass

**Added:** 19 `errors.*`/`sidebar.upload_help.*`/`log.*` keys + 64 `call_table.label.*` + 64 `call_table.desc.*` = **147 keys**.

**Deliberately out of scope (with one-line reasons — not gaps):**

- **Config-error page** (`R/app-bam.R:317` "No MaMaMIA configuration was supplied…", plus "Configuration error", "Suggested fix:", "The production registry is mounted read-only…"): this page is returned by `bam_launcher_ui()` when `!state$ok` and `bam_launcher_server()` returns an empty server function — no `language()` reactiveVal and no switcher exist on that path, so Russian mode is unreachable there.
- **`output$preflight_table`** (`R/app-bam.R:5147` "Upload a BAM or select the bundled example to see the preflight.", `:5154` "Preflight failed; see the calculation log."): dead output — no `tableOutput("preflight_table")` anywhere in the UI.
- **`output$assembly_summary`** (`R/app-bam.R:5171`, `:5178`, `:5187`, `:5203` "Not paired (%d)…"): dead output — no `uiOutput("assembly_summary")` anywhere; consequently the whole `bam_assembly_summary_ui()` block ("Assembly", "Donor contigs", "Recipient contigs", "Pairs", "Analysis windows", "a contig is too short", "N windows", "unavailable", `R/app-bam.R:1405-1470`) never renders.
- **JS `"Preparing download..."`** (`R/app-bam.R:3667`): transient fallback, always overridden by the server-set `data-busy-text` (`explorer_download_button`); never visible.
- **JS `"Download failed with ' + response.status`** (`R/app-bam.R:3723`): thrown error caught by `.catch` → `setError()`; message never displayed.
- **`bam_format_size` "unknown size"** (`R/app-bam.R:642`): internal degenerate fallback, not a normal UI path.
- **Cohort export card** (`R/cohort-export.R:554` "Run the analysis to enable the exports.", `:567` "Cohort"/"Sample summary", `:573-577` hover note, and `bam_cohort_table_header` column labels/descriptions incl. `:368` "Calls that pass calibration filters."): `R/cohort-export.R` is deferred in the spec §2 scope matrix (P3); these render in the Export tab but belong to the deferred cohort-export feature. Note: the in-app **empty state** and **card title** are genuinely user-visible chrome and should be picked up in the P3 follow-up.

**Sync note for the implementer (no edit made):** `R/i18n.R` already contains keys that are NOT in this copy file (`common.no`/`common.yes`, `errors.file_label*`/`errors.file_more`/`errors.prefix`/`errors.suggestion.*`, `results.class_*`/`copy_*`/`status_*`/`confidence_*`/`assignment_*`, `results.unplaced_note`, `sidebar.example_note`, `tooltip.confidence_*`/`copy_*`/`evidence_*`/`assigned`/`cross_group_candidate`/`group_only`/`total_gain`/`total_loss`/`margin`, `tooltip.donor`/`recipient`, `log.run_failed_stage`/`log.run_params_rejected`/`log.run_refused_stage`). These were added outside this document. **Do not regenerate `i18n.R` from the copy file wholesale** — merge, or those keys will be dropped. Also: this document's older `errors.short_contigs` row (4-arg) is superseded by the 2-arg `errors.short_contigs` + `errors.short_contigs_plural` split the engine already implemented (`R/app-bam.R:1362-1371`); keep the split.

**Biological wording corrected (report-only):** none of the newly-added strings required a factual correction — "breadth" → "широта покрытия", "component" → "компонент" (subgenome/copy-set), and the unplaced/deleted recipient chromosome-state rulings from the glossary were applied verbatim.
