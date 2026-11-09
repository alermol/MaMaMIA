## Dashboard theme, CSS and asset helpers. These live next to the
## explorer UI in dashboard.R; the launcher only wires them in.

## The bundled Latin subset keeps the interface independent of the host's
## installed fonts. Inter is distributed under the SIL Open Font License; the
## licence text ships next to the font.
#' @noRd
mamamia_font_dependency <- function() {
    path <- mamamia_path("www")
    if (!nzchar(path) || !dir.exists(path)) {
        return(NULL)
    }
    htmltools::htmlDependency(
        name = "mamamia-font",
        version = "1.0.0",
        src = c(file = path),
        stylesheet = "mamamia-font.css",
        all_files = TRUE
    )
}

#' @noRd
mamamia_app_version <- function() {
    tryCatch(mamamia_version(), error = function(e) "")
}

#' @noRd
mamamia_asset_url <- function(file) {
    dep <- mamamia_font_dependency()
    if (is.null(dep)) {
        return(file)
    }
    paste0(dep$name, "-", dep$version, "/", file)
}

#' @noRd
mamamia_favicon <- function() {
    shiny::tags$link(
        rel = "icon",
        type = "image/svg+xml",
        href = mamamia_asset_url("mamamia-logo.svg")
    )
}

## A warm, paper-like surface with a single deep-teal accent keeps the
## scientific content in front. The chart palette stays blue-red so scores and
## calls remain distinguishable from the interface accent.
#' @noRd
bam_dashboard_theme <- function() {
    bslib::bs_theme(
        version = 5,
        primary = "#0f5c53",
        secondary = "#6d7168",
        success = "#2f7d5d",
        info = "#2c6e8f",
        warning = "#b7791f",
        danger = "#b23a3a",
        bg = "#f6f5f0",
        fg = "#1b1e1a",
        base_font = paste(
            "'Inter'", "system-ui", "-apple-system", "'Segoe UI'", "Roboto",
            "'Helvetica Neue'", "Arial", "sans-serif",
            sep = ", "
        ),
        code_font = paste(
            "ui-monospace", "SFMono-Regular", "Menlo", "Consolas",
            "'Liberation Mono'", "monospace",
            sep = ", "
        ),
        "border-color" = "#e4e1d7",
        "border-radius" = ".5rem",
        "font-size-base" = ".95rem",
        "headings-font-weight" = "700"
    )
}

## The page is a three-column, three-row grid. The 2rem first column is a
## gutter; the header and footer span the sidebar and main columns. The grid
## keeps the log visible while the result tabs scroll independently.
#' @noRd
bam_dashboard_css <- function() {
    "
    :root {
        --mm-paper: #f6f5f0;
        --mm-surface: #ffffff;
        --mm-sidebar: #fbfaf6;
        --mm-ink: #1b1e1a;
        --mm-muted: #6d7168;
        --mm-line: #e4e1d7;
        --mm-line-strong: #d2cec0;
        --mm-accent: #0f5c53;
        --mm-accent-rgb: 15, 92, 83;
        --mm-accent-soft: rgba(15, 92, 83, .07);
        --mm-warm: #b7791f;
        --mm-shadow: 0 1px 2px rgba(27, 30, 26, .04),
            0 10px 30px rgba(27, 30, 26, .05);
    }
    .container-fluid {
        height: 100vh;
        padding: 0 !important;
    }
    html, body { height: 100%; }
    body {
        overflow: hidden;
        background-color: var(--mm-paper);
        background-image: radial-gradient(
            rgba(27, 30, 26, .028) .6px, transparent .6px
        );
        background-size: 18px 18px;
        color: var(--mm-ink);
        -webkit-font-smoothing: antialiased;
        text-rendering: optimizeLegibility;
    }

    .mamamia-header {
        position: relative;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 1rem;
        min-width: 0;
        height: 100%;
        padding: 0 .9rem;
        border-bottom: 1px solid var(--mm-line);
        background: var(--mm-surface);
    }
    .mamamia-header::after {
        content: '';
        position: absolute;
        right: 0;
        bottom: -1px;
        left: 0;
        height: 1px;
        background: linear-gradient(
            90deg,
            var(--mm-accent),
            rgba(183, 121, 31, .5) 32%,
            transparent 72%
        );
    }
    .mamamia-brand {
        display: flex;
        align-items: center;
        gap: .65rem;
        min-width: 0;
    }
    .mamamia-mark {
        display: grid;
        flex: 0 0 auto;
        place-items: center;
        width: 2.05rem;
        height: 2.05rem;
        border-radius: .62rem;
        background: var(--mm-ink);
        box-shadow: 0 2px 6px rgba(27, 30, 26, .18);
        color: var(--mm-paper);
    }
    .mamamia-mark svg {
        display: block;
        width: 1.3rem;
        height: 1.3rem;
    }
    .mamamia-brand-text { min-width: 0; }
    .mamamia-title-row {
        display: flex;
        align-items: baseline;
        gap: .45rem;
    }
    .mamamia-title {
        margin: 0;
        font-size: 1.55rem;
        font-weight: 800;
        line-height: 1;
        letter-spacing: -.035em;
    }
    .mamamia-version {
        padding: .05rem .35rem;
        border: 1px solid var(--mm-line);
        border-radius: 999px;
        color: var(--mm-muted);
        font-size: .58rem;
        font-weight: 650;
        letter-spacing: .04em;
        white-space: nowrap;
    }
    .mamamia-subtitle {
        margin: .2rem 0 0;
        color: var(--mm-muted);
        font-size: .7rem;
        letter-spacing: .012em;
    }
    .mamamia-header-note {
        display: flex;
        align-items: center;
        gap: .45rem;
        color: var(--mm-muted);
        font-size: .64rem;
        font-weight: 600;
        letter-spacing: .1em;
        text-transform: uppercase;
        white-space: nowrap;
    }
    .mamamia-header-note::before {
        content: '';
        width: .38rem;
        height: .38rem;
        border-radius: 50%;
        background: var(--mm-accent);
        box-shadow: 0 0 0 3px var(--mm-accent-soft);
    }
    @media (max-width: 980px) {
        .mamamia-header-note { display: none; }
    }

    .mamamia-sidebar {
        min-width: 0;
        min-height: 0;
        overflow-x: hidden;
        overflow-y: auto;
        padding: .75rem .85rem 1rem;
        border-right: 1px solid var(--mm-line);
        background: var(--mm-sidebar);
        font-size: .85rem;
        line-height: 1.35;
        scrollbar-width: thin;
        scrollbar-color: var(--mm-line-strong) transparent;
    }
    .mamamia-sidebar::-webkit-scrollbar { width: 7px; }
    .mamamia-sidebar::-webkit-scrollbar-thumb {
        border: 2px solid var(--mm-sidebar);
        border-radius: 999px;
        background: var(--mm-line-strong);
    }
    .mamamia-sidebar-section { margin-top: .15rem; }
    .mamamia-sidebar .mamamia-section-label {
        margin: .1rem 0 .45rem;
        color: var(--mm-muted);
        font-size: .6rem;
        font-weight: 750;
        letter-spacing: .13em;
        text-transform: uppercase;
    }
    .mamamia-sidebar .mamamia-sidebar-intro {
        margin-bottom: .55rem;
        color: var(--mm-muted);
        font-size: .7rem;
        line-height: 1.4;
    }
    .mamamia-sidebar .mamamia-sidebar-intro p {
        margin: 0;
    }
    .mamamia-sidebar .mamamia-recovery {
        margin: .1rem 0 .8rem;
        padding: .72rem .78rem .78rem;
        border: 1px solid rgba(183, 121, 31, .26);
        border-left: 3px solid var(--mm-warm);
        border-radius: .55rem;
        background: linear-gradient(
            180deg,
            rgba(183, 121, 31, .095),
            rgba(183, 121, 31, .035)
        );
        box-shadow: 0 1px 2px rgba(27, 30, 26, .03);
    }
    .mamamia-recovery-head {
        display: flex;
        align-items: center;
        gap: .42rem;
        margin-bottom: .35rem;
    }
    .mamamia-recovery-icon {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 1.25rem;
        height: 1.25rem;
        border-radius: 50%;
        background: rgba(183, 121, 31, .16);
        color: var(--mm-warm);
        font-size: .62rem;
    }
    .mamamia-recovery-title {
        color: var(--mm-ink);
        font-size: .74rem;
        font-weight: 700;
        line-height: 1.25;
    }
    .mamamia-recovery-message {
        margin: 0 0 .4rem;
        color: #675027;
        font-size: .67rem;
        line-height: 1.42;
        overflow-wrap: anywhere;
    }
    .mamamia-recovery-note {
        margin: 0 0 .45rem;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.4;
    }
    .mamamia-recovery-suggestion {
        margin: 0 0 .55rem;
        color: var(--mm-muted);
        font-size: .65rem;
        line-height: 1.4;
    }
    .mamamia-recovery-suggestion strong {
        color: #675027;
        font-weight: 700;
    }
    .mamamia-recovery-command {
        margin: 0 0 .55rem;
        padding: .42rem .5rem;
        border-radius: .35rem;
        background: rgba(255, 255, 255, .72);
        color: #5d4a25;
        font-size: .6rem;
        line-height: 1.35;
        white-space: pre-wrap;
        overflow-wrap: anywhere;
    }
    .mamamia-recovery-command code {
        color: inherit;
        font-size: inherit;
    }
    .mamamia-recovery-action {
        display: block;
        width: 100%;
        padding: .32rem .5rem;
        border: 1px solid rgba(183, 121, 31, .38);
        border-radius: .4rem;
        background: #fffaf0;
        color: #6b4f18;
        font-size: .68rem;
        font-weight: 650;
        line-height: 1.25;
        text-align: center;
        text-decoration: none;
        transition: background-color .15s ease, border-color .15s ease;
    }
    .mamamia-recovery-action:hover {
        border-color: rgba(183, 121, 31, .6);
        background: #fdf2dc;
        color: #563d10;
    }

    .mamamia-main {
        display: flex;
        flex-direction: column;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .6rem .7rem .65rem;
    }

    .mamamia-main > .tabbable,
    .mamamia-main > .card {
        flex: 1 1 auto;
        min-height: 0;
        margin-bottom: 0;
    }
    .mamamia-main > .card {
        overflow: hidden;
        border: 1px solid var(--mm-line);
        border-radius: .7rem;
        background: var(--mm-surface);
        box-shadow: var(--mm-shadow);
    }
    .mamamia-main .card-header {
        flex: 0 0 auto;
        padding: .3rem .6rem 0;
        border-bottom: 1px solid var(--mm-line);
        background: transparent;
    }
    .mamamia-main .card-body {
        min-height: 0;
        padding: 0;
        overflow: auto;
    }
    .mamamia-main .tab-content {
        min-height: 0;
        overflow: hidden;
    }
    /* bslib writes `padding:0` inline on fillable tab panes. Keep the tab
       content aligned with the tab labels instead of the card border. */
    .mamamia-main .tab-pane {
        height: 100%;
        padding: .65rem .7rem !important;
        overflow: auto;
    }
    /* DT is an html-fill-item: in the flex card body it can shrink to its
       control row while its table body overflows over following content. */
    .mamamia-main .card-body > .datatables.html-widget {
        flex: 0 0 auto !important;
        height: auto !important;
        min-height: 0;
    }
    .mamamia-main .nav-tabs {
        gap: .1rem;
        border-bottom: 0;
    }
    .mamamia-main .nav-link {
        padding: .48rem .72rem;
        border: 0;
        border-bottom: 2px solid transparent;
        border-radius: 0;
        color: var(--mm-muted);
        font-size: .8rem;
        font-weight: 600;
        line-height: 1.2;
    }
    .mamamia-main .nav-link:hover {
        border-bottom-color: var(--mm-line-strong);
        color: var(--mm-ink);
    }
    .mamamia-main .nav-link.active {
        border-bottom-color: var(--mm-accent);
        background: transparent;
        color: var(--mm-ink);
    }
    .mamamia-main .nav-link:focus-visible {
        outline: 2px solid var(--mm-accent);
        outline-offset: -2px;
    }

    .mamamia-footer {
        display: flex;
        flex-direction: column;
        gap: .2rem;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .25rem .85rem .35rem;
        border-top: 1px solid var(--mm-line);
        background: var(--mm-sidebar);
    }

    .mamamia-footer .shiny-text-output {
        flex: 1 1 auto;
        margin: 0;
        padding: .2rem .45rem;
        overflow-y: auto;
        border: 1px solid var(--mm-line);
        border-radius: .35rem;
        background: var(--mm-paper);
        color: #474a42;
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas,
            'Liberation Mono', monospace;
        font-size: .66rem;
        line-height: 1.35;
        white-space: pre-wrap;
    }

    .mamamia-footer .footer-status {
        display: flex;
        flex: 0 0 auto;
        align-items: center;
        gap: .5rem;
        min-height: 1.35rem;
        font-size: .72rem;
        line-height: 1.2;
    }

    .mamamia-footer .mamamia-footer-status {
        flex: 1 1 auto;
        min-width: 0;
    }
    .mamamia-footer .badge {
        padding: .32em .6em;
        border-radius: 999px;
        font-size: .58rem;
        font-weight: 750;
        letter-spacing: .07em;
        text-transform: uppercase;
        vertical-align: .08em;
    }

    .mamamia-footer .mamamia-footer-toggle {
        flex: 0 0 auto;
        padding: .15rem .45rem;
        border: 0;
        border-radius: .35rem;
        background: transparent;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.2;
    }
    .mamamia-footer .mamamia-footer-toggle:hover,
    .mamamia-footer .mamamia-footer-toggle:focus {
        background: var(--mm-paper);
        color: var(--mm-ink);
    }

    /* Collapsing the log gives its height back to the sidebar and main area. */
    #mamamia-grid.mamamia-footer-collapsed {
        grid-template-rows: 75px 1fr 2.25rem;
    }
    #mamamia-grid.mamamia-footer-collapsed .mamamia-footer .shiny-text-output {
        display: none;
    }
    @media (max-width: 500px) {
        #mamamia-grid.mamamia-footer-collapsed {
            grid-template-rows: 85px 350px 2.25rem 350px;
        }
    }

    .mamamia-sidebar .shiny-input-container {
        width: 100% !important;
        margin-bottom: .55rem;
    }
    .mamamia-sidebar .form-label,
    .mamamia-sidebar .control-label,
    .mamamia-sidebar .shiny-input-container > label {
        margin-bottom: .18rem;
        color: var(--mm-muted);
        font-size: .67rem;
        font-weight: 650;
        letter-spacing: .03em;
        line-height: 1.25;
    }
    .mamamia-sidebar .mamamia-profile-caption {
        margin: -.1rem 0 .45rem;
        color: var(--mm-muted);
        font-size: .67rem;
        line-height: 1.4;
    }
    .mamamia-sidebar .mamamia-param-heading {
        margin: .7rem 0 .35rem;
        color: var(--mm-muted);
        font-size: .6rem;
        font-weight: 750;
        letter-spacing: .12em;
        text-transform: uppercase;
    }
    .mamamia-sidebar .mamamia-param-heading:first-child {
        margin-top: .1rem;
    }
    .mamamia-sidebar .mamamia-param-actions {
        margin: .1rem 0 .55rem;
    }
    .mamamia-sidebar .mamamia-param-actions .btn {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: .35rem;
        width: 100%;
        height: 1.9rem;
        font-size: .72rem;
    }
    .mamamia-sidebar .mamamia-param-help {
        display: inline-flex;
        align-items: center;
        margin-left: .28rem;
        color: var(--mm-muted);
        font-size: .66rem;
        line-height: 1;
        cursor: help;
        vertical-align: baseline;
    }
    .mamamia-sidebar .mamamia-param-help:hover,
    .mamamia-sidebar .mamamia-param-help:focus {
        color: var(--mm-accent);
    }
    .mamamia-sidebar .mamamia-param-control {
        margin-bottom: .5rem;
    }
    .mamamia-sidebar .mamamia-param-control > .shiny-html-output {
        display: block;
        margin-bottom: .12rem;
    }
    .mamamia-sidebar .mamamia-param-control .shiny-label-null {
        display: none;
    }
    .mamamia-sidebar .mamamia-param-control > .form-group {
        margin-bottom: 0;
    }
    .mamamia-sidebar .control-label.mamamia-param-default {
        display: inline-block;
        color: var(--mm-ink);
        font-weight: 750;
    }
    .mamamia-sidebar .control-label.mamamia-param-changed {
        display: inline-block;
        padding: 0 .2rem;
        border-radius: .25rem;
        background: #f4e4c6;
        color: #7a4f0b;
        font-weight: 750;
    }
    .mamamia-sidebar .mamamia-param-control:has(.mamamia-param-changed) .form-control,
    .mamamia-sidebar .mamamia-param-control:has(.mamamia-param-changed) .form-select,
    .mamamia-sidebar .mamamia-param-control:has(.mamamia-param-changed) .shiny-input-select,
    .mamamia-sidebar .mamamia-resolution-slider:has(.mamamia-param-changed) .form-control {
        border-color: #d9a441;
        background-color: #fffaf0;
    }
    .mamamia-sidebar .mamamia-param-floors {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: .1rem .5rem;
        margin-bottom: .35rem;
    }
    .mamamia-sidebar .mamamia-param-floor {
        min-width: 0;
        margin-bottom: .1rem;
    }
    .mamamia-sidebar .mamamia-param-floor .control-label {
        letter-spacing: 0;
    }
    .mamamia-sidebar .mamamia-param-legend {
        margin: .55rem 0 0;
        color: var(--mm-muted);
        font-size: .62rem;
        line-height: 1.38;
    }
    .mamamia-sidebar .mamamia-param-legend .mamamia-param-changed {
        padding: 0 .16rem;
        border-radius: .2rem;
        background: #f4e4c6;
        color: #7a4f0b;
        font-weight: 750;
    }
    .mamamia-sidebar .form-control,
    .mamamia-sidebar .form-select {
        min-height: 2rem;
        height: auto;
        padding: .28rem .5rem;
        border-color: var(--mm-line);
        border-radius: .42rem;
        background-color: var(--mm-surface);
        font-size: .79rem;
        line-height: 1.3;
    }
    .mamamia-sidebar .form-control:focus,
    .mamamia-sidebar .form-select:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-sidebar .btn {
        padding: .25rem .5rem;
        font-size: .78rem;
        line-height: 1.25;
    }
    .mamamia-sidebar .mamamia-ghost-btn {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: .35rem;
        width: 100%;
        height: 1.95rem;
        margin-top: .45rem;
        border: 1px solid var(--mm-line);
        border-radius: .42rem;
        background: transparent;
        color: var(--mm-ink);
        font-size: .73rem;
        font-weight: 600;
    }
    .mamamia-sidebar .mamamia-ghost-btn:hover,
    .mamamia-sidebar .mamamia-ghost-btn:focus-visible {
        border-color: var(--mm-line-strong);
        background: var(--mm-paper);
        color: var(--mm-ink);
    }

    .mamamia-upload .form-group { margin-bottom: 0; }
    .mamamia-upload .control-label {
        position: absolute !important;
        width: 1px !important;
        height: 1px !important;
        padding: 0 !important;
        margin: -1px !important;
        overflow: hidden !important;
        clip: rect(0, 0, 0, 0) !important;
        white-space: nowrap !important;
        border: 0 !important;
    }
    .mamamia-upload-head {
        display: flex;
        align-items: center;
        gap: .5rem;
        margin-bottom: .45rem;
        cursor: pointer;
    }
    .mamamia-upload-copy {
        flex: 1 1 auto;
        min-width: 0;
    }
    .mamamia-upload-icon {
        display: grid;
        flex: 0 0 auto;
        place-items: center;
        width: 1.65rem;
        height: 1.65rem;
        border-radius: .45rem;
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
        font-size: .74rem;
    }
    .mamamia-upload-title {
        display: block;
        font-size: .76rem;
        font-weight: 650;
        line-height: 1.2;
    }
    .mamamia-upload-hint {
        display: block;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.3;
    }
    .mamamia-upload-title-row {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: .35rem;
        min-width: 0;
    }
    .mamamia-upload-limit {
        display: inline-flex;
        flex: 0 0 auto;
        align-items: center;
        gap: .12rem;
        color: var(--mm-muted);
        font-size: .6rem;
        font-variant-numeric: tabular-nums;
        white-space: nowrap;
    }
    .mamamia-upload-help {
        width: auto;
        height: auto;
        margin: 0;
        padding: 0;
        border: 0;
        border-radius: 50%;
        background: transparent;
        color: var(--mm-muted);
        font-size: .72rem;
        line-height: 1;
        box-shadow: none;
        opacity: .72;
        vertical-align: middle;
    }
    .mamamia-upload-help:hover,
    .mamamia-upload-help:focus {
        border: 0;
        background: transparent;
        color: var(--mm-accent);
        box-shadow: none;
        opacity: 1;
    }
    .mamamia-upload .input-group {
        gap: .4rem;
        padding: .5rem;
        border: 1px dashed var(--mm-line-strong);
        border-radius: .55rem;
        background: var(--mm-surface);
        cursor: pointer;
        transition: border-color .15s ease, background-color .15s ease;
    }
    .mamamia-upload .input-group:hover,
    .mamamia-upload.is-dragover .input-group {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
    }
    .mamamia-upload .input-group-btn {
        display: block;
        width: 100%;
        margin: 0;
    }
    .mamamia-upload .btn-file {
        display: flex;
        align-items: center;
        justify-content: center;
        width: 100%;
        height: 1.85rem;
        padding: 0 .55rem;
        border: 0;
        border-radius: .38rem;
        background: var(--mm-ink);
        color: var(--mm-paper);
        font-size: .7rem;
        font-weight: 650;
        letter-spacing: .01em;
        box-shadow: none;
    }
    .mamamia-upload .btn-file:hover,
    .mamamia-upload .btn-file:focus {
        background: var(--mm-accent);
        color: #fff;
    }
    .mamamia-upload input.form-control[readonly] {
        flex: 1 1 100%;
        min-height: 0;
        height: auto;
        padding: 0 .1rem;
        border: 0;
        border-radius: 0;
        background: transparent;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.3;
        text-align: center;
        box-shadow: none;
    }
    .mamamia-upload .progress {
        height: 3px;
        margin: .35rem 0 0;
        border-radius: 999px;
        background: var(--mm-line);
    }
    /* Shiny writes the file name (and then 'Upload complete') into the
       progress bar. Keep the bar as a pure fill indicator. */
    .mamamia-upload .progress-bar {
        overflow: hidden;
        background-color: var(--mm-accent);
        color: transparent;
        font-size: 0;
        line-height: 0;
        text-indent: -9999px;
        white-space: nowrap;
    }
    .mamamia-upload-timer {
        display: none;
        align-items: baseline;
        gap: .4rem;
        margin-top: .35rem;
        color: var(--mm-muted);
        font-size: .64rem;
        letter-spacing: .04em;
    }
    .mamamia-upload-timer.is-visible {
        display: flex;
    }
    .mamamia-upload-timer.is-running .mamamia-upload-timer-value {
        color: var(--mm-accent);
    }
    .mamamia-upload-timer.is-error .mamamia-upload-timer-value {
        color: var(--mm-warm);
    }
    .mamamia-upload-timer-label {
        font-weight: 600;
        text-transform: uppercase;
    }
    .mamamia-upload-timer-value {
        min-width: 4.2em;
        color: var(--mm-ink);
        font-weight: 650;
        font-variant-numeric: tabular-nums;
        letter-spacing: .02em;
    }

    .mamamia-sidebar .accordion {
        --bs-accordion-bg: transparent;
        --bs-accordion-btn-bg: transparent;
        --bs-accordion-active-bg: transparent;
        --bs-accordion-border-color: var(--mm-line);
        --bs-accordion-btn-focus-box-shadow: none;
        margin-top: .6rem;
    }
    .mamamia-sidebar .accordion-item {
        border: 0;
        border-top: 1px solid var(--mm-line);
        background: transparent;
    }
    .mamamia-sidebar .accordion-button {
        padding: .62rem .1rem;
        background: transparent;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 750;
        letter-spacing: .12em;
        line-height: 1.2;
        text-transform: uppercase;
        box-shadow: none;
    }
    .mamamia-sidebar .accordion-button:not(.collapsed) {
        color: var(--mm-ink);
    }
    .mamamia-sidebar .accordion-button::after {
        width: .62rem;
        height: .62rem;
        background-size: .62rem;
        opacity: .55;
    }
    .mamamia-sidebar .accordion-body {
        padding: .05rem .1rem .8rem;
    }
    .mamamia-sidebar .shiny-options-group {
        display: flex;
        flex-direction: column;
        gap: .2rem;
    }
    .mamamia-sidebar .radio-inline {
        display: flex;
        align-items: center;
        gap: .4rem;
        margin: 0;
        font-size: .76rem;
    }
    .mamamia-sidebar .radio-inline input,
    .mamamia-sidebar .form-check-input {
        accent-color: var(--mm-accent);
    }
    .mamamia-sidebar .form-check-label,
    .mamamia-sidebar .radio label,
    .mamamia-sidebar .checkbox label {
        font-size: .76rem;
        font-weight: 450;
        line-height: 1.35;
    }
    .mamamia-sidebar .text-muted.small,
    .mamamia-sidebar small {
        color: var(--mm-muted) !important;
        font-size: .69rem;
        line-height: 1.4;
    }
    .mamamia-sidebar .selectize-input {
        min-height: 2rem;
        padding: .25rem .4rem;
        border-color: var(--mm-line);
        border-radius: .42rem;
        background: var(--mm-surface);
        font-size: .74rem;
    }
    .mamamia-sidebar .selectize-control.multi .selectize-input > div {
        border: 0;
        border-radius: .3rem;
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
        font-size: .68rem;
    }
    .mamamia-run-panel {
        margin-top: .35rem;
        padding-top: .7rem;
        border-top: 1px solid var(--mm-line);
    }
    .mamamia-run-panel .btn {
        width: 100%;
    }
    .mamamia-run-panel #run {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: .45rem;
        height: 2.3rem;
        border: 0;
        border-radius: .45rem;
        background: var(--mm-accent);
        color: #ffffff;
        font-size: .76rem;
        font-weight: 650;
        letter-spacing: .015em;
        box-shadow: 0 2px 5px rgba(15, 92, 83, .22);
        transition: background-color .15s ease, transform .1s ease;
    }
    .mamamia-run-panel #run:hover:not(.disabled) {
        background: #0b4a43;
    }
    .mamamia-run-panel #run:active:not(.disabled) {
        transform: translateY(1px);
    }
    .mamamia-run-panel #run.disabled {
        box-shadow: none;
        opacity: .42;
    }

    .mamamia-main > .mamamia-filter-bar {
        flex: 0 0 auto;
    }
    .mamamia-filter-bar {
        display: flex;
        align-items: flex-end;
        gap: 0;
        padding: .1rem 0 .55rem;
        border: 0;
        border-bottom: 1px solid var(--mm-line);
        border-radius: 0;
        background: transparent;
    }
    .mamamia-sample-picker {
        flex: 0 0 auto;
    }
    /* Never hide the empty wrapper: Shiny suspends display:none outputs, so
       the picker could not render once a multi-sample BAM arrives. The
       populated input carries its own spacing so the empty state is flush. */
    .mamamia-sample-picker .shiny-input-container,
    .mamamia-sample-picker .form-group {
        width: 10.5rem;
        margin-bottom: 0;
        margin-right: .55rem;
    }
    .mamamia-sample-picker .control-label {
        margin-bottom: .15rem;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 700;
        letter-spacing: .1em;
        text-transform: uppercase;
    }
    .mamamia-sample-picker .form-select {
        min-height: 1.9rem;
        padding: .2rem .6rem;
        border-color: var(--mm-line);
        border-radius: .45rem;
        background-color: var(--mm-surface);
        font-size: .72rem;
    }
    .mamamia-filter-controls {
        display: flex;
        flex: 1 1 auto;
        flex-wrap: wrap;
        align-items: center;
        gap: .4rem;
        min-width: 0;
    }
    .mamamia-filter-pills {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: .45rem;
        width: 100%;
    }
    .mamamia-filter-pills .shiny-input-container {
        width: auto;
        min-width: 0;
        margin-bottom: 0;
        font-size: 0;
    }
    .mamamia-filter-pills .bootstrap-select {
        flex: 0 0 auto;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle {
        min-height: 1.9rem;
        height: 1.9rem;
        padding: 0 .78rem;
        border: 1px solid var(--mm-line);
        border-radius: .45rem;
        background: var(--mm-surface);
        color: var(--mm-muted);
        font-size: .72rem;
        font-weight: 600;
        box-shadow: none;
        transition: border-color .15s ease, color .15s ease,
            background-color .15s ease;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle .filter-option {
        display: flex;
        align-items: center;
        height: 100%;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle::after {
        width: 5px;
        height: 5px;
        margin: 0 0 0 .42rem;
        border: 0;
        border-right: 1.5px solid currentColor;
        border-bottom: 1.5px solid currentColor;
        transform: rotate(45deg) translateY(-1px);
        transition: transform .15s ease;
    }
    .mamamia-filter-pills .bootstrap-select.show > .dropdown-toggle::after {
        transform: rotate(225deg) translateY(-1px);
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:hover {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:focus {
        outline: none !important;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:focus-visible {
        outline: 2px solid var(--mm-accent) !important;
        outline-offset: 2px !important;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:not(.bs-placeholder) {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .bootstrap-select.show > .dropdown-toggle {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-filter-pills .bootstrap-select .filter-option-inner-inner {
        font-size: .72rem;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-menu {
        width: 10rem;
        min-width: 10rem;
        max-width: 10rem;
        padding: .28rem 0;
        border-color: var(--mm-line);
        border-radius: .45rem;
        box-shadow: 0 .45rem 1.2rem rgba(27, 30, 26, .13);
        font-size: .72rem;
    }
    .mamamia-filter-pills .bootstrap-select .dropdown-item {
        padding: .3rem .8rem;
    }
    .mamamia-filter-pills .bootstrap-select .dropdown-item.active,
    .mamamia-filter-pills .bootstrap-select .dropdown-item:active {
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .mamamia-filter-switch {
        display: flex;
        align-items: center;
        min-height: 1.9rem;
        margin-left: auto;
    }
    .mamamia-filter-pills .form-check.form-switch {
        display: flex;
        align-items: center;
        gap: .42rem;
        min-height: 1.9rem;
        margin: 0;
        padding-left: 2.4em;
        font-size: .72rem;
        white-space: nowrap;
    }
    .mamamia-filter-pills .form-check-input {
        margin-top: 0;
        border-color: var(--mm-line-strong);
        background-color: var(--mm-surface);
    }
    .mamamia-filter-pills .form-check-input:checked {
        border-color: var(--mm-accent);
        background-color: var(--mm-accent);
    }
    .mamamia-filter-pills .form-check-input:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-filter-pills .form-check-label {
        color: var(--mm-muted);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-segment-controls {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-end;
        gap: .5rem;
        margin-bottom: .55rem;
        padding-bottom: .55rem;
        border-bottom: 1px solid var(--mm-line);
    }
    .mamamia-segment-controls .shiny-input-container {
        margin-bottom: 0;
    }
    .mamamia-segment-controls .form-group {
        margin-bottom: 0;
    }
    .mamamia-segment-controls .control-label {
        margin-bottom: .18rem;
        color: var(--mm-muted);
        font-size: .64rem;
        font-weight: 650;
        letter-spacing: .03em;
    }
    .mamamia-segment-controls .form-control {
        height: 2rem;
        padding: .2rem .5rem;
        border-color: var(--mm-line);
        border-radius: .4rem;
        font-size: .76rem;
    }
    .mamamia-segment-controls .mamamia-call-columns {
        flex: 0 1 11rem;
        min-width: 9rem;
    }
    .mamamia-call-columns .bootstrap-select {
        display: block;
        width: 100% !important;
        padding: 0;
        border: 0;
        background: transparent;
        box-shadow: none;
    }
    .mamamia-call-columns .dropdown-toggle {
        width: 100%;
        text-align: left;
    }
    .mamamia-call-columns .bootstrap-select .dropdown-toggle .filter-option {
        display: flex;
        align-items: center;
        height: 100%;
    }
    .mamamia-call-columns .bootstrap-select .dropdown-toggle .filter-option-inner {
        display: block;
        width: 100%;
    }
    .mamamia-segment-controls .btn {
        height: 2rem;
        padding: 0 .65rem;
        border: 1px solid var(--mm-line);
        border-radius: .4rem;
        background: var(--mm-surface);
        color: var(--mm-ink);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-segment-controls .btn:hover {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-segment-controls .d-flex { gap: .4rem !important; }
    .mamamia-export-panel {
        display: flex;
        flex-direction: column;
        gap: .9rem;
        max-width: none;
    }
    .mamamia-export-actions {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(17rem, 1fr));
        gap: .8rem;
    }
    @media (min-width: 56rem) and (max-width: 74.99rem) {
        .mamamia-export-actions {
            grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        .mamamia-export-actions .mamamia-export-card:last-child {
            display: grid;
            grid-column: 1 / -1;
            grid-template-columns: minmax(0, 1fr) auto;
            align-content: center;
            align-items: center;
            column-gap: 1rem;
            row-gap: .1rem;
            min-height: 0;
            padding: .6rem .7rem;
        }
        .mamamia-export-actions .mamamia-export-card:last-child h5 {
            grid-column: 1;
            margin-bottom: 0;
        }
        .mamamia-export-actions .mamamia-export-card:last-child .shiny-text-output {
            grid-column: 1;
            margin: 0;
        }
        .mamamia-export-actions .mamamia-export-card:last-child .mamamia-download-row {
            grid-column: 2;
            grid-row: 1 / span 2;
            margin: 0;
            align-self: center;
        }
    }
    .mamamia-export-card {
        display: flex;
        flex-direction: column;
        min-width: 0;
        padding: .85rem .7rem;
        border: 1px solid var(--mm-line);
        border-radius: .65rem;
        background: var(--mm-surface);
    }
    .mamamia-export-actions .mamamia-export-card {
        min-height: 8.5rem;
    }
    .mamamia-export-card .mamamia-export-lead {
        max-width: 55ch;
    }
    .mamamia-export-actions .mamamia-export-card .btn {
        margin-top: auto;
        align-self: flex-start;
    }
    .mamamia-download-row {
        display: flex;
        align-items: center;
        gap: .55rem;
        flex-wrap: wrap;
    }
    .mamamia-export-actions .mamamia-export-card .mamamia-download-row {
        margin-top: auto;
    }
    .mamamia-export-actions .mamamia-export-card .mamamia-download-row .btn {
        margin-top: 0;
        align-self: center;
    }
    .mamamia-download-row.is-busy .btn {
        opacity: .7;
        cursor: progress;
        pointer-events: none;
    }
    .mamamia-download-status {
        display: none;
        align-items: center;
        gap: .42rem;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.25;
    }
    .mamamia-download-status.is-visible {
        display: inline-flex;
    }
    .mamamia-download-status.is-error {
        color: #a13a2f;
    }
    .mamamia-download-spinner {
        width: .8rem;
        height: .8rem;
        flex: 0 0 auto;
        border: .12rem solid var(--mm-line);
        border-top-color: var(--mm-accent);
        border-radius: 50%;
        animation: mamamia-download-spin .8s linear infinite;
    }
    .mamamia-download-status.is-error .mamamia-download-spinner {
        display: none;
    }
    @keyframes mamamia-download-spin {
        to { transform: rotate(360deg); }
    }
    @media (prefers-reduced-motion: reduce) {
        .mamamia-download-spinner {
            animation-duration: 2.4s;
        }
    }
    .mamamia-export-card-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: .8rem;
        flex-wrap: wrap;
    }
    .mamamia-export-card-head h5 {
        margin: 0;
    }
    .mamamia-export-card-foot {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: .8rem;
        flex-wrap: wrap;
        margin-top: .7rem;
    }
    .mamamia-export-note {
        margin: 0;
        max-width: none;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.4;
    }
    .mamamia-export-actions .mamamia-export-section + .mamamia-export-section {
        margin-top: 0;
    }
    .mamamia-export-panel #cohort_ui {
        margin-top: 0;
    }
    .mamamia-export-lead {
        margin: 0 0 .55rem;
        color: var(--mm-muted);
        font-size: .7rem;
        line-height: 1.42;
    }
    .mamamia-call-note {
        margin: .2rem 0 .5rem;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.35;
    }
    .mamamia-export-panel h5 {
        margin: 1rem 0 .4rem;
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 750;
        letter-spacing: .12em;
        text-transform: uppercase;
    }
    .mamamia-export-card h5 {
        margin: 0 0 .35rem;
    }
    .mamamia-export-panel .mamamia-export-card-head h5 {
        margin: 0;
    }
    .mamamia-export-panel .mamamia-export-section h5 {
        margin-top: 0;
    }
    .mamamia-export-panel .mamamia-call-columns {
        max-width: 16rem;
        margin-bottom: .55rem;
    }
    .mamamia-export-panel .mamamia-export-card-head .mamamia-call-columns {
        max-width: none;
        margin-bottom: 0;
    }
    .mamamia-export-panel .mamamia-export-card-head .mamamia-call-columns .form-group {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: flex-end;
        gap: .55rem;
        margin: 0;
    }
    .mamamia-export-panel .mamamia-export-card-head .mamamia-call-columns .control-label {
        margin: 0;
        white-space: nowrap;
    }
    .mamamia-export-panel .mamamia-export-card-head .mamamia-call-columns .bootstrap-select {
        width: auto !important;
        min-width: 9rem;
    }
    .mamamia-export-panel .mamamia-export-card-head .mamamia-call-columns .dropdown-toggle .filter-option {
        justify-content: center;
        text-align: center;
    }
    .mamamia-export-panel .mamamia-call-columns .control-label {
        margin-bottom: .18rem;
        color: var(--mm-muted);
        font-size: .64rem;
        font-weight: 650;
        letter-spacing: .03em;
    }
    .mamamia-export-panel .mamamia-call-columns .dropdown-toggle {
        height: 2rem;
        padding: 0 .65rem;
        border: 1px solid var(--mm-line);
        border-radius: .4rem;
        background: var(--mm-surface);
        color: var(--mm-ink);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-export-panel .mamamia-call-columns .dropdown-toggle:hover,
    .mamamia-export-panel .mamamia-call-columns .dropdown-toggle:focus {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-export-panel h5:first-child { margin-top: 0; }
    .mamamia-export-panel .btn {
        border-radius: .4rem;
        font-size: .78rem;
    }
    .mamamia-export-panel .btn.disabled,
    .mamamia-export-panel .btn:disabled {
        opacity: .45;
        pointer-events: none;
    }
    .mamamia-export-panel #explorer-selection_status {
        margin-bottom: .35rem;
        color: var(--mm-muted);
        font-size: .74rem;
    }
    .mamamia-export-panel .table-responsive {
        margin-bottom: .7rem;
        scrollbar-width: thin;
        scrollbar-color: var(--mm-line-strong) transparent;
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar {
        height: 7px;
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar-thumb {
        border-radius: 999px;
        background: var(--mm-line-strong);
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar-track {
        background: transparent;
    }
    .mamamia-export-panel .table {
        margin-bottom: 0;
        font-size: .72rem;
        white-space: nowrap;
    }
    .mamamia-export-panel .table th {
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 700;
        letter-spacing: .04em;
    }
    .mamamia-export-panel .table td,
    .mamamia-export-panel .table th {
        padding: .35rem .5rem;
    }
    .mamamia-karyotype {
        display: flex;
        flex-direction: column;
        height: 100%;
        min-height: 0;
    }
    .mamamia-karyotype > div:first-child {
        flex: 0 0 auto;
        justify-content: space-between;
        gap: .5rem !important;
        padding: 0 0 .5rem !important;
    }
    .mamamia-karyotype-plot {
        position: relative;
        display: flex;
        flex: 1 1 auto;
        flex-direction: column;
        min-height: 0;
    }
    .mamamia-karyotype-plot > .html-widget {
        flex: 1 1 auto;
        min-height: 0;
        height: 100% !important;
    }
    .mamamia-analysis-loading {
        position: absolute;
        inset: 0;
        z-index: 4;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        gap: .85rem;
        padding: 1.5rem;
        text-align: center;
        background: var(--mm-surface);
        opacity: 0;
        visibility: hidden;
        transition: opacity .25s ease, visibility .25s ease;
        pointer-events: none;
    }
    .mamamia-analysis-loading.is-running {
        opacity: 1;
        visibility: visible;
    }
    .mamamia-analysis-loader {
        position: relative;
        width: 92px;
        height: 92px;
    }
    .mamamia-analysis-loader::before {
        content: '';
        position: absolute;
        inset: 12%;
        border-radius: 50%;
        background: var(--mm-accent-soft);
        animation: mamamia-loader-halo 2.4s ease-in-out infinite;
    }
    .mamamia-analysis-loader svg {
        position: relative;
        display: block;
        width: 100%;
        height: 100%;
        overflow: visible;
    }
    .mamamia-analysis-orbit {
        fill: none;
        stroke: var(--mm-line-strong);
        stroke-width: 1.4;
        stroke-dasharray: 4 7;
        transform-box: view-box;
        transform-origin: 40px 40px;
        animation: mamamia-loader-orbit 10s linear infinite;
    }
    .mamamia-analysis-strand {
        fill: none;
        stroke: var(--mm-accent);
        stroke-width: 2.6;
        stroke-linecap: round;
        animation: mamamia-loader-strand 2.4s ease-in-out infinite;
    }
    .mamamia-analysis-strand-warm {
        stroke: var(--mm-warm);
        animation-delay: .3s;
    }
    .mamamia-analysis-chiasma {
        fill: var(--mm-accent);
        transform-box: view-box;
        transform-origin: 40px 40px;
        animation: mamamia-loader-chiasma 1.6s ease-in-out infinite;
    }
    .mamamia-analysis-eyebrow {
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 700;
        letter-spacing: .14em;
        text-transform: uppercase;
    }
    .mamamia-analysis-stage {
        color: var(--mm-ink);
        font-size: 1rem;
        font-weight: 650;
        line-height: 1.3;
    }
    .mamamia-analysis-timer {
        display: flex;
        align-items: baseline;
        gap: .4rem;
        color: var(--mm-muted);
        font-size: .72rem;
        letter-spacing: .04em;
    }
    .mamamia-analysis-timer-label {
        font-weight: 600;
        text-transform: uppercase;
    }
    .mamamia-analysis-timer-value {
        min-width: 4.2em;
        color: var(--mm-ink);
        font-weight: 650;
        font-variant-numeric: tabular-nums;
        letter-spacing: .02em;
        text-align: left;
    }
    .mamamia-analysis-track {
        width: min(240px, 72%);
        height: 3px;
        overflow: hidden;
        border-radius: 999px;
        background: var(--mm-line);
    }
    .mamamia-analysis-fill {
        display: block;
        width: 0;
        height: 100%;
        border-radius: inherit;
        background: linear-gradient(90deg, var(--mm-accent), var(--mm-warm));
        transition: width .25s ease;
    }
    @keyframes mamamia-loader-orbit {
        to { transform: rotate(360deg); }
    }
    @keyframes mamamia-loader-strand {
        0%, 100% { opacity: .55; }
        50% { opacity: 1; }
    }
    @keyframes mamamia-loader-chiasma {
        0%, 100% { transform: scale(.75); opacity: .55; }
        50% { transform: scale(1.25); opacity: 1; }
    }
    @keyframes mamamia-loader-halo {
        0%, 100% { transform: scale(.85); opacity: .45; }
        50% { transform: scale(1.1); opacity: .9; }
    }
    @media (prefers-reduced-motion: reduce) {
        .mamamia-analysis-orbit,
        .mamamia-analysis-strand,
        .mamamia-analysis-chiasma,
        .mamamia-analysis-loader::before {
            animation: none;
        }
        .mamamia-analysis-fill { transition: none; }
    }
    .mamamia-karyotype .shiny-input-container {
        width: auto;
        margin-bottom: 0;
    }
    .mamamia-karyotype .control-label {
        position: absolute !important;
        width: 1px !important;
        height: 1px !important;
        padding: 0 !important;
        margin: -1px !important;
        overflow: hidden !important;
        clip: rect(0, 0, 0, 0) !important;
        white-space: nowrap !important;
        border: 0 !important;
    }
    .mamamia-karyotype .shiny-options-group {
        display: inline-flex;
        gap: 2px;
        padding: 2px;
        border: 1px solid var(--mm-line);
        border-radius: .45rem;
        background: var(--mm-paper);
        font-size: 0;
    }
    .mamamia-karyotype .radio-inline {
        margin: 0;
        padding: 0;
        font-size: 0;
    }
    .mamamia-karyotype .radio-inline input {
        position: absolute;
        width: 1px;
        height: 1px;
        margin: 0;
        opacity: 0;
        pointer-events: none;
    }
    .mamamia-karyotype .radio-inline span {
        display: block;
        padding: .3rem .8rem;
        border-radius: .32rem;
        color: var(--mm-muted);
        font-size: .7rem;
        font-weight: 600;
        line-height: 1.25;
        text-align: center;
        cursor: pointer;
        transition: background-color .15s ease, color .15s ease;
    }
    .mamamia-karyotype .radio-inline input:checked + span {
        background: var(--mm-surface);
        box-shadow: 0 1px 2px rgba(27, 30, 26, .12);
        color: var(--mm-ink);
    }
    .mamamia-karyotype .radio-inline input:focus-visible + span {
        outline: 2px solid var(--mm-accent);
        outline-offset: 1px;
    }
    .mamamia-karyotype-plot-controls {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: .75rem;
    }
    .mamamia-karyotype-deletions {
        display: flex;
        align-items: center;
        gap: .45rem;
    }
    .mamamia-karyotype-deletions .form-check.form-switch {
        display: flex;
        align-items: center;
        gap: .42rem;
        min-height: 1.9rem;
        margin: 0;
        padding-left: 2.4em;
        font-size: .72rem;
        white-space: nowrap;
    }
    .mamamia-karyotype-deletions .form-check-input {
        margin-top: 0;
        border-color: var(--mm-line-strong);
        background-color: var(--mm-surface);
    }
    .mamamia-karyotype-deletions .form-check-input:checked {
        border-color: var(--mm-accent);
        background-color: var(--mm-accent);
    }
    .mamamia-karyotype-deletions .form-check-input:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-karyotype-deletions .form-check-label {
        color: var(--mm-muted);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-karyotype-deletions .shiny-text-output {
        color: var(--mm-muted);
        font-size: .7rem;
        font-weight: 600;
        white-space: nowrap;
    }
    .mamamia-karyotype .html-widget {
        flex: 1 1 auto;
        min-height: 0;
        height: 100% !important;
    }
    .mamamia-tooltip-modal.show {
        display: flex !important;
        align-items: center;
        justify-content: center;
        height: 100vh;
        height: 100dvh;
        padding: 1rem;
        overflow: hidden;
    }
    .mamamia-tooltip-modal .modal-dialog {
        width: min(800px, calc(100vw - 2rem));
        max-width: calc(100vw - 2rem);
        max-height: calc(100vh - 2rem);
        margin: 0;
    }
    .mamamia-tooltip-modal .modal-content {
        display: flex;
        flex-direction: column;
        max-height: calc(100vh - 2rem);
    }
    .mamamia-tooltip-modal .modal-header,
    .mamamia-tooltip-modal .modal-footer {
        flex: 0 0 auto;
    }
    .mamamia-tooltip-modal .modal-body {
        min-height: 0;
        overflow-y: auto;
        overscroll-behavior: contain;
    }
    .mamamia-tooltip-help {
        padding-right: .35rem;
    }
    .mamamia-help-lead {
        margin: 0 0 .55rem;
        color: var(--mm-ink);
        font-size: .82rem;
        line-height: 1.5;
    }
    .mamamia-help-callout {
        margin: 0 0 1rem;
        padding: .55rem .7rem;
        border-left: 3px solid var(--mm-accent);
        border-radius: .35rem;
        background: var(--mm-accent-soft);
        color: var(--mm-ink);
        font-size: .78rem;
        line-height: 1.45;
    }
    .mamamia-help-section + .mamamia-help-section {
        margin-top: 1rem;
    }
    .mamamia-help-intro {
        margin: 0 0 .4rem;
        color: var(--mm-muted);
        font-size: .78rem;
        line-height: 1.45;
    }
    .mamamia-help-section h6 {
        margin: 0 0 .3rem;
        color: var(--mm-muted);
        font-size: .66rem;
        font-weight: 700;
        letter-spacing: .07em;
        text-transform: uppercase;
    }
    .mamamia-help-list {
        display: grid;
        grid-template-columns: 10rem minmax(0, 1fr);
        align-items: start;
        gap: .28rem .9rem;
        margin: 0;
        text-align: left;
    }
    .mamamia-help-list dt {
        color: var(--mm-ink);
        font-size: .78rem;
        font-weight: 650;
    }
    .mamamia-help-list dd {
        margin: 0;
        color: var(--mm-muted);
        font-size: .78rem;
        line-height: 1.4;
    }
    .mamamia-help-list dd strong {
        color: var(--mm-ink);
        font-weight: 650;
    }
    .mamamia-help-facts {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: .5rem;
        margin: 0 0 1rem;
    }
    .mamamia-help-facts > div {
        padding: .45rem .55rem;
        border: 1px solid var(--mm-line);
        border-radius: .45rem;
        background: var(--mm-surface);
    }
    .mamamia-help-facts span {
        display: block;
        color: var(--mm-muted);
        font-size: .58rem;
        font-weight: 650;
        letter-spacing: .06em;
        text-transform: uppercase;
    }
    .mamamia-help-facts strong {
        display: block;
        margin-top: .08rem;
        color: var(--mm-ink);
        font-size: .78rem;
        font-weight: 700;
    }
    .mamamia-help-steps {
        margin: 0 0 .65rem;
        padding-left: 1.15rem;
        color: var(--mm-muted);
        font-size: .78rem;
        line-height: 1.5;
    }
    .mamamia-help-steps li + li {
        margin-top: .3rem;
    }
    .mamamia-help-steps strong {
        color: var(--mm-ink);
        font-weight: 650;
    }
    .mamamia-help-command {
        display: block;
        margin: .45rem 0 .65rem;
        padding: .55rem .65rem;
        border: 1px solid rgba(var(--mm-accent-rgb), .35);
        border-left: 3px solid var(--mm-accent);
        border-radius: .45rem;
        background: var(--mm-ink);
        color: var(--mm-paper);
        font-family: ui-monospace, 'SFMono-Regular', 'Cascadia Code',
            'Liberation Mono', Menlo, Consolas, monospace;
        font-size: .74rem;
        line-height: 1.45;
        letter-spacing: .01em;
        overflow-x: auto;
        text-align: left;
        white-space: pre;
        tab-size: 2;
    }
    .mamamia-help-command:focus-visible {
        outline: 2px solid var(--mm-accent);
        outline-offset: 2px;
    }
    .mamamia-tooltip-modal .modal-footer {
        display: flex;
        align-items: center;
        justify-content: flex-end;
        gap: .45rem;
    }
    .mamamia-tooltip-modal .modal-footer .btn {
        min-height: 2rem;
        margin: 0;
        padding: .35rem .75rem;
        font-size: .78rem;
        line-height: 1.25;
    }
    .mamamia-sidebar .form-check {
        margin-bottom: .25rem;
        min-height: 1.25rem;
    }
    .mamamia-sidebar .radio-inline + .radio-inline {
        margin-left: .6rem;
    }
    .mamamia-resolution-slider {
        position: relative;
    }
    .mamamia-resolution-ticks {
        position: relative;
        height: .9rem;
        margin: -.15rem .7rem .15rem;
        color: var(--mm-muted);
        font-size: .62rem;
        line-height: 1;
    }
    .mamamia-resolution-ticks span {
        position: absolute;
        white-space: nowrap;
    }
    .mamamia-resolution-ticks span.last {
        transform: translateX(-100%);
    }
    .mamamia-resolution-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: .45rem;
        margin-bottom: .1rem;
    }
    .mamamia-resolution-head > .shiny-html-output {
        flex: 0 0 auto;
    }
    .mamamia-resolution-head .control-label {
        display: inline-block;
        margin: 0;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 700;
        letter-spacing: .1em;
        text-transform: uppercase;
    }
    .mamamia-resolution-head .form-group {
        flex: 0 0 6rem;
        width: 6rem !important;
        margin-bottom: 0;
    }
    .mamamia-resolution-head .form-control {
        min-height: 1.7rem;
        padding: .15rem .45rem;
        border-color: var(--mm-line);
        border-radius: .45rem;
        background-color: var(--mm-surface);
        font-size: .72rem;
    }
    .mamamia-resolution-slider .shiny-label-null {
        display: none;
    }
    .mamamia-resolution-locked {
        position: relative;
    }
    .mamamia-resolution-locked .shiny-input-container {
        margin-bottom: 0;
        opacity: .55;
        filter: grayscale(1);
    }
    .mamamia-resolution-locked .irs {
        pointer-events: none;
    }
    .mamamia-resolution-lock {
        position: absolute;
        top: 0;
        right: 0;
        bottom: 0;
        left: 0;
        width: 100%;
        height: 100%;
        padding: 0;
        border: 0;
        border-radius: .4rem;
        background: transparent !important;
        color: transparent !important;
        cursor: not-allowed !important;
        transition: none !important;
        z-index: 5;
    }
    .mamamia-resolution-lock:hover,
    .mamamia-resolution-lock:focus,
    .mamamia-resolution-lock:active,
    .mamamia-resolution-lock.active {
        background: transparent !important;
        border-color: transparent !important;
        box-shadow: none !important;
        color: transparent !important;
        cursor: not-allowed !important;
    }
    .mamamia-resolution-lock:focus-visible {
        outline: 2px solid rgba(var(--mm-accent-rgb), .45);
        outline-offset: 2px;
    }
    .mamamia-resolution-note {
        margin: .1rem 0 0;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.35;
    }
    .mamamia-fixed-control {
        display: flex;
        align-items: center;
        min-height: 2rem;
        padding: .28rem .5rem;
        border: 1px solid var(--mm-line);
        border-radius: .42rem;
        background: var(--mm-paper);
        color: var(--mm-muted);
        font-size: .79rem;
        line-height: 1.3;
        cursor: not-allowed;
    }
    .mamamia-fixed-note {
        margin: .15rem 0 .55rem;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.35;
    }
    .mamamia-sidebar .table {
        margin-bottom: .25rem;
        font-size: .68rem;
    }

    /* Call table */
    .mamamia-main .dataTables_wrapper {
        font-size: .75rem;
        overflow-x: hidden;
        overflow-x: clip;
    }
    .mamamia-main .dataTables_wrapper .row:first-child {
        align-items: center;
        margin-bottom: .45rem;
    }
    .mamamia-main .dataTables_length label,
    .mamamia-main .dataTables_filter label {
        color: var(--mm-muted);
        font-size: .68rem;
        font-weight: 550;
    }
    .mamamia-main .dataTables_length .form-select {
        height: 1.85rem;
        padding: .1rem 1.6rem .1rem .45rem;
        border-color: var(--mm-line);
        border-radius: .38rem;
        font-size: .72rem;
    }
    .mamamia-main .dataTables_wrapper .dataTables_filter input {
        height: 1.85rem;
        padding: .1rem .6rem;
        border: 1px solid var(--mm-line);
        border-radius: 999px;
        background: var(--mm-surface);
        font-size: .72rem;
        box-shadow: none;
    }
    .mamamia-main .dataTables_wrapper .dataTables_filter input:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-main table.dataTable {
        border-collapse: separate;
        border-spacing: 0;
    }
    .mamamia-main table.dataTable thead th {
        padding: .45rem .55rem;
        border-bottom: 1px solid var(--mm-line-strong);
        background: var(--mm-paper);
        color: var(--mm-muted);
        font-size: .6rem;
        font-weight: 750;
        letter-spacing: .06em;
        line-height: 1.2;
        text-transform: uppercase;
        vertical-align: bottom;
        white-space: nowrap;
    }
    .mamamia-main table.dataTable thead > tr > th[class*='sorting'] {
        padding-right: 1.5rem;
    }
    @media (max-width: 1320px) {
        .mamamia-main table.dataTable thead > tr > th[class*='sorting'] {
            padding-right: 14px;
        }
        .mamamia-main table.dataTable thead > tr > th[class*='sorting']::before,
        .mamamia-main table.dataTable thead > tr > th[class*='sorting']::after {
            right: 2px;
        }
    }
    .mamamia-main table.dataTable thead th[title] {
        cursor: help;
    }
    .mamamia-main table.dataTable thead th .mamamia-call-label {
        border-bottom: 1px dotted var(--mm-line-strong);
    }
    .mamamia-main table.dataTable tbody td {
        padding: .4rem .55rem;
        border-bottom: 1px solid var(--mm-line);
        font-size: .72rem;
        vertical-align: middle;
        white-space: nowrap;
    }
    .mamamia-main table.dataTable.table-striped > tbody > tr:nth-of-type(odd) > * {
        --bs-table-bg-type: var(--mm-paper);
    }
    .mamamia-main table.dataTable tbody tr:hover > * {
        --bs-table-bg-state: var(--mm-accent-soft);
    }
    .mamamia-main table.dataTable tbody tr.selected > * {
        --bs-table-bg-state: var(--mm-accent-soft);
        background-color: rgba(var(--mm-accent-rgb), .12) !important;
        box-shadow: none !important;
        color: var(--mm-ink) !important;
    }
    .mamamia-main table.dataTable tbody tr.selected > td:first-child {
        box-shadow: inset 3px 0 0 var(--mm-accent) !important;
    }
    .mamamia-main table.dataTable thead tr.filters th {
        padding: .2rem .3rem .35rem;
        background: var(--mm-paper);
    }
    .mamamia-main table.dataTable thead tr.filters input {
        height: 1.65rem;
        padding: .1rem .35rem;
        border: 1px solid var(--mm-line);
        border-radius: .3rem;
        font-size: .66rem;
    }
    .mamamia-main .dataTables_info,
    .mamamia-main .dataTables_paginate {
        margin-top: .55rem;
        color: var(--mm-muted);
        font-size: .68rem;
    }
    .mamamia-main .page-link {
        padding: .2rem .5rem;
        border-color: var(--mm-line);
        color: var(--mm-muted);
        font-size: .7rem;
    }
    .mamamia-main .active > .page-link,
    .mamamia-main .page-link.active {
        border-color: var(--mm-accent);
        background: var(--mm-accent);
    }
    .mamamia-main .dataTables_empty {
        padding: 2rem 1rem !important;
        color: var(--mm-muted);
        font-size: .8rem;
    }

    .mamamia-main .shiny-output-error-validation {
        display: flex;
        align-items: center;
        justify-content: center;
        height: 100%;
        min-height: 12rem;
        padding: 2rem;
        color: var(--mm-muted);
        font-size: .85rem;
        text-align: center;
    }

    /* Analysis progress */
    .shiny-progress-notification {
        padding: .6rem .8rem;
        border: 1px solid var(--mm-line);
        border-radius: .6rem;
        background: var(--mm-surface);
        box-shadow: 0 .6rem 1.6rem rgba(27, 30, 26, .16);
    }
    .shiny-progress-notification .progress {
        height: 4px;
        margin: .3rem 0 .1rem;
        border-radius: 999px;
        background: var(--mm-line);
    }
    .shiny-progress-notification .progress-bar {
        background-color: var(--mm-accent);
    }
    .shiny-progress-notification .progress-text {
        color: var(--mm-ink);
        font-size: .78rem;
    }

    /* Toasts. bslib doubles the class and id selectors, so match that
       specificity to keep the compact spacing. */
    #shiny-notification-panel#shiny-notification-panel {
        right: 1rem;
        bottom: 1rem;
        width: min(23rem, calc(100vw - 2rem));
        padding: 0;
        z-index: 99999;
    }
    .shiny-notification.shiny-notification {
        padding: .5rem 2rem .5rem .75rem;
        margin: 0 0 .4rem;
        border: 1px solid var(--mm-line);
        border-left: 3px solid var(--mm-accent);
        border-radius: .6rem;
        background: var(--mm-surface);
        box-shadow: 0 .6rem 1.6rem rgba(27, 30, 26, .16);
        color: var(--mm-ink);
        font-size: .78rem;
        line-height: 1.4;
    }
    .shiny-notification.shiny-notification:last-child {
        margin-bottom: 0;
    }
    .shiny-notification.shiny-notification .shiny-progress-notification {
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
    }
    .shiny-notification.shiny-notification .shiny-notification-close {
        top: .25rem;
        right: .45rem;
        width: 1.5rem;
        height: 1.5rem;
        padding: 0;
        color: var(--mm-muted);
        font-size: 1rem;
        font-weight: 400;
    }
    .shiny-notification.shiny-notification .shiny-notification-close:hover {
        color: var(--mm-ink);
        font-weight: 400;
    }

    ::selection {
        background: rgba(var(--mm-accent-rgb), .18);
    }
    .btn:focus-visible,
    .form-control:focus-visible,
    .form-select:focus-visible {
        outline: 2px solid var(--mm-accent);
        outline-offset: 2px;
    }
    "
}
