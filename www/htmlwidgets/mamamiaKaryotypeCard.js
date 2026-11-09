HTMLWidgets.widget({
  name: "mamamiaKaryotypeCard",
  type: "output",

  factory: function (el, width, height) {
    var payload = null;
    var plotFont = "'Inter', system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif";
    var legendFontSize = 10;
    var legendRowHeight = 19;
    var legendSwatch = 11;
    var legendGap = 6;
    var legendItemGap = 16;

    function svgEl(name, attrs) {
      var node = document.createElementNS("http://www.w3.org/2000/svg", name);
      if (attrs) {
        Object.keys(attrs).forEach(function (key) {
          if (attrs[key] !== null && attrs[key] !== undefined) {
            node.setAttribute(key, attrs[key]);
          }
        });
      }
      return node;
    }

    function textEl(text, attrs) {
      var node = svgEl("text", attrs);
      node.textContent = text;
      if (!node.getAttribute("font-family")) {
        node.setAttribute("font-family", plotFont);
      }
      return node;
    }

    function sceneSpec() {
      return (payload && payload.scene) || {};
    }

    function sceneLegend() {
      return sceneSpec().legend || {};
    }

    function callLegendEntries() {
      // Legend policy is resolved on the server; this is an adapter accessor.
      return sceneLegend().entries || [];
    }

    function legendItemWidth(entry) {
      return legendSwatch + legendGap + entry.label.length * 5.6;
    }

    function layoutLegend(entries, maxWidth) {
      var rows = [];
      var row = [];
      var width = 0;
      entries.forEach(function (entry) {
        var itemWidth = legendItemWidth(entry);
        var gap = row.length ? legendItemGap : 0;
        if (row.length && width + gap + itemWidth > maxWidth) {
          rows.push({ items: row, width: width });
          row = [];
          width = 0;
          gap = 0;
        }
        row.push({ entry: entry, x: width + gap, width: itemWidth });
        width += gap + itemWidth;
      });
      if (row.length) {
        rows.push({ items: row, width: width });
      }
      return rows;
    }

    function escapeHtml(value) {
      return String(value == null ? "" : value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;");
    }

    function isNumber(value) {
      return typeof value === "number" && Number.isFinite(value);
    }

    function addLegend(svg, x, y, maxWidth) {
      var palette = sceneSpec().palette || payload.palette;
      var spec = sceneLegend();
      var legend = svgEl("g", { "class": "karyotype-legend" });
      svg.appendChild(legend);
      if (spec.kind === "heatmap") {
        var heatmap = spec.heatmap || {};
        var gradientId = "mamamia-karyotype-gradient-" +
          Math.random().toString(36).slice(2);
        var defs = svgEl("defs", {});
        var gradient = svgEl("linearGradient", {
          id: gradientId, x1: "0%", x2: "100%", y1: "0%", y2: "0%"
        });
        [
          ["0%", heatmap.negative],
          ["50%", heatmap.zero],
          ["100%", heatmap.positive]
        ].forEach(function (stop) {
          gradient.appendChild(svgEl("stop", {
            offset: stop[0], "stop-color": stop[1]
          }));
        });
        defs.appendChild(gradient);
        legend.appendChild(defs);
        legend.appendChild(textEl("Score", {
          x: x, y: y, "font-size": 10, "font-weight": 600,
          fill: palette.text
        }));
        legend.appendChild(svgEl("rect", {
          x: x + 36, y: y - 9, width: 126, height: 10,
          fill: "url(#" + gradientId + ")", rx: 1.5
        }));
        var labels = heatmap.labels || {};
        legend.appendChild(textEl(labels.min || "", {
          x: x + 36, y: y + 13, "font-size": 9, "font-weight": 500,
          fill: palette.text
        }));
        legend.appendChild(textEl(labels.zero || "0", {
          x: x + 99, y: y + 13, "font-size": 9, "font-weight": 500,
          fill: palette.text, "text-anchor": "middle"
        }));
        legend.appendChild(textEl(labels.max || "", {
          x: x + 162, y: y + 13, "font-size": 9, "font-weight": 500,
          fill: palette.text, "text-anchor": "end"
        }));
      } else {
        var rows = layoutLegend(callLegendEntries(), maxWidth);
        rows.forEach(function (row, rowIndex) {
          var rowY = y + rowIndex * legendRowHeight;
          row.items.forEach(function (item) {
            var entry = item.entry;
            var itemX = x + item.x;
            legend.appendChild(svgEl("rect", {
              x: itemX, y: rowY - 9, width: legendSwatch,
              height: legendSwatch, fill: entry.fill,
              stroke: entry.stroke,
              "stroke-width": entry.stroke_width,
              "stroke-dasharray": entry.dash,
              rx: 1.5
            }));
            legend.appendChild(textEl(entry.label, {
              x: itemX + legendSwatch + legendGap,
              y: rowY,
              "font-size": legendFontSize,
              "font-weight": 500,
              fill: palette.text
            }));
          });
        });
      }
      var note = sceneSpec().omitted_note;
      if (note && note.text) {
        legend.appendChild(textEl(note.text, {
          x: x + maxWidth, y: y, "text-anchor": "end",
          "font-size": 9, "font-weight": 600, fill: palette.text
        }));
      }
      return legend;
    }

    function hasValue(value) {
      return value !== null && value !== undefined &&
        String(value) !== "" && String(value) !== "NA";
    }

    function drawHeatmapBins(bar, bins) {
      // Bins, including their colour, are resolved in the payload scene by
      // the shared R bin policy. This is only a geometry adapter.
      (bins || []).forEach(function (bin) {
        var start = Number(bin.start);
        var end = Number(bin.end);
        if (!isFinite(start) || !isFinite(end) || end <= start) {
          return;
        }
        bar.group.appendChild(svgEl("rect", {
          x: bar.x,
          y: bar.y + start * bar.height,
          width: bar.width,
          height: Math.max(0.5, (end - start) * bar.height + 0.6),
          fill: bin.fill,
          "shape-rendering": "crispEdges",
          "pointer-events": "none"
        }));
      });
    }

    function drawCallsForBar(bar, calls) {
      // Calls are scene primitives with colour, opacity, stroke and dash
      // already resolved.
      calls.forEach(function (call) {
        var y0 = bar.y + Math.max(0, Math.min(1, call.start)) * bar.height;
        var y1 = bar.y + Math.max(0, Math.min(1, call.end)) * bar.height;
        bar.group.appendChild(svgEl("rect", {
          x: bar.x,
          y: y0,
          width: bar.width,
          height: Math.max(1, y1 - y0),
          fill: call.fill,
          opacity: call.opacity,
          stroke: call.stroke,
          "stroke-width": call.stroke_width,
          "stroke-dasharray": call.dash,
          "pointer-events": "none"
        }));
      });
    }

    function safeColour(value, fallback) {
      return /^#[0-9a-fA-F]{6}$/.test(String(value || ""))
        ? String(value)
        : fallback;
    }

    function fmtSigned(value, digits) {
      if (!isNumber(value)) {
        return "NA";
      }
      var number = Number(value);
      return (number > 0 ? "+" : "") + number.toFixed(digits || 2);
    }

    function fmtBp(value) {
      if (!isNumber(value) || Number(value) < 0) {
        return "NA";
      }
      var number = Number(value);
      if (number >= 1e6) {
        return (number / 1e6).toFixed(1) + " Mb";
      }
      if (number >= 1e3) {
        return (number / 1e3).toFixed(1) + " kb";
      }
      return Math.round(number) + " bp";
    }

    function tipRow(label, value) {
      if (!hasValue(value)) {
        return "";
      }
      return '<div class="mamamia-tip-key">' + escapeHtml(label) +
        '</div><div class="mamamia-tip-value">' +
        escapeHtml(value) + "</div>";
    }

    function tipGrid(rows) {
      var body = (rows || []).filter(function (row) { return row; }).join("");
      return body ? '<div class="mamamia-tip-grid">' + body + "</div>" : "";
    }

    function tipSection(label, body) {
      if (!body) {
        return "";
      }
      return '<div class="mamamia-tip-section">' +
        '<div class="mamamia-tip-label">' + escapeHtml(label) + "</div>" +
        body + "</div>";
    }

    function tipChip(chip) {
      if (!chip || !hasValue(chip.text)) {
        return "";
      }
      return '<span class="mamamia-tip-chip ' +
        escapeHtml(chip.kind || "neutral") + '">' +
        escapeHtml(chip.text) + "</span>";
    }

    function tipMetric(metric) {
      return '<div class="mamamia-tip-metric"><span>' +
        escapeHtml(metric.label) + '</span><strong class="' +
        escapeHtml(metric.class || "muted") + '">' +
        escapeHtml(metric.value) + "</strong></div>";
    }

    function tooltipHeader(tooltip) {
      var accent = safeColour(tooltip.accent, "#6c757d");
      var status = tooltip.status || {};
      return '<div class="mamamia-tip" style="border-left-color:' +
        accent + ";--tip-accent:" + accent + '">' +
        '<div class="mamamia-tip-head">' +
        '<div class="mamamia-tip-title"><span class="mamamia-tip-class">' +
        escapeHtml(tooltip.title || "call") +
        '</span><span class="mamamia-tip-status ' +
        escapeHtml(status.class || "mamamia-tip-status-neutral") + '">' +
        escapeHtml(status.text || "call") + "</span></div>" +
        '<div class="mamamia-tip-sub">' +
        escapeHtml(tooltip.subtitle || "") + "</div></div>";
    }

    function tooltipSections(sections) {
      var html = "";
      (sections || []).forEach(function (section) {
        var body = "";
        if (section.rows) {
          body = tipGrid(section.rows.map(function (row) {
            return tipRow(row.label, row.value);
          }));
        } else if (section.metrics) {
          body = '<div class="mamamia-tip-metrics">' +
            section.metrics.map(tipMetric).join("") + "</div>";
          var chips = (section.chips || []).map(tipChip).join("");
          if (chips) {
            body += '<div class="mamamia-tip-chips">' + chips + "</div>";
          }
        }
        if (!body && section.empty) {
          body = '<div class="mamamia-tip-empty">' +
            escapeHtml(section.empty) + "</div>";
        }
        html += tipSection(section.label, body);
      });
      return html;
    }

    function sourceFraction(fraction, chromosome) {
      var removed = chromosome.deleted_intervals || [];
      if (!removed.length) {
        return fraction;
      }
      var surviving = 1;
      removed.forEach(function (interval) {
        surviving -= Math.max(0, Number(interval.end) - Number(interval.start));
      });
      if (!(surviving > 0)) {
        return fraction;
      }
      var target = fraction * surviving;
      var cursor = 0;
      var removedTotal = 0;
      for (var i = 0; i < removed.length; i += 1) {
        var from = Number(removed[i].start);
        var to = Number(removed[i].end);
        if (!isNumber(from) || !isNumber(to)) {
          continue;
        }
        if (from > cursor) {
          var displayEnd = (from - removedTotal) / surviving;
          if (fraction <= displayEnd) {
            return Math.max(cursor, Math.min(from, target + removedTotal));
          }
        }
        removedTotal += Math.max(0, to - from);
        cursor = Math.max(cursor, to);
      }
      return Math.max(cursor, Math.min(1, target + removedTotal));
    }

    function describe(chromosome, fraction, calls) {
      var cache = calls && calls.__tipCache;
      if (cache && cache.call &&
          fraction >= cache.call.start && fraction <= cache.call.end) {
        return cache.html;
      }
      var matched = null;
      (calls || []).forEach(function (call) {
        if (fraction >= call.start && fraction <= call.end) {
          matched = call;
        }
      });
      var html;
      if (matched && matched.tooltip) {
        // The server resolves every display string for the matched-call card.
        // This adapter only renders the tooltip model.
        html = tooltipHeader(matched.tooltip) +
          tooltipSections(matched.tooltip.sections);
      } else if (matched) {
        // Deliberately minimal fallback for a legacy payload without a
        // tooltip. Canonical call fields are not reconstructed here.
        var legacyParts = [];
        if (hasValue(payload.sample)) {
          legacyParts.push(payload.sample);
        }
        if (hasValue(chromosome.group)) {
          legacyParts.push("group " + chromosome.group);
        }
        html = tooltipHeader({
          accent: "#6c757d",
          title: "call",
          status: {
            text: "call",
            class: "mamamia-tip-status-neutral"
          },
          subtitle: legacyParts.join(" · ")
        }) + tipSection(
          "Call",
          '<div class="mamamia-tip-empty">' +
            "Call details unavailable in this legacy payload.</div>"
        );
      } else {
        var subParts = [];
        if (hasValue(payload.sample)) {
          subParts.push(payload.sample);
        }
        if (hasValue(chromosome.group)) {
          subParts.push("group " + chromosome.group);
        }
        if (isNumber(chromosome.length)) {
          subParts.push(fmtBp(sourceFraction(fraction, chromosome) *
            chromosome.length));
        }
        var className = chromosome.unplaced
          ? "unplaced donor"
          : chromosome.deleted
            ? "deleted recipient"
            : chromosome.added
              ? "added donor"
              : "recipient chromosome";
        html = tooltipHeader({
          accent: "#6c757d",
          title: className,
          status: {
            text: "no call",
            class: "mamamia-tip-status-neutral"
          },
          subtitle: subParts.join(" · ")
        });
        var chromosomeRows = [];
        if (hasValue(chromosome.component)) {
          chromosomeRows.push(tipRow("Component", chromosome.component));
        }
        if (chromosome.unplaced === true) {
          chromosomeRows.push(tipRow(
            "State", "unplaced donor chromosome"
          ));
        } else if (chromosome.deleted === true) {
          chromosomeRows.push(tipRow(
            "State", "deleted recipient (no coverage)"
          ));
        } else if (chromosome.added === true) {
          chromosomeRows.push(tipRow("State", "added donor chromosome"));
        } else {
          chromosomeRows.push(tipRow("State", "recipient chromosome"));
        }
        if (isNumber(chromosome.length)) {
          chromosomeRows.push(tipRow("Length", fmtBp(chromosome.length)));
        }
        html += tipSection("Chromosome", tipGrid(chromosomeRows));
        if (payload.mode === "calls") {
          html += tipSection(
            "Call",
            '<div class="mamamia-tip-empty">' +
              "No assigned donor segment at this position.</div>"
          );
        }
      }

      if (payload.mode === "heatmap") {
        var scores = chromosome.scores || [];
        var index = scores.length
          ? Math.max(0, Math.min(
            scores.length - 1,
            Math.floor(fraction * scores.length)
          ))
          : -1;
        var position = [
          tipRow("Position", fmtBp(
            sourceFraction(fraction, chromosome) * (chromosome.length || 0)
          )),
          tipRow(
            "Score",
            index >= 0 ? fmtSigned(scores[index]) : "NA"
          )
        ];
        html += tipSection("Position", tipGrid(position));
      }

      html += '<div class="mamamia-tip-footer">' +
        "Click to filter this group</div></div>";
      if (matched && calls) {
        calls.__tipCache = {call: matched, html: html};
      }
      return html;
    }

    function sendClick(group) {
      if (window.Shiny && el.id) {
        window.Shiny.setInputValue(
          el.id + "_click",
          { group: group },
          { priority: "event" }
        );
      }
    }

    function exportFileName(kind) {
      var sample = payload && payload.sample
        ? String(payload.sample)
        : "sample";
      sample = sample.replace(/[^A-Za-z0-9._-]+/g, "-");
      var mode = payload && payload.mode ? String(payload.mode) : "karyotype";
      return "mamamia-karyotype-" + sample + "-" + mode + "-" + kind + ".svg";
    }

    function svgWithBackground(svg) {
      var clone = svg.cloneNode(true);
      clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      var svgWidth = Number(clone.getAttribute("width"));
      var svgHeight = Number(clone.getAttribute("height"));
      if (!(svgWidth > 0)) {
        svgWidth = 800;
      }
      if (!(svgHeight > 0)) {
        svgHeight = 520;
      }
      clone.insertBefore(svgEl("rect", {
        x: 0,
        y: 0,
        width: svgWidth,
        height: svgHeight,
        fill: "#ffffff"
      }), clone.firstChild);
      return clone;
    }

    function downloadSvg(svg, filename) {
      var source = '<?xml version="1.0" encoding="UTF-8"?>\n' +
        new XMLSerializer().serializeToString(svgWithBackground(svg));
      var blob = new Blob([source], { type: "image/svg+xml;charset=utf-8" });
      var url = URL.createObjectURL(blob);
      var link = document.createElement("a");
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      link.remove();
      window.setTimeout(function () {
        URL.revokeObjectURL(url);
      }, 1000);
    }

    function fullExportPayload() {
      // The current contract requires a self-contained full view. A payload
      // without one has no full export; positional side channels are not
      // reconstructed.
      var view = payload.views && payload.views.full;
      if (!view || !(view.chromosomes || []).length) {
        return null;
      }
      var full = {};
      Object.keys(view).forEach(function (key) {
        full[key] = view[key];
      });
      full.chromosomes = view.chromosomes.map(function (chromosome) {
        var copy = {};
        Object.keys(chromosome).forEach(function (key) {
          copy[key] = chromosome[key];
        });
        copy.selected = true;
        return copy;
      });
      full.hide_deletions = false;
      full.hidden_deletions = { segments: 0, bases: 0 };
      full.views = null;
      return full;
    }

    function renderFullSvg() {
      var fullPayload = fullExportPayload();
      if (!fullPayload) {
        return null;
      }
      var fullWidth = Math.max(1600, Number(width) || 0);
      var fullHeight = Math.max(900, Number(height) || 0);
      var holder = document.createElement("div");
      holder.style.position = "fixed";
      holder.style.left = "-100000px";
      holder.style.top = "0";
      holder.style.width = fullWidth + "px";
      holder.style.height = fullHeight + "px";
      document.body.appendChild(holder);
      var savedEl = el;
      var savedPayload = payload;
      var savedWidth = width;
      var savedHeight = height;
      el = holder;
      payload = fullPayload;
      width = fullWidth;
      height = fullHeight;
      var svg = null;
      try {
        render();
        svg = holder.querySelector("svg");
      } finally {
        el = savedEl;
        payload = savedPayload;
        width = savedWidth;
        height = savedHeight;
        holder.remove();
      }
      return svg ? svg.cloneNode(true) : null;
    }

    function downloadKaryotype(kind) {
      if (!payload || !(payload.chromosomes || []).length) {
        return false;
      }
      var svg = kind === "full"
        ? renderFullSvg()
        : el.querySelector("svg");
      if (!svg) {
        return false;
      }
      downloadSvg(svg, exportFileName(kind === "full" ? "full" : "current"));
      return true;
    }

    function render() {
      if (!payload) {
        return;
      }
      var hasChromosomes = (payload.chromosomes || []).length > 0;
      var buttonPrefix = el.id ? el.id.replace(/-karyotype$/, "") : "";
      ["current", "full"].forEach(function (kind) {
        var button = document.getElementById(
          buttonPrefix + "-download_karyotype_" + kind
        );
        if (button) {
          button.disabled = !hasChromosomes;
        }
      });
      el.innerHTML = "";
      var root = document.createElement("div");
      root.className = "mamamia-karyotype-root";
      root.style.width = "100%";
      root.style.height = "100%";
      root.style.overflow = "auto";
      el.appendChild(root);

      var tooltip = document.createElement("div");
      tooltip.className = "mamamia-karyotype-tooltip";
      tooltip.style.display = "none";
      root.appendChild(tooltip);

      var chromosomes = payload.chromosomes || [];
      if (!chromosomes.length) {
        var empty = document.createElement("div");
        empty.className = "mamamia-karyotype-empty";

        var mark = document.createElement("div");
        mark.className = "mamamia-empty-mark";
        mark.setAttribute("aria-hidden", "true");
        mark.innerHTML =
          '<svg viewBox="0 0 32 32" fill="none" ' +
          'xmlns="http://www.w3.org/2000/svg">' +
          '<path d="M11 6 C 11 9, 11 12, 16 16 C 21 20, 21 23, 21 26" ' +
          'stroke="currentColor" stroke-width="2.2" ' +
          'stroke-linecap="round"/>' +
          '<path d="M21 6 C 21 9, 21 12, 16 16 C 11 20, 11 23, 11 26" ' +
          'stroke="#b7791f" stroke-width="2.2" ' +
          'stroke-linecap="round"/></svg>';
        empty.appendChild(mark);

        var title = document.createElement("p");
        title.className = "mamamia-empty-title";
        title.textContent = payload.message ||
          "No chromosomes in this result.";
        empty.appendChild(title);

        if (payload.hint) {
          var hint = document.createElement("p");
          hint.className = "mamamia-empty-hint";
          hint.textContent = payload.hint;
          empty.appendChild(hint);
        }

        root.appendChild(empty);
        return;
      }

      var box = el.getBoundingClientRect();
      var availableWidth = Math.max(width || box.width || 800, 480);
      var availableHeight = Math.max(height || box.height || 520, 420);
      var sideMargin = 22;
      var legendWidth = Math.max(0, availableWidth - sideMargin * 2);
      var legendRows = layoutLegend(
        callLegendEntries(), legendWidth
      ).length;
      var margin = {
        top: 58 + Math.max(0, legendRows - 1) * legendRowHeight,
        right: sideMargin,
        bottom: 36,
        left: sideMargin
      };
      var maxBarHeight = Math.max(
        220,
        availableHeight - margin.top - margin.bottom
      );
      var totalHeight = margin.top + maxBarHeight + margin.bottom;
      var sceneChromosomes = sceneSpec().chromosomes || [];
      var maxLength = 0;
      chromosomes.forEach(function (chromosome, index) {
        var sceneChromosome = sceneChromosomes[index] || {};
        var length = Number(sceneChromosome.display_length);
        if (isFinite(length) && length > maxLength) {
          maxLength = length;
        }
      });
      if (!(maxLength > 0)) {
        maxLength = 1;
      }
      var pairGap = 3;
      var chromosomeGap = 8;
      var maxCopies = 1;
      var nPairs = 0;
      chromosomes.forEach(function (chromosome, index) {
        var sceneChromosome = sceneChromosomes[index] || {};
        var copies = Math.max(1, Number(sceneChromosome.copies) || 1);
        if (copies > maxCopies) {
          maxCopies = copies;
        }
        if (copies === 1 || Number(sceneChromosome.copy) === 1) {
          nPairs += 1;
        }
      });
      var slotWidth = (availableWidth - margin.left - margin.right) /
        Math.max(1, nPairs);
      var barWidth = Math.max(
        7,
        Math.min(
          34,
          (slotWidth - Math.max(0, maxCopies - 1) * pairGap -
            chromosomeGap) / maxCopies
        )
      );
      var barsWidth = 0;
      chromosomes.forEach(function (chromosome, index) {
        barsWidth += barWidth;
        if (index < chromosomes.length - 1) {
          var sceneChromosome = sceneChromosomes[index] || {};
          var copies = Math.max(1, Number(sceneChromosome.copies) || 1);
          var copy = Math.max(1, Number(sceneChromosome.copy) || 1);
          barsWidth += (copies > 1 && copy < copies)
            ? pairGap
            : chromosomeGap;
        }
      });
      var plotWidth = Math.max(
        availableWidth - margin.left - margin.right, barsWidth
      );
      var totalWidth = margin.left + plotWidth + margin.right;
      var svg = svgEl("svg", {
        width: totalWidth,
        height: totalHeight,
        viewBox: "0 0 " + totalWidth + " " + totalHeight,
        role: "img",
        "aria-label": "Karyotype in assembly order"
      });
      svg.style.display = "block";
      root.appendChild(svg);

      addLegend(svg, margin.left, 16, legendWidth);

      var cursorX = margin.left + (plotWidth - barsWidth) / 2;
      var y = margin.top;
      chromosomes.forEach(function (chromosome, index) {
        var sceneChromosome = sceneChromosomes[index] || {};
        var added = sceneChromosome.added === true;
        var deleted = sceneChromosome.deleted === true;
        var barCalls = chromosome.calls || [];
        var styledCalls = sceneChromosome.calls || [];
        var copies = Math.max(1, Number(sceneChromosome.copies) || 1);
        var copy = Math.max(1, Number(sceneChromosome.copy) || 1);
        var height = Math.max(
          3, maxBarHeight * (Number(
            sceneChromosome.display_length
          ) || 0) / maxLength
        );
        var x = cursorX;
        var pairCenter = x + (copies > 1
          ? barWidth + pairGap / 2
          : barWidth / 2);

        var chromGroup = svgEl("g", {
          "class": "karyotype-chromosome",
          "data-chromosome": chromosome.id,
          "data-group": chromosome.group,
          "data-copy": copy
        });
        if (sceneChromosome.selected === false) {
          chromGroup.setAttribute("opacity", "0.32");
        }
        svg.appendChild(chromGroup);

        var baseSpec = sceneChromosome.base || {};
        var base = svgEl("rect", {
          x: x,
          y: y,
          width: barWidth,
          height: height,
          fill: baseSpec.fill,
          stroke: baseSpec.stroke,
          "stroke-width": baseSpec.stroke_width,
          "stroke-dasharray": baseSpec.dash,
          rx: 2
        });
        chromGroup.appendChild(base);

        var barBox = {
          group: chromGroup,
          x: x,
          y: y,
          width: barWidth,
          height: height
        };
        if (barCalls.length) {
          if (sceneSpec().mode === "heatmap" && !deleted) {
            drawHeatmapBins(barBox, sceneChromosome.heatmap_bins);
            var crossCalls = styledCalls.filter(function (call) {
              return call.cross_group === true;
            });
            if (crossCalls.length) {
              drawCallsForBar(barBox, crossCalls);
            }
          } else {
            drawCallsForBar(barBox, styledCalls);
          }
        }

        if (sceneChromosome.caption) {
          chromGroup.appendChild(textEl(sceneChromosome.caption.text, {
            x: pairCenter,
            y: margin.top - 7,
            "text-anchor": "middle",
            "font-size": 9,
            "font-weight": 700,
            "letter-spacing": "0.03em",
            fill: sceneChromosome.caption.fill
          }));
        }

        if (sceneChromosome.show_label === true) {
          chromGroup.appendChild(textEl(
            sceneChromosome.label || "",
            {
              x: pairCenter,
              y: margin.top + maxBarHeight + 17,
              "text-anchor": "middle",
              "font-size": 11,
              "font-weight": 600,
              fill: (sceneSpec().palette || payload.palette).text
            }
          ));
        }

        var hit = svgEl("rect", {
          x: x - 2,
          y: y - 2,
          width: barWidth + 4,
          height: height + 4,
          fill: "transparent",
          "class": "karyotype-hit"
        });
        hit.style.cursor = "pointer";
        chromGroup.appendChild(hit);
        hit.addEventListener("mousemove", function (event) {
          var rect = svg.getBoundingClientRect();
          var pointX = event.clientX - rect.left;
          var pointY = event.clientY - rect.top;
          var fraction = Math.max(
            0, Math.min(1, (pointY - y) / height)
          );
          tooltip.innerHTML = describe(chromosome, fraction, barCalls);
          tooltip.style.display = "block";
          var pad = 10;
          var gap = 14;
          var tipWidth = tooltip.offsetWidth || 260;
          var tipHeight = tooltip.offsetHeight || 120;
          var viewX = event.clientX;
          var viewY = event.clientY;
          var left = viewX + gap;
          var top = viewY + gap;
          if (left + tipWidth > window.innerWidth - pad) {
            left = viewX - tipWidth - gap;
          }
          if (top + tipHeight > window.innerHeight - pad) {
            top = viewY - tipHeight - gap;
          }
          var rootRect = root.getBoundingClientRect();
          var minLeft = rootRect.left + pad;
          var maxLeft = rootRect.right - tipWidth - pad;
          var minTop = rootRect.top + pad;
          var maxTop = rootRect.bottom - tipHeight - pad;
          left = Math.max(minLeft, Math.min(left, Math.max(minLeft, maxLeft)));
          top = Math.max(minTop, Math.min(top, Math.max(minTop, maxTop)));
          tooltip.style.left = left + "px";
          tooltip.style.top = top + "px";
        });
        hit.addEventListener("mouseleave", function () {
          tooltip.style.display = "none";
        });
        hit.addEventListener("click", function (event) {
          event.stopPropagation();
          sendClick(chromosome.group);
        });

        if (index < chromosomes.length - 1) {
          cursorX += barWidth + ((copies > 1 && copy < copies)
            ? pairGap
            : chromosomeGap);
        }
      });

      svg.addEventListener("click", function () {
        sendClick(null);
      });
    }

    el.mamamiaDownloadKaryotype = downloadKaryotype;

    return {
      renderValue: function (x) {
        payload = x;
        render();
      },
      resize: function (newWidth, newHeight) {
        width = newWidth;
        height = newHeight;
        render();
      }
    };
  }
});
