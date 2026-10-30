HTMLWidgets.widget({
  name: "mamamiaKaryotypeCard",
  type: "output",

  factory: function (el, width, height) {
    var payload = null;

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
      return node;
    }

    function escapeHtml(value) {
      return String(value == null ? "" : value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;");
    }

    function hexToRgb(hex) {
      var value = hex.replace("#", "");
      if (value.length === 3) {
        value = value[0] + value[0] + value[1] + value[1] + value[2] + value[2];
      }
      return [
        parseInt(value.slice(0, 2), 16),
        parseInt(value.slice(2, 4), 16),
        parseInt(value.slice(4, 6), 16)
      ];
    }

    function mixColour(from, to, t) {
      var a = hexToRgb(from);
      var b = hexToRgb(to);
      var c = a.map(function (v, i) {
        return Math.round(v + (b[i] - v) * Math.max(0, Math.min(1, t)));
      });
      return "rgb(" + c[0] + "," + c[1] + "," + c[2] + ")";
    }

    function scoreColour(value, limit, palette) {
      if (!isNumber(value)) {
        return "#eeeeee";
      }
      var t = Math.max(-1, Math.min(1, value / limit));
      if (t < 0) {
        return mixColour(palette.zero, palette.negative, -t);
      }
      return mixColour(palette.zero, palette.positive, t);
    }

    function isNumber(value) {
      return typeof value === "number" && Number.isFinite(value);
    }

    function fmt(value) {
      return isNumber(value) ? Number(value).toFixed(2) : "NA";
    }

    function addLegend(svg, x, y, payload) {
      var palette = payload.palette;
      var legend = svgEl("g", { "class": "karyotype-legend" });
      svg.appendChild(legend);
      if (payload.mode === "heatmap") {
        var gradientId = "mamamia-karyotype-gradient-" +
          Math.random().toString(36).slice(2);
        var defs = svgEl("defs", {});
        var gradient = svgEl("linearGradient", {
          id: gradientId, x1: "0%", x2: "100%", y1: "0%", y2: "0%"
        });
        [
          ["0%", palette.negative],
          ["50%", palette.zero],
          ["100%", palette.positive]
        ].forEach(function (stop) {
          gradient.appendChild(svgEl("stop", {
            offset: stop[0], "stop-color": stop[1]
          }));
        });
        defs.appendChild(gradient);
        legend.appendChild(defs);
        legend.appendChild(textEl("Score", {
          x: x, y: y, "font-size": 9, fill: palette.text
        }));
        legend.appendChild(svgEl("rect", {
          x: x + 34, y: y - 8, width: 120, height: 9,
          fill: "url(#" + gradientId + ")", rx: 1
        }));
        legend.appendChild(textEl("-" + Number(payload.score_limit).toFixed(1), {
          x: x + 34, y: y + 12, "font-size": 8, fill: palette.text
        }));
        legend.appendChild(textEl("0", {
          x: x + 94, y: y + 12, "font-size": 8, fill: palette.text,
          "text-anchor": "middle"
        }));
        legend.appendChild(textEl("+" + Number(payload.score_limit).toFixed(1), {
          x: x + 154, y: y + 12, "font-size": 8, fill: palette.text,
          "text-anchor": "end"
        }));
      } else {
        var entries = [
          ["substitution", "substitution"],
          ["addition", "addition"],
          ["deletion", "deletion"],
          ["complex", "complex"],
          ["other", "other"],
          ["deleted", "deleted"],
          ["added", "added donor"]
        ];
        var cursor = x;
        entries.forEach(function (entry) {
          var colour = palette[entry[0]];
          legend.appendChild(svgEl("rect", {
            x: cursor, y: y - 8, width: 10, height: 10,
            fill: colour,
            stroke: entry[0] === "deleted" ? palette.negative : palette.text,
            "stroke-width": entry[0] === "deleted" ? 1.2 : 0.6,
            "stroke-dasharray": entry[0] === "deleted" ? "3 2" : null,
            rx: 1
          }));
          legend.appendChild(textEl(entry[1], {
            x: cursor + 13, y: y, "font-size": 8, fill: palette.text
          }));
          cursor += 20 + entry[1].length * 5;
        });
      }
      return legend;
    }

    function hasValue(value) {
      return value !== null && value !== undefined &&
        String(value) !== "" && String(value) !== "NA";
    }

    function callOpacity(call) {
      var opacity = call.confidence === "high" ? 0.95 :
        (call.confidence === "moderate" ? 0.72 : 0.45);
      if (call.pass === false) {
        opacity *= 0.45;
      }
      if (call.recipient_assignment === "ambiguous") {
        opacity *= 0.55;
      }
      if (call.recipient_assignment === "no_recipient_loss") {
        opacity *= 0.35;
      }
      return opacity;
    }

    function drawHeatmapForBar(chromosome, bar, calls) {
      var scores = chromosome.scores || [];
      calls = calls.filter(function (call) {
        return call.cross_group !== true;
      });
      if (!scores.length || !calls.length) {
        return;
      }
      var bins = Math.min(
        scores.length, Math.max(40, Math.floor(bar.height / 2))
      );
      var binHeight = bar.height / bins;
      var intervals = calls.map(function (call) {
        return {
          start: call.start || 0,
          end: call.end === null || call.end === undefined ? 1 : call.end
        };
      });
      for (var i = 0; i < bins; i += 1) {
        var centre = (i + 0.5) / bins;
        var inside = intervals.some(function (interval) {
          return centre >= interval.start && centre <= interval.end;
        });
        if (!inside) {
          continue;
        }
        var from = Math.floor(i * scores.length / bins);
        var to = Math.max(
          from + 1, Math.floor((i + 1) * scores.length / bins)
        );
        var sum = 0;
        var count = 0;
        for (var j = from; j < to; j += 1) {
          var value = scores[j];
          if (value !== null && value !== undefined && isFinite(value)) {
            sum += value;
            count += 1;
          }
        }
        var score = count ? sum / count : NaN;
        bar.group.appendChild(svgEl("rect", {
          x: bar.x,
          y: bar.y + i * binHeight,
          width: bar.width,
          height: binHeight + 0.6,
          fill: scoreColour(score, payload.score_limit, payload.palette),
          "shape-rendering": "crispEdges",
          "pointer-events": "none"
        }));
      }
    }

    function drawCallsForBar(bar, calls) {
      calls.forEach(function (call) {
        var y0 = bar.y + Math.max(0, Math.min(1, call.start)) * bar.height;
        var y1 = bar.y + Math.max(0, Math.min(1, call.end)) * bar.height;
        var colour = payload.palette[call.class] || payload.palette.other;
        var cross = call.cross_group === true;
        bar.group.appendChild(svgEl("rect", {
          x: bar.x,
          y: y0,
          width: bar.width,
          height: Math.max(1, y1 - y0),
          fill: colour,
          opacity: callOpacity(call),
          stroke: cross ? "#333333" : "#ffffff",
          "stroke-width": cross ? 1.3 : 0.5,
          "stroke-dasharray": cross ? "3 2" : null,
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

    function signClass(value) {
      if (!isNumber(value)) {
        return "muted";
      }
      return Number(value) > 0 ? "pos" : Number(value) < 0 ? "neg" : "muted";
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

    function fmtSpan(start, end) {
      if (isNumber(start) && isNumber(end)) {
        return fmtBp(start) + " – " + fmtBp(end);
      }
      if (isNumber(start)) {
        return "from " + fmtBp(start);
      }
      if (isNumber(end)) {
        return "to " + fmtBp(end);
      }
      return "";
    }

    function fmtPercent(value) {
      return isNumber(value) ? Math.round(100 * Number(value)) + "%" : "NA";
    }

    function statusClass(value) {
      var key = String(value || "").toLowerCase();
      if (key === "confident") {
        return "mamamia-tip-status-confident";
      }
      if (key === "ambiguous") {
        return "mamamia-tip-status-ambiguous";
      }
      if (key === "needs review" || key === "needs_review") {
        return "mamamia-tip-status-review";
      }
      return "mamamia-tip-status-neutral";
    }

    function resolvedStatus(call) {
      if (hasValue(call.class_status)) {
        return call.class_status;
      }
      if (call.class_confidence === "high") {
        return "confident";
      }
      if (call.class_confidence === "moderate") {
        return "ambiguous";
      }
      if (call.class_confidence === "low") {
        return "needs_review";
      }
      return "call";
    }

    function statusLabel(value) {
      return String(value || "call").replace(/_/g, " ");
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

    function tipChip(value, kind) {
      if (!hasValue(value)) {
        return "";
      }
      return '<span class="mamamia-tip-chip ' +
        escapeHtml(kind || "neutral") + '">' +
        escapeHtml(value) + "</span>";
    }

    function tipMetric(label, value) {
      return '<div class="mamamia-tip-metric"><span>' +
        escapeHtml(label) + '</span><strong class="' +
        signClass(value) + '">' + escapeHtml(fmtSigned(value)) +
        "</strong></div>";
    }

    function assignmentLabel(value) {
      var key = String(value || "");
      if (key === "assigned") {
        return "assigned";
      }
      if (key === "cross_group") {
        return "cross-group candidate";
      }
      if (key === "ambiguous") {
        return "ambiguous";
      }
      if (key === "no_recipient_loss") {
        return "no host";
      }
      return key || "group only";
    }

    function describe(chromosome, fraction, calls) {
      var matched = null;
      (calls || []).forEach(function (call) {
        if (fraction >= call.start && fraction <= call.end) {
          matched = call;
        }
      });
      var accent = safeColour(
        matched ? payload.palette[matched.class] : null,
        "#6c757d"
      );
      var className = matched
        ? (matched.class || "call")
        : chromosome.deleted
          ? "deleted recipient"
          : chromosome.added
            ? "added donor"
            : "recipient chromosome";
      var status = matched ? resolvedStatus(matched) : "no call";
      var subParts = [];
      if (hasValue(payload.sample)) {
        subParts.push(payload.sample);
      }
      if (hasValue(chromosome.group)) {
        subParts.push("group " + chromosome.group);
      }
      if (matched) {
        var callSpan = fmtSpan(
          matched.donor_start_bp,
          matched.donor_end_bp
        );
        var spanLabel = "donor ";
        if (!callSpan) {
          callSpan = fmtSpan(matched.start_bp, matched.end_bp);
          spanLabel = "";
        }
        if (callSpan) {
          subParts.push(spanLabel + callSpan);
        }
      } else if (isNumber(chromosome.length)) {
        subParts.push(fmtBp(fraction * chromosome.length));
      }

      var html = '<div class="mamamia-tip" style="border-left-color:' +
        accent + ";--tip-accent:" + accent + '">' +
        '<div class="mamamia-tip-head">' +
        '<div class="mamamia-tip-title"><span class="mamamia-tip-class">' +
        escapeHtml(className) + '</span><span class="mamamia-tip-status ' +
        statusClass(status) + '">' +
        escapeHtml(statusLabel(status)) + "</span></div>" +
        '<div class="mamamia-tip-sub">' +
        escapeHtml(subParts.join(" · ")) + "</div></div>";

      if (matched) {
        var recipient = hasValue(matched.recipient_chr_name)
          ? matched.recipient_chr_name
          : (matched.recipient_assignment === "no_recipient_loss"
            ? "no host"
            : "unassigned");
        if (hasValue(matched.recipient_group)) {
          recipient += " · group " + matched.recipient_group;
        }
        var recipientSpan = fmtSpan(
          matched.recipient_start_bp,
          matched.recipient_end_bp
        );
        if (recipientSpan) {
          recipient += " · " + recipientSpan;
        }
        var donor = hasValue(matched.donor_chr_name)
          ? matched.donor_chr_name
          : (hasValue(matched.donor_group)
            ? "group " + matched.donor_group
            : "unknown");
        if (hasValue(matched.donor_chr_name) &&
            hasValue(matched.donor_group)) {
          donor += " · group " + matched.donor_group;
        }
        var donorSpan = fmtSpan(matched.donor_start_bp, matched.donor_end_bp);
        if (donorSpan) {
          donor += " · " + donorSpan;
        }
        var placement = [
          tipRow("Recipient", recipient),
          tipRow("Donor", donor),
          tipRow(
            "Relation",
            matched.recipient_assignment === "no_recipient_loss"
              ? "not placed"
              : (matched.cross_group ? "cross-group" : "same group")
          ),
          tipRow(
            "Assignment",
            assignmentLabel(matched.recipient_assignment) +
              (isNumber(matched.recipient_margin)
                ? " · margin " + fmtSigned(matched.recipient_margin)
                : "")
          )
        ];
        if (hasValue(matched.class_group) &&
            matched.class_group !== matched.class) {
          placement.push(tipRow(
            "Group class",
            matched.class_group + " → " + matched.class
          ));
        }
        if (matched.recipient_assignment === "ambiguous" &&
            hasValue(matched.recipient_alternative)) {
          placement.push(tipRow(
            "Alternative",
            matched.recipient_alternative
          ));
        }
        html += tipSection("Placement", tipGrid(placement));

        var metrics = '<div class="mamamia-tip-metrics">' +
          tipMetric("Score", matched.score) +
          tipMetric("Donor", matched.donor_mean) +
          tipMetric("Recipient", matched.recipient_mean) +
          tipMetric("Total", matched.total_mean) +
          "</div>";
        var chips = [];
        if (hasValue(matched.class_confidence)) {
          chips.push(tipChip(
            "class " + matched.class_confidence, "confidence"
          ));
        }
        if (hasValue(matched.class_evidence)) {
          chips.push(tipChip(matched.class_evidence, "evidence"));
        }
        if (hasValue(matched.copy_state)) {
          chips.push(tipChip(matched.copy_state, "copy"));
        }
        if (hasValue(matched.total_evidence) &&
            matched.total_evidence !== "neutral") {
          chips.push(tipChip(
            "total " + matched.total_evidence, "copy"
          ));
        }
        if (matched.component_status === "heterogeneous") {
          chips.push(tipChip("components heterogeneous", "warn"));
        }
        if (matched.mappability_status === "low") {
          chips.push(tipChip("low mappability", "warn"));
        }
        if (matched.pass === false) {
          chips.push(tipChip("not passing", "warn"));
        }
        if (chips.length) {
          metrics += '<div class="mamamia-tip-chips">' +
            chips.join("") + "</div>";
        }
        html += tipSection("Evidence", metrics);

        var dosage = [];
        if (isNumber(matched.affected_copies)) {
          var copies = matched.affected_copies +
            (isNumber(matched.recipient_copy_number)
              ? " of " + matched.recipient_copy_number + " copies"
              : " copies");
          dosage.push(tipRow("Copies", copies));
        } else {
          dosage.push(tipRow("Copies", "unavailable"));
        }
        if (isNumber(matched.affected_fraction)) {
          dosage.push(tipRow(
            "Fraction", fmtPercent(matched.affected_fraction)
          ));
        }
        if (isNumber(matched.donor_copies)) {
          dosage.push(tipRow("Donor copies", fmt(matched.donor_copies)));
        }
        if (hasValue(matched.dosage_confidence)) {
          dosage.push(tipRow("Confidence", matched.dosage_confidence));
        }
        html += tipSection("Dosage", tipGrid(dosage));

        var geometry = [];
        if (isNumber(matched.width_bp)) {
          geometry.push(tipRow("Donor width", fmtBp(matched.width_bp)));
        }
        if (isNumber(matched.recipient_deletion_width_bp)) {
          geometry.push(tipRow(
            "Recipient run", fmtBp(matched.recipient_deletion_width_bp)
          ));
        }
        if (isNumber(matched.recipient_width_overlap)) {
          geometry.push(tipRow(
            "Overlap",
            fmtPercent(matched.recipient_width_overlap)
          ));
        }
        if (isNumber(matched.recipient_fraction_ratio)) {
          geometry.push(tipRow(
            "Width ratio", fmt(matched.recipient_fraction_ratio)
          ));
        }
        html += tipSection("Geometry", tipGrid(geometry));
      } else {
        var chromosomeRows = [];
        if (hasValue(chromosome.component)) {
          chromosomeRows.push(tipRow("Component", chromosome.component));
        }
        if (chromosome.deleted === true) {
          chromosomeRows.push(tipRow(
            "State", "deleted recipient (no coverage)"
          ));
        } else if (chromosome.added === true) {
          chromosomeRows.push(tipRow("State", "added donor chromosome"));
        } else {
          chromosomeRows.push(tipRow("State", "recipient chromosome"));
        }
        if (Number(chromosome.copies) > 1) {
          chromosomeRows.push(tipRow(
            "Homolog", chromosome.copy + " of " + chromosome.copies
          ));
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
          tipRow("Position", fmtBp(fraction * (chromosome.length || 0))),
          tipRow(
            "Score",
            index >= 0 ? fmtSigned(scores[index]) : "NA"
          )
        ];
        html += tipSection("Position", tipGrid(position));
      }

      html += '<div class="mamamia-tip-footer">' +
        "Click to filter this group</div></div>";
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

    function render() {
      if (!payload) {
        return;
      }
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
        var empty = document.createElement("p");
        empty.className = "mamamia-karyotype-empty";
        empty.textContent = payload.message ||
          "No chromosomes in this result.";
        root.appendChild(empty);
        return;
      }

      var box = el.getBoundingClientRect();
      var availableWidth = Math.max(width || box.width || 800, 480);
      var availableHeight = Math.max(height || box.height || 520, 420);
      var margin = { top: 58, right: 22, bottom: 34, left: 22 };
      var maxBarHeight = Math.max(
        220,
        availableHeight - margin.top - margin.bottom
      );
      var totalHeight = margin.top + maxBarHeight + margin.bottom;
      var maxLength = 0;
      chromosomes.forEach(function (chromosome) {
        var length = Number(chromosome.length);
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
      chromosomes.forEach(function (chromosome) {
        var copies = Math.max(1, Number(chromosome.copies) || 1);
        if (copies > maxCopies) {
          maxCopies = copies;
        }
        if (copies === 1 || Number(chromosome.copy) === 1) {
          nPairs += 1;
        }
      });
      var slotWidth = (availableWidth - margin.left - margin.right) /
        Math.max(1, nPairs);
      var barWidth = Math.max(
        6,
        Math.min(
          24,
          (slotWidth - Math.max(0, maxCopies - 1) * pairGap -
            chromosomeGap) / maxCopies
        )
      );
      var barsWidth = 0;
      chromosomes.forEach(function (chromosome, index) {
        barsWidth += barWidth;
        if (index < chromosomes.length - 1) {
          var copies = Math.max(1, Number(chromosome.copies) || 1);
          var copy = Math.max(1, Number(chromosome.copy) || 1);
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

      addLegend(svg, margin.left, 14, payload);

      var cursorX = margin.left + (plotWidth - barsWidth) / 2;
      var y = margin.top;
      chromosomes.forEach(function (chromosome, index) {
        var added = chromosome.added === true;
        var deleted = !added && chromosome.deleted === true;
        var barCalls = chromosome.calls || [];
        var copies = Math.max(1, Number(chromosome.copies) || 1);
        var copy = Math.max(1, Number(chromosome.copy) || 1);
        var height = Math.max(
          3, maxBarHeight * (Number(chromosome.length) || 0) / maxLength
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
        if (chromosome.selected === false) {
          chromGroup.setAttribute("opacity", "0.32");
        }
        svg.appendChild(chromGroup);

        var baseFill = deleted ? payload.palette.deleted :
          (added ? payload.palette.added : "#f4f4f4");
        var baseStroke = deleted ? payload.palette.negative :
          (added ? payload.palette.addition : "#666666");
        var base = svgEl("rect", {
          x: x,
          y: y,
          width: barWidth,
          height: height,
          fill: baseFill,
          stroke: baseStroke,
          "stroke-width": deleted ? 1.4 : (added ? 1.8 : 0.7),
          "stroke-dasharray": deleted ? "4 3" : null,
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
          if (payload.mode === "heatmap" && !deleted) {
            drawHeatmapForBar(chromosome, barBox, barCalls);
            var crossCalls = barCalls.filter(function (call) {
              return call.cross_group === true;
            });
            if (crossCalls.length) {
              drawCallsForBar(barBox, crossCalls);
            }
          } else {
            drawCallsForBar(barBox, barCalls);
          }
        }

        if (deleted && copy === 1) {
          chromGroup.appendChild(textEl("deleted", {
            x: pairCenter,
            y: margin.top + 10,
            "text-anchor": "middle",
            "font-size": 7,
            "font-weight": 600,
            fill: payload.palette.negative
          }));
        }
        if (added && copy === 1) {
          chromGroup.appendChild(textEl("added", {
            x: pairCenter,
            y: margin.top + 10,
            "text-anchor": "middle",
            "font-size": 7,
            "font-weight": 600,
            fill: payload.palette.addition
          }));
        }

        if (copy === 1 || copies === 1) {
          chromGroup.appendChild(textEl(
            chromosome.label || chromosome.name || chromosome.component || "",
            {
              x: pairCenter,
              y: margin.top + maxBarHeight + 14,
              "text-anchor": "middle",
              "font-size": 9,
              fill: payload.palette.text
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
