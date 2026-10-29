HTMLWidgets.widget({
  name: "mamamiaKaryotype",
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
      if (value === null || value === undefined || !isFinite(value)) {
        return "#eeeeee";
      }
      var t = Math.max(-1, Math.min(1, value / limit));
      if (t < 0) {
        return mixColour(palette.zero, palette.negative, -t);
      }
      return mixColour(palette.zero, palette.positive, t);
    }

    function fmt(value) {
      return isFinite(value) ? Number(value).toFixed(2) : "NA";
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

    function describe(chromosome, fraction, calls) {
      var name = escapeHtml(
        chromosome.label || chromosome.name || chromosome.component || ""
      );
      var mb = fraction * (chromosome.length || 0) / 1e6;
      var lines = [
        "<strong>" + name + "</strong>",
        escapeHtml(chromosome.component || "")
      ];
      if (Number(chromosome.copies) > 1) {
        lines.push("homolog " + chromosome.copy + " of " +
          chromosome.copies);
      }
      if (chromosome.deleted === true) {
        lines.push("recipient chromosome deleted (no coverage)");
      }
      lines.push(chromosome.added === true
        ? "added donor chromosome"
        : "recipient chromosome");
      var matched = null;
      (calls || []).forEach(function (call) {
        if (fraction >= call.start && fraction <= call.end) {
          matched = call;
        }
      });
      if (matched) {
        lines.push(
          escapeHtml(matched.class) + " · " +
          escapeHtml(matched.confidence)
        );
        if (hasValue(matched.recipient_assignment) &&
            matched.recipient_assignment !== "assigned") {
          lines.push("assignment: " +
            escapeHtml(matched.recipient_assignment));
        }
        if (isFinite(matched.recipient_weight)) {
          lines.push("weight " + fmt(matched.recipient_weight));
        }
        if (hasValue(matched.recipient_chr_name)) {
          lines.push("assigned to " + escapeHtml(matched.recipient_chr_name));
        }
        if (isFinite(matched.affected_copies)) {
          var total = isFinite(matched.recipient_copy_number)
            ? " of " + matched.recipient_copy_number
            : "";
          lines.push("affected copies " + matched.affected_copies + total);
        } else if (isFinite(matched.affected_fraction)) {
          lines.push("recipient loss " +
            Math.round(100 * matched.affected_fraction) + "%");
        }
        if (hasValue(matched.dosage_confidence) &&
            matched.dosage_confidence !== "unavailable") {
          lines.push("dosage " + escapeHtml(matched.dosage_confidence));
        }
        if (matched.cross_group === true) {
          var donor = hasValue(matched.donor_chr_name)
            ? matched.donor_chr_name
            : ("group " + matched.donor_group);
          lines.push("cross-group: donor " + escapeHtml(donor) +
            " → recipient group " +
            escapeHtml(matched.recipient_group || ""));
        }
        if (hasValue(matched.donor_chr_name)) {
          lines.push("donor " + escapeHtml(matched.donor_chr_name));
        }
        lines.push("score " + fmt(matched.score));
      } else if (payload.mode === "calls") {
        lines.push("no assigned donor segment at this position");
      }
      if (payload.mode === "heatmap" && matched &&
          matched.cross_group !== true) {
        var scores = chromosome.scores || [];
        var index = scores.length
          ? Math.max(0, Math.min(scores.length - 1,
            Math.floor(fraction * scores.length)))
          : -1;
        lines.push(mb.toFixed(1) + " Mb");
        lines.push("score " + fmt(index >= 0 ? scores[index] : NaN));
      }
      return lines.join("<br>");
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
