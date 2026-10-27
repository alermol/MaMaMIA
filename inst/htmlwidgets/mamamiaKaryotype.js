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

    function drawHeatmap(group, bar) {
      var scores = group.scores || [];
      if (!scores.length) {
        return;
      }
      var bins = Math.min(scores.length, Math.max(40, Math.floor(bar.height / 2)));
      var binHeight = bar.height / bins;
      for (var i = 0; i < bins; i += 1) {
        var from = Math.floor(i * scores.length / bins);
        var to = Math.max(from + 1, Math.floor((i + 1) * scores.length / bins));
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

    function drawCalls(group, bar) {
      var calls = group.calls || [];
      calls.forEach(function (call) {
        var y0 = bar.y + Math.max(0, Math.min(1, call.start)) * bar.height;
        var y1 = bar.y + Math.max(0, Math.min(1, call.end)) * bar.height;
        var colour = payload.palette[call.class] || payload.palette.other;
        var opacity = call.confidence === "high" ? 0.95 :
          (call.confidence === "moderate" ? 0.72 : 0.45);
        if (call.pass === false) {
          opacity *= 0.45;
        }
        bar.group.appendChild(svgEl("rect", {
          x: bar.x,
          y: y0,
          width: bar.width,
          height: Math.max(1, y1 - y0),
          fill: colour,
          opacity: opacity,
          stroke: "#ffffff",
          "stroke-width": 0.5,
          "pointer-events": "none"
        }));
      });
    }

    function describe(group, chromosome, added, deleted, fraction) {
      var name = escapeHtml(chromosome.name || chromosome.subgenome || "");
      var mb = fraction * (chromosome.length || group.length || 0) / 1e6;
      var lines = [
        "<strong>" + name + "</strong>",
        escapeHtml(chromosome.subgenome || "")
      ];
      if (deleted) {
        lines.push("recipient chromosome deleted (no coverage)");
        return lines.join("<br>");
      }
      if (added) {
        lines.push("added donor chromosome");
      }
      if (payload.mode === "calls") {
        var matched = null;
        (group.calls || []).forEach(function (call) {
          if (fraction >= call.start && fraction <= call.end) {
            matched = call;
          }
        });
        if (matched) {
          lines.push(
            escapeHtml(matched.class) + " · " +
            escapeHtml(matched.confidence)
          );
          lines.push("score " + fmt(matched.score));
        } else {
          lines.push("no call at this position");
        }
      } else {
        var scores = group.scores || [];
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

      var groups = payload.groups || [];
      if (!groups.length) {
        var empty = document.createElement("p");
        empty.className = "mamamia-karyotype-empty";
        empty.textContent = payload.message ||
          "No chromosome groups in this result.";
        root.appendChild(empty);
        return;
      }

      var box = el.getBoundingClientRect();
      var availableWidth = Math.max(width || box.width || 800, 480);
      var availableHeight = Math.max(height || box.height || 520, 420);
      var margin = { top: 58, right: 22, bottom: 34, left: 22 };
      var columnWidth = Math.max(
        96,
        (availableWidth - margin.left - margin.right) / groups.length
      );
      var totalWidth = Math.max(
        availableWidth,
        margin.left + groups.length * columnWidth + margin.right
      );
      var barHeight = Math.max(
        260,
        availableHeight - margin.top - margin.bottom
      );
      var totalHeight = margin.top + barHeight + margin.bottom;
      var svg = svgEl("svg", {
        width: totalWidth,
        height: totalHeight,
        viewBox: "0 0 " + totalWidth + " " + totalHeight,
        role: "img",
        "aria-label": "Recipient karyotype"
      });
      svg.style.display = "block";
      root.appendChild(svg);

      addLegend(svg, margin.left, 14, payload);

      groups.forEach(function (group, groupIndex) {
        var columnX = margin.left + groupIndex * columnWidth;
        var column = svgEl("g", {
          "class": "karyotype-group",
          "data-group": group.group
        });
        if (!group.selected) {
          column.setAttribute("opacity", "0.32");
        }
        svg.appendChild(column);

        column.appendChild(svgEl("line", {
          x1: columnX, x2: columnX,
          y1: margin.top - 8, y2: margin.top + barHeight + 18,
          stroke: payload.palette.grid,
          "stroke-width": 1,
          "pointer-events": "none"
        }));
        column.appendChild(textEl(group.label || ("Group " + group.group), {
          x: columnX + columnWidth / 2,
          y: 38,
          "text-anchor": "middle",
          "font-size": 12,
          "font-weight": 600,
          fill: payload.palette.text
        }));

        var bars = [];
        (group.recipient || []).forEach(function (chromosome) {
          bars.push({ kind: "recipient", chromosome: chromosome });
        });
        if ((group.added || []).length) {
          bars.push({ kind: "gap" });
          (group.added || []).forEach(function (chromosome) {
            bars.push({ kind: "added", chromosome: chromosome });
          });
        }
        var gap = 5;
        var gapWidth = 10;
        var barWidth = Math.max(
          9,
          Math.min(
            26,
            (columnWidth - 20) / Math.max(1, bars.length) - gap
          )
        );
        var widthWithoutGaps = 0;
        bars.forEach(function (bar) {
          widthWithoutGaps += bar.kind === "gap" ? gapWidth : barWidth + gap;
        });
        var cursorX = columnX + (columnWidth - widthWithoutGaps + gap) / 2;

        bars.forEach(function (bar) {
          if (bar.kind === "gap") {
            cursorX += gapWidth;
            return;
          }
          var chromosome = bar.chromosome;
          var added = bar.kind === "added";
          var deleted = !added && chromosome.deleted === true;
          var x = cursorX;
          cursorX += barWidth + gap;
          var y = margin.top;

          var baseFill = deleted ? payload.palette.deleted :
            (added ? payload.palette.added : "#f4f4f4");
          var baseStroke = deleted ? payload.palette.negative :
            (added ? payload.palette.addition : "#666666");
          var base = svgEl("rect", {
            x: x,
            y: y,
            width: barWidth,
            height: barHeight,
            fill: baseFill,
            stroke: baseStroke,
            "stroke-width": deleted ? 1.4 : (added ? 1.8 : 0.7),
            "stroke-dasharray": deleted ? "4 3" : null,
            rx: 2
          });
          column.appendChild(base);

          var barBox = {
            group: column,
            x: x,
            y: y,
            width: barWidth,
            height: barHeight
          };
          if (!deleted) {
            if (payload.mode === "heatmap") {
              drawHeatmap(group, barBox);
            } else {
              drawCalls(group, barBox);
            }
          }

          if (deleted) {
            column.appendChild(textEl("deleted", {
              x: x + barWidth / 2,
              y: 49,
              "text-anchor": "middle",
              "font-size": 8,
              "font-weight": 600,
              fill: payload.palette.negative
            }));
          }
          if (added) {
            column.appendChild(textEl("added", {
              x: x + barWidth / 2,
              y: 49,
              "text-anchor": "middle",
              "font-size": 8,
              "font-weight": 600,
              fill: payload.palette.addition
            }));
          }

          var subgenome = chromosome.subgenome || "";
          if (added && chromosome.copies > 1) {
            subgenome += " " + chromosome.copy;
          }
          column.appendChild(textEl(subgenome, {
            x: x + barWidth / 2,
            y: margin.top + barHeight + 14,
            "text-anchor": "middle",
            "font-size": 9,
            fill: payload.palette.text
          }));

          var hit = svgEl("rect", {
            x: x - 2,
            y: y - 2,
            width: barWidth + 4,
            height: barHeight + 4,
            fill: "transparent",
            "class": "karyotype-hit"
          });
          hit.style.cursor = "pointer";
          column.appendChild(hit);
          hit.addEventListener("mousemove", function (event) {
            var rect = svg.getBoundingClientRect();
            var pointX = event.clientX - rect.left;
            var pointY = event.clientY - rect.top;
            var fraction = Math.max(
              0, Math.min(1, (pointY - y) / barHeight)
            );
            tooltip.innerHTML = describe(
              group, chromosome, added, deleted, fraction
            );
            tooltip.style.display = "block";
            tooltip.style.left = Math.min(
              Math.max(4, pointX + 12),
              Math.max(4, root.clientWidth - 230)
            ) + "px";
            tooltip.style.top = Math.max(4, pointY + 12) + "px";
          });
          hit.addEventListener("mouseleave", function () {
            tooltip.style.display = "none";
          });
          hit.addEventListener("click", function (event) {
            event.stopPropagation();
            sendClick(group.group);
          });
        });
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
