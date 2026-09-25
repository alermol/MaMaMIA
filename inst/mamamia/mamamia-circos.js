// Interactive layer for the MaMaMIA introgression circos.
//
// BioCircos gives straight tangential genome labels, a label flip test that can
// never fire (d.angle > 2*pi && d.angle < 0), no handlers on the genome ring,
// no bridge to Shiny, and tooltips that stay under the pointer. This file adds
// those behaviours and nothing else; every geometry it needs is read from the
// rendered SVG, so it survives any widget size.
//
// Selection is applied here rather than by re-rendering the widget, because a
// rebuild serialises thousands of arcs. A picked call simply turns from red to
// blue, on its ribbon and on its two chromosome marks. The server stays
// authoritative: it echoes the table selection back as "mamamia-selection", and
// a click on a band reports the call it hit.
//
// Entry point, called from htmlwidgets::onRender():
//   window.mamamiaCircosEnhance(el, x);

window.mamamiaCircosInstances = window.mamamiaCircosInstances || [];

window.mamamiaCircosPaint = function(instance, wanted) {
  var svg = d3.select(instance.el).select('svg');
  var ribbonRows = instance.x.ribbon_rows || [];
  var spanRows = instance.x.span_rows || [];
  svg.selectAll(instance.ribbonSelector).each(function(d, i) {
    var on = wanted.indexOf(Number(ribbonRows[i])) >= 0;
    d3.select(this).style('stroke', on ? instance.pickedColor : null);
  });
  svg.selectAll(instance.spanSelector).each(function(d, i) {
    var on = wanted.indexOf(Number(spanRows[i])) >= 0;
    d3.select(this).style('fill', on ? instance.pickedColor : null);
  });
  instance.selected = wanted;
};

window.mamamiaCircosApplySelection = function(ids) {
  var wanted = (ids || []).map(Number);
  window.mamamiaCircosInstances.forEach(function(instance) {
    if (document.body.contains(instance.el)) {
      window.mamamiaCircosPaint(instance, wanted);
    }
  });
};

window.mamamiaCircosEnhance = function(el, x) {
  var summary = x.chromosome_summary || {};
  var donors = x.donor_chr_ids || [];
  var ids = x.chromosome_ids || [];
  var rows = x.hit_area_rows || [];
  var svg = d3.select(el).select('svg');
  var svgNode = svg.node();

  var instance = {
    el: el,
    x: x,
    selected: [],
    pickedColor: '#2166ac',
    ribbonSelector: 'path.BioCircosLINK[stroke="#cb181d55"]',
    spanSelector: 'path.BioCircosCNV[fill="#cb181d"]',
    hitSelector: 'path.BioCircosLINK[stroke="rgba(255,255,255,0.01)"]'
  };

  // BioCircos renders once at whatever size the container had at that moment.
  // If the container has grown since, the SVG keeps its old coordinate system
  // while being displayed stretched, so the pointer no longer lines up with the
  // drawing. Ask for a single re-render at the size actually laid out.
  if (svgNode && !el.getAttribute('data-mamamia-sized')) {
    var box = el.getBoundingClientRect();
    var declared = parseFloat(svgNode.getAttribute('width')) || 0;
    if (box.width > 20 && declared > 0 && Math.abs(declared - box.width) > 2) {
      el.setAttribute('data-mamamia-sized', '1');
      setTimeout(function() { window.dispatchEvent(new Event('resize')); }, 200);
    }
  }

  // No text in the widget should read or behave as selectable text. BioCircos
  // keeps thousands of hidden dragText nodes; the stylesheet covers this too,
  // but a widget can render before its stylesheet is applied.
  svg.style('user-select', 'none').style('-webkit-user-select', 'none');
  svg.selectAll('text')
    .style('pointer-events', 'none')
    .style('user-select', 'none');
  d3.selectAll('[id^="BioCircos"][id$="Tooltip"]')
    .style('pointer-events', 'none')
    .style('user-select', 'none');

  var arcs = svg.selectAll('path[name]')[0];
  var defs = svg.append('defs');

  // Only the ideogram labels hold a chromosome name; they are not siblings of
  // the bands, so they are matched on content.
  var texts = svg.selectAll('text').filter(function() {
    return ids.indexOf(this.textContent) >= 0;
  });

  var names = [];
  texts.each(function(d, i) {
    var t = d3.select(this);
    var name = t.text();
    names.push(name);

    // Readable against either band: dark ink on the light recipient, white on
    // the dark donor.
    t.style('fill', donors.indexOf(name) >= 0 ? '#ffffff' : '#2f3437');

    var node = arcs[i];
    if (!node) { return; }

    var transform = t.attr('transform') || '';
    var open = transform.indexOf('rotate(');
    var close = transform.indexOf(')', open);
    var deg = open >= 0 ? parseFloat(transform.substring(open + 7, close)) : 0;

    // Sampling the sector outline yields the band radius and angular span at
    // any widget size, unlike genomeLabelDy, which is a fixed pixel offset.
    var len = node.getTotalLength();
    var rmin = Infinity, rmax = 0, base = null, amin = Infinity, amax = -Infinity;
    for (var s = 0; s <= 32; s++) {
      var pt = node.getPointAtLength(len * s / 32);
      var rr = Math.sqrt(pt.x * pt.x + pt.y * pt.y);
      if (rr < rmin) { rmin = rr; }
      if (rr > rmax) { rmax = rr; }
      var ang = Math.atan2(pt.y, pt.x);
      if (base === null) { base = ang; }
      var rel = ang - base;
      while (rel > Math.PI) { rel = rel - 2 * Math.PI; }
      while (rel < -Math.PI) { rel = rel + 2 * Math.PI; }
      if (rel < amin) { amin = rel; }
      if (rel > amax) { amax = rel; }
    }
    var radius = (rmin + rmax) / 2;
    if (!(radius > 0) || !(amax > amin)) { return; }

    // Centre the name on its own segment, bent along the band.
    var mid = base + (amin + amax) / 2 + Math.PI / 2;
    var half = (amax - amin) / 2 * 0.92;
    var lower = deg > 90 && deg < 270;
    var a0 = lower ? mid + half : mid - half;
    var a1 = lower ? mid - half : mid + half;
    var id = 'mamamia-arc-' + i;
    defs.append('path').attr('id', id).attr('fill', 'none').attr('d',
      'M' + (Math.sin(a0) * radius) + ',' + (-Math.cos(a0) * radius) +
      ' A' + radius + ',' + radius + ' 0 0 ' + (lower ? 0 : 1) + ' ' +
      (Math.sin(a1) * radius) + ',' + (-Math.cos(a1) * radius));

    t.attr('transform', null);
    var path = t.text('').append('textPath');
    path.node().setAttributeNS('http://www.w3.org/1999/xlink', 'xlink:href', '#' + id);
    path.attr('startOffset', '50%').text(name);

    // Names are never cut off: an over-long one shrinks to its arc instead.
    var arcLength = radius * 2 * half;
    var textLength = 0;
    try { textLength = path.node().getComputedTextLength(); } catch (e) { textLength = 0; }
    if (textLength > 0 && arcLength > 0 && textLength > arcLength) {
      var size = parseFloat(window.getComputedStyle(t.node()).fontSize) || 6;
      t.style('font-size', (size * arcLength / textLength).toFixed(2) + 'px');
    }
  });

  // Hover on a chromosome: the calls that touch it, with locations and sizes.
  var tip = d3.select('body').append('div')
    .style('position', 'absolute')
    .style('pointer-events', 'none')
    .style('user-select', 'none')
    .style('opacity', 0)
    .style('background', 'rgba(20,24,28,0.94)')
    .style('color', '#fff')
    .style('padding', '6px 8px')
    .style('border-radius', '4px')
    .style('font-size', '11px')
    .style('line-height', '1.45')
    .style('max-width', '360px')
    .style('z-index', '10000');
  svg.selectAll('path[name]')
    .on('mouseover', function() {
      var idx = parseInt(d3.select(this).attr('name'), 10) - 1;
      tip.html(summary[names[idx]] || names[idx] || '').style('opacity', 1);
    })
    .on('mousemove', function() {
      tip.style('left', (d3.event.pageX + 14) + 'px')
         .style('top', (d3.event.pageY + 14) + 'px');
    })
    .on('mouseout', function() { tip.style('opacity', 0); });

  // Click on a band: repaint at once and report the call, so the table selection
  // follows without the widget being rebuilt.
  var hits = svg.selectAll(instance.hitSelector);
  hits.on('click', function() {
    var idx = Array.prototype.indexOf.call(hits[0], this);
    if (idx < 0 || idx >= rows.length) { return; }
    var row = Number(rows[idx]);
    var wanted = instance.selected.slice();
    var at = wanted.indexOf(row);
    if (at >= 0) { wanted.splice(at, 1); } else { wanted.push(row); }
    window.mamamiaCircosPaint(instance, wanted);
    if (window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue('circos_clicked_row', row, {priority: 'event'});
    }
  });

  // Remember this instance so selection messages reach it without a rebuild.
  var known = window.mamamiaCircosInstances.filter(function(other) {
    return other.el === el;
  });
  if (!known.length) {
    window.mamamiaCircosInstances.push(instance);
  }
  window.mamamiaCircosPaint(instance, instance.selected);

  if (window.Shiny && Shiny.addCustomMessageHandler && !window.mamamiaSelectionBound) {
    window.mamamiaSelectionBound = true;
    Shiny.addCustomMessageHandler('mamamia-selection', function(selection) {
      window.mamamiaCircosApplySelection(selection);
    });
  }
};
