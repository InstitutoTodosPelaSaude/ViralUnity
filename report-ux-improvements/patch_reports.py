#!/usr/bin/env python3
"""Apply the agreed visual-QA fixes to all four ViralUnity report HTMLs.

The four files share a byte-identical template (CSS block + custom JS tail), so
one patch applies to all. Every change is client-side; the pre-baked Plotly JSON
is never rewritten in place. Two charts are rebuilt client-side from the embedded
STATS data (average depth) or adjusted via the Plotly API (reads, aggregated
coverage) at runtime.
"""
import sys, pathlib

FILES = ["real_illumina.html", "real_nanopore.html", "real_influenza.html", "real_guaroa.html"]
ROOT = pathlib.Path("/Users/filiperomero/Desktop/report_previews")

# --- 1) CSS: anchor on the unique .chart-note rule, append new styles ------------
CSS_ANCHOR = ".chart-note { font-size: 0.8rem; color: var(--text-muted); margin: 8px 4px 0; }"
CSS_ADD = CSS_ANCHOR + """
.scale-toggle { display: inline-flex; border: 1px solid var(--border); border-radius: 7px; overflow: hidden; vertical-align: middle; margin: 0 0 14px; }
.scale-toggle button { border: none; background: var(--card-bg); color: var(--text-secondary); font: inherit; font-size: 0.8rem; padding: 6px 12px; cursor: pointer; }
.scale-toggle button + button { border-left: 1px solid var(--border); }
.scale-toggle button.active { background: var(--accent); color: #fff; }
.scale-toggle button:hover:not(.active) { background: var(--hover); }
.scale-toggle button:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
.chart-controls { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 16px; margin-bottom: 6px; }
.chart-controls .scale-toggle { margin: 0; }
table.stats-table th[role="button"]:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
/* Separate the per-sample stats table from the coverage plot(s) below it. The
   table already carries a hairline under its last row, so spacing alone reads
   cleanly here (no second line). */
#sample-coverage, #segment-coverage { margin-top: 28px; }
/* Stacked coverage plots (one per segment, or one per sample in the per-segment
   view) get breathing room and a hairline divider so each reads as its own panel. */
#sample-coverage > div + div, #segment-coverage > div + div { margin-top: 32px; padding-top: 28px; border-top: 1px solid var(--border); }
/* Per-sample detail is collapsed by default so run-level info leads; only those
   who want per-sample coverage expand it. */
.sample-details > summary { cursor: pointer; list-style: none; display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.sample-details > summary::-webkit-details-marker { display: none; }
.sample-details > summary::before { content: '\\25B8'; color: var(--text-muted); font-size: 0.85em; transition: transform 0.15s ease; }
.sample-details[open] > summary::before { transform: rotate(90deg); }
.sample-details > summary:focus-visible { outline: 2px solid var(--accent); outline-offset: 3px; border-radius: 4px; }
.sample-summary-title { font-size: 1rem; font-weight: 600; }
.sample-summary-hint { font-size: 0.8rem; color: var(--text-muted); font-weight: 400; }
.sample-details-body { margin-top: 18px; }
@media (prefers-reduced-motion: reduce) { .sample-details > summary::before { transition: none; } }"""

# --- 2) JS: exact current block (from the ORIGINAL file) -------------------------
JS_OLD = '''/* ---- lazy per-sample coverage plots ---- */
// Log10 range for a depth axis — must match Python's _log_depth_range: floor at
// 1 (10^0), top = 1.5x the data max but never below 150, so 20x/100x always show.
function logDepthRange(maxDepth) {
  const top = Math.max(150, (maxDepth || 0) * 1.5);
  return [0, Math.log10(top)];
}
function coverageLayout(title, maxDepth) {
  const dark = isDark();
  return {
    width: FIG_WIDTH, height: FIG_HEIGHT, title: {text: title},
    margin: {l: 60, r: 30, t: 50, b: 60},
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: {family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 13,
           color: dark ? '#c3c2b7' : '#52514e'},
    xaxis: {title: {text: 'Genome position'}, gridcolor: 'rgba(137,135,129,0.25)', zeroline: false},
    yaxis: {title: {text: 'Depth (log)'}, type: 'log', range: logDepthRange(maxDepth),
            gridcolor: 'rgba(137,135,129,0.25)', zeroline: false},
    hovermode: 'x unified', showlegend: false,
    shapes: [
      {type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 20, y1: 20, line: {color: '#898781', width: 1, dash: 'dot'}},
      {type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 100, y1: 100, line: {color: '#898781', width: 1, dash: 'dot'}}
    ],
    annotations: [
      {xref: 'paper', x: 1, y: 20, yref: 'y', text: '20x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}},
      {xref: 'paper', x: 1, y: 100, yref: 'y', text: '100x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}}
    ]
  };
}
function showSample(sample) {
  const entries = COVERAGE[sample] || [];
  const holder = document.getElementById('sample-coverage');
  holder.innerHTML = '';
  entries.forEach(function (e, i) {
    const div = document.createElement('div');
    div.id = 'sample-cov-' + i;
    holder.appendChild(div);
    const label = e.label ? (sample + ' — ' + e.label) : sample;
    const maxDepth = e.y.length ? Math.max.apply(null, e.y) : 0;
    // Clamp zeros to 1 so the log axis can render them (shown at the floor).
    const y = e.y.map(function (v) { return Math.max(1, v); });
    const trace = {x: e.x, y: y, mode: 'lines', line: {width: 2, color: '#2a78d6'}, name: label};
    Plotly.newPlot(div, [trace], coverageLayout(label, maxDepth), {responsive: false, displayModeBar: false});
  });
  renderSampleStats(sample);
}
function renderSampleStats(sample) {
  const rows = STATS[sample] || [];
  const box = document.getElementById('sample-stats');
  if (!rows.length) { box.innerHTML = ''; return; }
  const cols = Object.keys(rows[0]);
  let html = '<table class="stats-table"><thead><tr>';
  cols.forEach(function (c) { html += '<th>' + c + '</th>'; });
  html += '</tr></thead><tbody>';
  rows.forEach(function (r) {
    html += '<tr>';
    cols.forEach(function (c) { html += '<td>' + r[c] + '</td>'; });
    html += '</tr>';
  });
  html += '</tbody></table>';
  box.innerHTML = html;
}

/* ---- init ---- */
window.addEventListener('load', function () {
  themeCharts();
  const first = document.getElementById('sampleSelect');
  if (first && first.value) showSample(first.value);
});
'''

JS_NEW = '''/* ---- number formatting (shared) ---- */
// Group digits with thousands separators. Pinned to en-US so the separator is a
// comma on every viewer's locale (a bare toLocaleString() renders "46.341" in
// de-DE etc., which reads as a decimal). Only the visible text is grouped; the
// data-sort attribute keeps the raw value so column sorting is unaffected.
var GROUP = new Intl.NumberFormat('en-US');
function fmtNumberString(text) {
  var t = String(text).trim();
  if (!/^-?\\d+(\\.\\d+)?$/.test(t)) return text;   // leave %, ids, labels alone
  var neg = t.charAt(0) === '-';
  var body = neg ? t.slice(1) : t;
  var parts = body.split('.');
  var grouped = GROUP.format(parseInt(parts[0], 10));
  if (parts.length > 1) grouped += '.' + parts[1];   // keep fraction digits verbatim
  return (neg ? '-' : '') + grouped;
}
function fmtPct(p) {
  if (!isFinite(p)) return '\\u2014';
  if (p >= 1)    return p.toFixed(1) + '%';
  if (p >= 0.01) return p.toFixed(2) + '%';
  if (p > 0)     return p.toFixed(3) + '%';
  return '0%';
}
// Mapping rate = mapped reads / QC-passed reads, made unit-consistent by library
// layout. ViralUnity's Illumina path is paired-end and runs a QC filter (so
// Total > QC-passed) while counting Mapped as individual reads; its Nanopore
// path is single-end and runs no filter (Total == QC-passed). We therefore treat
// "Total > QC-passed" as the paired-end signal and double the QC-passed
// denominator so it is in the same unit (individual reads) as Mapped.
// NOTE: this is a heuristic standing in for a real library-layout field - see
// RECOMMENDATIONS.md. It is correct for the pipeline's current behaviour but can
// misclassify an unfiltered paired run or a filtered single-end run.
function mappedPct(total, qc, mapped) {
  if (!(isFinite(total) && isFinite(qc) && isFinite(mapped)) || qc <= 0) return null;
  var paired = total > qc;
  return {pct: mapped / (paired ? qc * 2 : qc) * 100, paired: paired};
}

/* ---- main statistics table: number formatting + mapped-as-percentage ---- */
function statsColMap() {
  var map = {};
  document.querySelectorAll('table#stats-table thead th').forEach(function (th) {
    map[th.textContent.replace(/[\\u25B2\\u25BC]/g, '').trim()] = parseInt(th.getAttribute('data-col'), 10);
  });
  return map;
}
function formatStatsTable() {
  var table = document.getElementById('stats-table');
  if (!table || !table.tBodies.length) return;
  var col = statsColMap();
  var mCol = col['Mapped reads'], tCol = col['Total reads'], qCol = col['QC-passed reads'];
  Array.from(table.tBodies[0].rows).forEach(function (tr) {
    if (mCol != null && tCol != null && qCol != null) {
      var total = parseFloat(tr.cells[tCol].getAttribute('data-sort'));
      var qc = parseFloat(tr.cells[qCol].getAttribute('data-sort'));
      var mapped = parseFloat(tr.cells[mCol].getAttribute('data-sort'));
      var r = mappedPct(total, qc, mapped);
      if (r) {
        var cell = tr.cells[mCol];
        cell.textContent = fmtPct(r.pct);
        cell.setAttribute('data-sort', String(r.pct));   // sort by rate, matching display
        cell.setAttribute('title', GROUP.format(mapped) + ' mapped reads / ' +
          (r.paired ? '2\\u00d7 ' : '') + 'QC-passed reads');
      }
    }
    Array.from(tr.cells).forEach(function (td, idx) {
      if (idx === mCol) return;                          // already replaced with %
      td.textContent = fmtNumberString(td.textContent);
    });
  });
  if (mCol != null) {                                    // relabel header, keep the sort arrow + onclick
    var th = table.querySelector('thead th[data-col="' + mCol + '"]');
    if (th) {
      var arrow = th.querySelector('.sort-arrow');
      th.textContent = 'Mapped %';
      if (arrow) th.appendChild(arrow);
    }
  }
}

/* ---- sortable-header accessibility ---- */
function enhanceSortHeaders() {
  document.querySelectorAll('table.stats-table th[onclick]').forEach(function (th) {
    th.setAttribute('role', 'button');
    th.setAttribute('tabindex', '0');
    if (!th.hasAttribute('aria-sort')) th.setAttribute('aria-sort', 'none');
    th.addEventListener('keydown', function (ev) {
      if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); th.click(); }
    });
    th.addEventListener('click', function () {           // runs after the inline onclick
      th.closest('table').querySelectorAll('th[aria-sort]').forEach(function (o) {
        o.setAttribute('aria-sort', 'none');
      });
      th.setAttribute('aria-sort', th.getAttribute('data-asc') === 'true' ? 'ascending' : 'descending');
    });
  });
}

/* ---- Linear/Log10 toggle control ---- */
function scaleToggleEl(initial, onChange) {
  var wrap = document.createElement('div');
  wrap.className = 'scale-toggle';
  wrap.setAttribute('role', 'group');
  wrap.setAttribute('aria-label', 'Y-axis scale');
  ['linear', 'log'].forEach(function (s) {
    var b = document.createElement('button');
    b.type = 'button';
    b.textContent = s === 'linear' ? 'Linear' : 'Log10';
    b.setAttribute('data-scale', s);
    if (s === initial) b.className = 'active';
    b.addEventListener('click', function () {
      wrap.querySelectorAll('button').forEach(function (x) { x.classList.remove('active'); });
      b.classList.add('active');
      onChange(s);
    });
    wrap.appendChild(b);
  });
  return wrap;
}
function themeFont() {
  return {family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 13,
          color: isDark() ? '#c3c2b7' : '#52514e'};
}
function maxOfFullData(gd) {
  var mx = 0;
  (gd._fullData || []).forEach(function (t) {
    var y = t.y;
    if (y && y.length) for (var i = 0; i < y.length; i++) { if (y[i] > mx) mx = y[i]; }
  });
  return mx;
}
// Decoded y-values for a trace (Plotly's binary bdata form lives decoded in _fullData).
function decodedY(gd, idx) {
  var fd = gd._fullData && gd._fullData[idx];
  if (fd && fd.y && fd.y.length) return Array.prototype.slice.call(fd.y);
  var d = gd.data && gd.data[idx];
  return (d && d.y && d.y.length) ? Array.prototype.slice.call(d.y) : [];
}
// Clean layout for the aggregated-coverage figure, rebuilt via Plotly.react on
// each scale toggle. This also drops the embedded title + legend and pins the
// threshold markers to 20/100 (no dependence on the baked chart's buggy state).
function aggregatedLayout(scale, mx) {
  var isLog = scale === 'log';
  return {
    autosize: true, height: FIG_HEIGHT,
    margin: {l: 60, r: 30, t: 24, b: 60},
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: themeFont(),
    xaxis: {title: {text: 'Genome position'}, gridcolor: 'rgba(137,135,129,0.25)', zeroline: false},
    yaxis: isLog
      ? {title: {text: 'Depth'}, type: 'log', range: [0, Math.log10(Math.max(150, mx * 1.5))], dtick: 1,
         gridcolor: 'rgba(137,135,129,0.25)', zeroline: false}
      : {title: {text: 'Depth'}, type: 'linear', range: [0, mx * 1.08],
         gridcolor: 'rgba(137,135,129,0.25)', zeroline: false},
    hovermode: 'x unified', showlegend: false,
    shapes: [
      {type: 'line', xref: 'x domain', x0: 0, x1: 1, y0: 20, y1: 20, line: {color: '#898781', width: 1, dash: 'dot'}},
      {type: 'line', xref: 'x domain', x0: 0, x1: 1, y0: 100, y1: 100, line: {color: '#898781', width: 1, dash: 'dot'}}
    ],
    annotations: [
      {xref: 'paper', x: 1, y: 20, yref: 'y', text: '20x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}},
      {xref: 'paper', x: 1, y: 100, yref: 'y', text: '100x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}}
    ]
  };
}

/* ---- pre-baked Plotly charts: adjust at runtime via the Plotly API ---- */
function enhanceBakedCharts() {
  if (!window.Plotly) return;
  document.querySelectorAll('main .card').forEach(function (card) {
    var h2 = card.querySelector('h2');
    if (!h2) return;
    var title = h2.textContent.trim();
    var gds = card.querySelectorAll('.plotly-graph-div');
    if (!gds.length) return;

    if (/^Average depth/i.test(title)) {
      card.remove();                      // per-sample mean depth already lives in the stats tables
      return;
    }

    // Reflow: neutralize the hard-coded 900px wrapper so charts shrink to the card.
    gds.forEach(function (gd) {
      var wrap = gd.closest('[style*="width"]');
      if (wrap) { wrap.style.width = '100%'; wrap.style.maxWidth = '900px'; }
      try { Plotly.relayout(gd, {autosize: true, width: null}); Plotly.Plots.resize(gd); } catch (e) {}
    });

    try {
      if (/^Reads per sample/i.test(title)) {
        h2.textContent = 'Sequencing throughput';
        gds.forEach(function (gd) {
          var data = gd.data || [];
          var totalIdx = -1, qcIdx = -1, mapIdx = -1;
          data.forEach(function (t, idx) {
            var name = t.name || '';
            var c = /total/i.test(name) ? '#2a78d6'
                  : /qc/i.test(name)    ? '#1baf7a'
                  :                       '#eda100';
            Plotly.restyle(gd, {'marker.color': c, 'opacity': 1}, [idx]);
            if (/total/i.test(name)) totalIdx = idx;
            else if (/qc/i.test(name)) qcIdx = idx;
            else if (/mapped/i.test(name)) mapIdx = idx;
          });
          // Show the lower (Mapped) panel as % mapped, matching the table, instead
          // of an absolute count on a second read-count scale.
          if (totalIdx >= 0 && qcIdx >= 0 && mapIdx >= 0) {
            var totY = decodedY(gd, totalIdx), qcY = decodedY(gd, qcIdx), mapY = decodedY(gd, mapIdx);
            var pct = mapY.map(function (m, i) { var r = mappedPct(totY[i], qcY[i], m); return r ? r.pct : null; });
            Plotly.restyle(gd, {y: [pct], name: ['Mapped %'],
              hovertemplate: ['%{x}<br>%{y:.2f}% mapped<extra></extra>']}, [mapIdx]);
            // Widen the gap between the two panels (0.14 -> 0.22 of the plot height)
            // and re-anchor the lower subplot title to the new domain edge.
            Plotly.relayout(gd, {
              'yaxis2.title.text': 'Mapped %', 'yaxis2.range': [0, 100], 'yaxis2.autorange': false,
              'yaxis.domain': [0.60, 1.0], 'yaxis2.domain': [0.0, 0.38]
            });
            var anns = (gd.layout.annotations || []).map(function (a) {
              return /mapped/i.test(a.text || '') ? Object.assign({}, a, {text: 'Mapped %', y: 0.38}) : a;
            });
            Plotly.relayout(gd, {annotations: anns});
          }
          Plotly.relayout(gd, {showlegend: true});
        });
      } else if (/^Aggregated coverage/i.test(title)) {
        h2.textContent = title.replace(/\\s*\\(log scale\\)/i, '');
        // Rebuild with Plotly.react (not relayout) on toggle. A relayout that flips
        // yaxis.type on these ~2000-point line figures hits a Plotly slow path
        // (~6 s, freezing the tab); react re-renders the same data in ~15 ms. It
        // also yields a clean layout each time - drops the embedded title + legend,
        // pins the 20x/100x markers, and sets a data-fit range - with no reliance
        // on the baked chart's buggy annotation state.
        var setScale = function (s) {
          gds.forEach(function (gd) {
            var mx = maxOfFullData(gd) || 1;
            Plotly.react(gd, gd.data, aggregatedLayout(s, mx), {responsive: true, displayModeBar: false});
          });
        };
        h2.insertAdjacentElement('afterend', scaleToggleEl('linear', setScale));
        setScale('linear');
      }
    } catch (e) { /* leave the baked chart untouched if the API shape surprises us */ }
  });
}

/* ---- lazy per-sample coverage plots ---- */
var covScale = 'linear';      // 'linear' (default) | 'log'
var currentSample = null;
// Log10 range for a depth axis - matches Python's _log_depth_range: floor at
// 1 (10^0), top = 1.5x the data max but never below 150, so 20x/100x always show.
function logDepthRange(maxDepth) {
  const top = Math.max(150, (maxDepth || 0) * 1.5);
  return [0, Math.log10(top)];
}
function coverageLayout(title, maxDepth, scale) {
  const dark = isDark();
  const isLog = scale === 'log';
  // Y-axis caption is just "Depth" in both modes; the Linear/Log10 toggle states
  // the scale, and log ticks are pinned to powers of ten (dtick 1) so there are
  // no bare "2"/"5" minor labels to mistake for a linear scale.
  const yaxis = isLog
    ? {title: {text: 'Depth'}, type: 'log', range: logDepthRange(maxDepth), dtick: 1,
       gridcolor: 'rgba(137,135,129,0.25)', zeroline: false}
    : {title: {text: 'Depth'}, type: 'linear', rangemode: 'tozero',
       gridcolor: 'rgba(137,135,129,0.25)', zeroline: false};
  return {
    autosize: true, height: FIG_HEIGHT, title: {text: title},
    margin: {l: 60, r: 30, t: 50, b: 60},
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: {family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 13,
           color: dark ? '#c3c2b7' : '#52514e'},
    xaxis: {title: {text: 'Genome position'}, gridcolor: 'rgba(137,135,129,0.25)', zeroline: false},
    yaxis: yaxis,
    hovermode: 'x unified', showlegend: false,
    shapes: [
      {type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 20, y1: 20, line: {color: '#898781', width: 1, dash: 'dot'}},
      {type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 100, y1: 100, line: {color: '#898781', width: 1, dash: 'dot'}}
    ],
    annotations: [
      {xref: 'paper', x: 1, y: 20, yref: 'y', text: '20x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}},
      {xref: 'paper', x: 1, y: 100, yref: 'y', text: '100x', showarrow: false, xanchor: 'right', font: {color: '#898781', size: 11}}
    ]
  };
}
function showSample(sample) {
  currentSample = sample;
  const entries = COVERAGE[sample] || [];
  const holder = document.getElementById('sample-coverage');
  holder.innerHTML = '';
  entries.forEach(function (e, i) {
    const div = document.createElement('div');
    div.id = 'sample-cov-' + i;
    holder.appendChild(div);
    const label = e.label ? (sample + ' — ' + e.label) : sample;
    const maxDepth = e.y.length ? Math.max.apply(null, e.y) : 0;
    // Honest zeros: on a linear axis plot raw depth (0 sits on the baseline);
    // on a log axis break the line at true zeros (null) so no-coverage reads as a
    // gap, never as a false depth of 1.
    const y = covScale === 'log'
      ? e.y.map(function (v) { return v > 0 ? v : null; })
      : e.y;
    const trace = {x: e.x, y: y, mode: 'lines', line: {width: 2, color: '#2a78d6'},
                   name: label, connectgaps: false};
    Plotly.newPlot(div, [trace], coverageLayout(label, maxDepth, covScale),
                   {responsive: true, displayModeBar: false});
  });
  renderSampleStats(sample);
}
// Wrap the "Coverage by sample" card body in a collapsed <details> so run-level
// cards lead and per-sample detail is opt-in. Returns the <details> (or null).
function makeSampleSectionCollapsible() {
  var card = Array.from(document.querySelectorAll('main .card')).find(function (c) {
    var h = c.querySelector('h2');
    return h && /Coverage by sample/i.test(h.textContent);
  });
  if (!card) return null;
  var h2 = card.querySelector('h2');
  if (!h2) return null;
  var details = document.createElement('details');
  details.className = 'sample-details';
  var summary = document.createElement('summary');
  var title = document.createElement('span');
  title.className = 'sample-summary-title';
  title.textContent = 'By sample';
  var hint = document.createElement('span');
  hint.className = 'sample-summary-hint';
  hint.textContent = 'Coverage and depth for one sample';
  summary.appendChild(title);
  summary.appendChild(hint);
  details.appendChild(summary);
  var body = document.createElement('div');
  body.className = 'sample-details-body';
  var node = h2.nextSibling;                 // move everything after the heading into the body
  while (node) { var next = node.nextSibling; body.appendChild(node); node = next; }
  details.appendChild(body);
  h2.remove();
  card.appendChild(details);
  return details;
}
function injectCoverageToggle() {
  var sel = document.getElementById('sampleSelect');
  if (!sel) return;
  var toggle = scaleToggleEl('linear', function (s) {
    covScale = s;
    if (currentSample) showSample(currentSample);
  });
  toggle.style.margin = '0 0 14px 8px';
  sel.insertAdjacentElement('afterend', toggle);
}
// Shared stats-table renderer (Mapped column shown as a rate; other numbers grouped).
function statsTableHTML(rows) {
  if (!rows.length) return '';
  var cols = Object.keys(rows[0]);
  var html = '<table class="stats-table"><thead><tr>';
  cols.forEach(function (c) { html += '<th>' + (c === 'Mapped reads' ? 'Mapped %' : c) + '</th>'; });
  html += '</tr></thead><tbody>';
  rows.forEach(function (r) {
    html += '<tr>';
    cols.forEach(function (c) {
      var val;
      if (c === 'Mapped reads') {
        var res = mappedPct(parseFloat(r['Total reads']), parseFloat(r['QC-passed reads']), parseFloat(r['Mapped reads']));
        val = res ? fmtPct(res.pct) : fmtNumberString(r[c]);
      } else {
        val = fmtNumberString(r[c]);
      }
      html += '<td>' + val + '</td>';
    });
    html += '</tr>';
  });
  return html + '</tbody></table>';
}
function renderSampleStats(sample) {
  var box = document.getElementById('sample-stats');
  if (box) box.innerHTML = statsTableHTML(STATS[sample] || []);
}

/* ---- per-segment view: transpose of per-sample (pick a segment, see it across
        every sample). Only built for segmented genomes. ---- */
var segCovScale = 'linear';
var currentSegment = null;
function showSegment(segment) {
  currentSegment = segment;
  var holder = document.getElementById('segment-coverage');
  if (!holder) return;
  holder.innerHTML = '';
  var i = 0;
  Object.keys(COVERAGE).forEach(function (sample) {
    var entry = (COVERAGE[sample] || []).filter(function (e) { return (e.label || '') === segment; })[0];
    if (!entry) return;
    var div = document.createElement('div');
    div.id = 'segment-cov-' + i; i++;
    holder.appendChild(div);
    var label = segment ? (sample + ' — ' + segment) : sample;   // unsegmented: sample only
    var maxDepth = entry.y.length ? Math.max.apply(null, entry.y) : 0;
    var y = segCovScale === 'log' ? entry.y.map(function (v) { return v > 0 ? v : null; }) : entry.y;
    var trace = {x: entry.x, y: y, mode: 'lines', line: {width: 2, color: '#2a78d6'}, name: label, connectgaps: false};
    Plotly.newPlot(div, [trace], coverageLayout(label, maxDepth, segCovScale), {responsive: true, displayModeBar: false});
  });
  // For a segmented run the stats table (one row per sample for this segment) is a
  // useful focused subset; for an unsegmented run it would just duplicate the
  // Assembly-statistics table, so it is omitted.
  var box = document.getElementById('segment-stats');
  if (box) box.innerHTML = segment ? statsTableHTML(rowsForSegment(segment)) : '';
}
function rowsForSegment(segment) {
  var rows = [];
  Object.keys(STATS).forEach(function (sample) {
    (STATS[sample] || []).forEach(function (r) { if ((r.Segment || '') === segment) rows.push(r); });
  });
  return rows;
}
// Build the "By segment" accordion (mirror of "By sample"), inserted right after
// the per-sample card. For segmented genomes it has a Segment selector and shows
// the chosen segment across every sample; for unsegmented genomes there is one
// (whole-genome) segment, so the selector is omitted and it simply stacks every
// sample's coverage.
function buildSegmentSection() {
  var samples = Object.keys(COVERAGE);
  if (!samples.length) return;
  var labels = (COVERAGE[samples[0]] || []).map(function (e) { return e.label; });
  if (!labels.length) return;                          // no coverage entries at all
  var labeled = labels.filter(function (l) { return l && l.length; });
  var isSegmented = labeled.length > 0;
  var segList = isSegmented ? labeled : [''];          // [''] = single whole genome
  var segChoice = segList[0];

  var sampleCard = document.querySelector('.sample-details');
  sampleCard = sampleCard ? sampleCard.closest('.card') : null;
  if (!sampleCard) return;

  var card = document.createElement('div');
  card.className = 'card';
  var details = document.createElement('details');
  details.className = 'sample-details';
  var summary = document.createElement('summary');
  var title = document.createElement('span'); title.className = 'sample-summary-title'; title.textContent = 'By segment';
  var hint = document.createElement('span'); hint.className = 'sample-summary-hint';
  hint.textContent = isSegmented ? 'Coverage and depth for one segment, across samples'
                                 : 'Coverage and depth across all samples';
  summary.appendChild(title); summary.appendChild(hint);
  details.appendChild(summary);

  var body = document.createElement('div'); body.className = 'sample-details-body';
  var controls = document.createElement('div');
  if (segList.length > 1) {                            // only offer a selector when there's a choice
    var lab = document.createElement('label'); lab.className = 'control'; lab.setAttribute('for', 'segmentSelect'); lab.textContent = 'Segment';
    var sel = document.createElement('select'); sel.id = 'segmentSelect';
    segList.forEach(function (sg) { var o = document.createElement('option'); o.value = sg; o.textContent = sg; sel.appendChild(o); });
    sel.addEventListener('change', function () { segChoice = sel.value; showSegment(segChoice); });
    controls.appendChild(lab); controls.appendChild(sel);
  }
  var toggle = scaleToggleEl('linear', function (s) { segCovScale = s; showSegment(segChoice); });
  toggle.style.margin = '0 0 14px 8px';
  controls.appendChild(toggle);
  body.appendChild(controls);
  var stats = document.createElement('div'); stats.id = 'segment-stats'; body.appendChild(stats);
  var cov = document.createElement('div'); cov.id = 'segment-coverage'; cov.className = 'plot-holder'; body.appendChild(cov);
  details.appendChild(body);
  card.appendChild(details);
  sampleCard.insertAdjacentElement('afterend', card);

  var rendered = false;
  details.addEventListener('toggle', function () {
    if (details.open && !rendered) { rendered = true; showSegment(segChoice); }
  });
}

/* ---- init ---- */
window.addEventListener('load', function () {
  themeCharts();
  formatStatsTable();
  enhanceSortHeaders();
  enhanceBakedCharts();
  // The eyebrow over the collapsible section reads "Details" (it now holds both a
  // by-sample and a by-segment view, not only per-sample).
  document.querySelectorAll('.section-title').forEach(function (s) {
    if (/per[-\\s]?sample/i.test(s.textContent)) s.textContent = 'Details';
  });
  var details = makeSampleSectionCollapsible();
  injectCoverageToggle();
  buildSegmentSection();
  const first = document.getElementById('sampleSelect');
  if (details) {
    // Lazy-render: only build the per-sample plot when the section is first opened.
    var rendered = false;
    details.addEventListener('toggle', function () {
      if (details.open && !rendered && first && first.value) { rendered = true; showSample(first.value); }
    });
  } else if (first && first.value) {
    showSample(first.value);
  }
});
'''

def patch(text, fn):
    assert text.count(CSS_ANCHOR) == 1, f"{fn}: CSS anchor count != 1"
    assert text.count(JS_OLD) == 1, f"{fn}: JS block count != 1 (is the file pristine?)"
    text = text.replace(CSS_ANCHOR, CSS_ADD, 1)
    text = text.replace(JS_OLD, JS_NEW, 1)
    return text

def main():
    write = "--write" in sys.argv
    for fn in FILES:
        p = ROOT / fn
        text = p.read_text(encoding="utf-8")
        new = patch(text, fn)
        assert CSS_ADD in new and JS_NEW in new, f"{fn}: replacement missing"
        print(f"{fn}: OK  ({len(text)} -> {len(new)} bytes)")
        if write:
            p.write_text(new, encoding="utf-8")
            print("  written")
    print("DRY RUN (pass --write to apply)" if not write else "ALL WRITTEN")

if __name__ == "__main__":
    main()
