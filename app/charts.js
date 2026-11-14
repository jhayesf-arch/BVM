// Small dependency-free canvas line charts with threshold lines and a probe cursor.

const C = {
  axis: "#3a4658", grid: "#222b37", text: "#8a97a8", line: "#4cc2ff", line2: "#8a97a8", probe: "#ffffff",
};

function niceTicks(lo, hi, count = 4) {
  const span = hi - lo || Math.abs(hi) || 1;
  const step0 = span / count, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= count) || 10 * mag;
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) ticks.push(+v.toPrecision(12));
  return ticks;
}
export function fmt(v, digits = 3) {
  if (v == null || !isFinite(v)) return "—";
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(1);
  if (a >= 100) return v.toFixed(0);
  return (+v.toPrecision(digits)).toString();
}

export class LineChart {
  /**
   * opts: { title, unit, thresholds: [{y, label, color}], onProbe(frac) }
   */
  constructor(canvas, opts) {
    this.canvas = canvas;
    this.opts = opts;
    this.ctx = canvas.getContext("2d");
    this.series = [];
    this.probe = null;
    this.marker = null;
    new ResizeObserver(() => this.draw()).observe(canvas);
    if (opts.onProbe) {
      const handler = (e) => {
        if (e.type === "pointermove" && e.buttons !== 1 && !opts.hoverProbe) return;
        const r = canvas.getBoundingClientRect(), pad = this._pad();
        const frac = (e.clientX - r.left - pad.l) / (r.width - pad.l - pad.r);
        if (frac >= 0 && frac <= 1) opts.onProbe(frac);
      };
      canvas.addEventListener("pointerdown", handler);
      canvas.addEventListener("pointermove", handler);
    }
  }

  _pad() { return { l: 46, r: 10, t: 22, b: 26 }; }

  set(series, { xLabel, probe = null, marker = null, xRange } = {}) {
    this.series = series; // [{x: [], y: [], color?, dashed?, label?}]
    this.xLabel = xLabel;
    this.probe = probe;
    this.marker = marker;
    this.xRange = xRange;
    this.draw();
  }

  setProbe(frac) { this.probe = frac; this.draw(); }

  draw() {
    const cv = this.canvas, dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    if (!w || !h) return;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    }
    const ctx = this.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const pad = this._pad(), pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
    ctx.font = "11px Inter, system-ui, sans-serif";

    ctx.fillStyle = "#cfd8e3";
    ctx.textAlign = "left";
    ctx.fillText(`${this.opts.title}${this.opts.unit ? ` (${this.opts.unit})` : ""}`, pad.l, 14);
    if (!this.series.length) return;

    const xs = this.series.flatMap((s) => Array.from(s.x)), ys = this.series.flatMap((s) => Array.from(s.y));
    let [x0, x1] = this.xRange || [Math.min(...xs), Math.max(...xs)];
    let y0 = Math.min(...ys), y1 = Math.max(...ys);
    for (const t of this.opts.thresholds || []) if (t.alwaysShow) { y0 = Math.min(y0, t.y); y1 = Math.max(y1, t.y); }
    if (this.opts.zeroBased) y0 = Math.min(0, y0);
    if (y1 - y0 < 1e-9 * Math.max(1, Math.abs(y1))) { const m = Math.abs(y1) * 0.1 || 1; y0 -= m; y1 += m; }
    const yPad = (y1 - y0) * 0.08;
    y1 += yPad; if (!this.opts.zeroBased || y0 < 0) y0 -= yPad;
    const X = (x) => pad.l + ((x - x0) / (x1 - x0 || 1)) * pw;
    const Y = (y) => pad.t + ph - ((y - y0) / (y1 - y0)) * ph;
    this.map = { X, Y, x0, x1 };

    // grid + ticks
    ctx.strokeStyle = C.grid; ctx.lineWidth = 1; ctx.fillStyle = C.text; ctx.textAlign = "right";
    for (const t of niceTicks(y0, y1)) {
      ctx.beginPath(); ctx.moveTo(pad.l, Y(t)); ctx.lineTo(pad.l + pw, Y(t)); ctx.stroke();
      ctx.fillText(fmt(t), pad.l - 6, Y(t) + 4);
    }
    ctx.textAlign = "center";
    for (const t of niceTicks(x0, x1, 5)) ctx.fillText(fmt(t), X(t), h - pad.b + 14);
    if (this.xLabel) { ctx.textAlign = "right"; ctx.fillText(this.xLabel, pad.l + pw, h - 2); }

    // threshold lines
    for (const t of this.opts.thresholds || []) {
      if (t.y < y0 || t.y > y1) continue;
      ctx.save(); ctx.setLineDash([4, 4]); ctx.strokeStyle = t.color; ctx.beginPath();
      ctx.moveTo(pad.l, Y(t.y)); ctx.lineTo(pad.l + pw, Y(t.y)); ctx.stroke(); ctx.restore();
      ctx.fillStyle = t.color; ctx.textAlign = "left"; ctx.fillText(t.label, pad.l + 4, Y(t.y) - 4);
    }
    // bands (e.g. ischemia ranges), drawn under lines
    for (const b of this.opts.bands || []) {
      const ya = Math.max(b.from, y0), yb = Math.min(b.to, y1);
      if (yb <= ya) continue;
      ctx.fillStyle = b.color; ctx.fillRect(pad.l, Y(yb), pw, Y(ya) - Y(yb));
    }

    // series
    for (const s of this.series) {
      ctx.save();
      ctx.strokeStyle = s.color || C.line; ctx.lineWidth = s.width || 1.8;
      if (s.dashed) ctx.setLineDash([5, 4]);
      ctx.beginPath();
      for (let i = 0; i < s.x.length; i++) (i ? ctx.lineTo : ctx.moveTo).call(ctx, X(s.x[i]), Y(s.y[i]));
      ctx.stroke();
      if (s.dots) { ctx.fillStyle = s.color || C.line; for (let i = 0; i < s.x.length; i++) { ctx.beginPath(); ctx.arc(X(s.x[i]), Y(s.y[i]), 2.2, 0, 7); ctx.fill(); } }
      ctx.restore();
    }
    // legend
    const labelled = this.series.filter((s) => s.label);
    let ly = pad.t + ph - 6;
    ctx.textAlign = "right";
    for (const s of labelled) {
      ctx.fillStyle = s.color || C.line; ctx.fillText(s.label, pad.l + pw - 4, ly);
      ly -= 13;
    }

    // probe cursor
    if (this.probe != null) {
      const px = pad.l + this.probe * pw;
      ctx.strokeStyle = C.probe; ctx.globalAlpha = 0.6; ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      ctx.globalAlpha = 1;
    }
    // current-configuration marker
    if (this.marker) {
      ctx.fillStyle = "#fff"; ctx.strokeStyle = C.line; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(X(this.marker.x), Y(this.marker.y), 5, 0, 7); ctx.fill(); ctx.stroke();
    }
  }
}
