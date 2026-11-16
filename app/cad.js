// CAD import: read a vessel surface (STL, OBJ, STEP/IGES) and reduce it to the 1D model's inputs,
// a centerline and an equivalent-circle lumen radius at each station.
//
// Works for a lumen (fluid-domain) solid, an open lumen surface, or a hollow vessel-wall solid (inner and
// outer surfaces: the inner one is the lumen, and the gap gives the wall thickness). The method:
//   1. weld vertices; find the two ends (open boundary rims, or the flat end caps of a closed solid);
//   2. march from the inlet end: slice the mesh with a plane normal to the current direction, keep the
//      cross-section loop that contains the marching point (the smallest one, so the lumen and not the outer
//      wall), step to its centroid, and turn the direction toward the newest centroids;
//   3. smooth the centroid path, resample it, and slice again normal to the smoothed tangent to get the
//      area A of each section; the model uses the equivalent radius √(A/π).
// Only the main path is followed: side branches in the CAD are cut across and do not take flow.
// Pure JS (no three.js) so the self-test can run it; STEP/IGES are tessellated by occt-import-js (OpenCascade,
// LGPL), loaded from a CDN only when such a file is opened.

const OCCT_URL = "https://cdn.jsdelivr.net/npm/occt-import-js@0.0.23/dist/";
export const UNIT_SCALE = { mm: 1e-3, cm: 1e-2, m: 1, in: 0.0254, um: 1e-6 };

// ───────────────────────── file readers (all return a triangle soup, 9 floats per triangle) ─────────────────────────
export function parseSTL(buf) {
  const dv = new DataView(buf);
  if (buf.byteLength >= 84) {
    const nTri = dv.getUint32(80, true);
    if (84 + nTri * 50 === buf.byteLength) {
      const out = new Float32Array(nTri * 9);
      for (let i = 0; i < nTri; i++) for (let k = 0; k < 9; k++) out[9 * i + k] = dv.getFloat32(84 + 50 * i + 12 + 4 * k, true);
      return out;
    }
  }
  const text = new TextDecoder().decode(buf);
  const nums = [];
  const re = /vertex\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)/g;
  let m;
  while ((m = re.exec(text))) nums.push(+m[1], +m[2], +m[3]);
  if (!nums.length || nums.length % 9) throw new Error("Not a readable STL file (no complete facets found).");
  return Float32Array.from(nums);
}

export function parseOBJ(text) {
  const v = [], out = [];
  for (const line of text.split(/\r?\n/)) {
    const t = line.trim().split(/\s+/);
    if (t[0] === "v") v.push([+t[1], +t[2], +t[3]]);
    else if (t[0] === "f") {
      const ids = t.slice(1).map((w) => { const i = parseInt(w, 10); return i < 0 ? v.length + i : i - 1; });
      for (let k = 1; k + 1 < ids.length; k++) for (const i of [ids[0], ids[k], ids[k + 1]]) out.push(...v[i]);
    }
  }
  if (!out.length) throw new Error("No faces found in the OBJ file.");
  return Float32Array.from(out);
}

let occtPromise = null;
function loadOcct() {
  occtPromise ??= new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = `${OCCT_URL}occt-import-js.js`;
    s.onload = () => window.occtimportjs({ locateFile: (f) => OCCT_URL + f }).then(resolve, reject);
    s.onerror = () => reject(new Error("Could not load the STEP/IGES reader (needs an internet connection). Export the design as STL instead."));
    document.head.appendChild(s);
  });
  return occtPromise;
}

async function parseBrep(buf, ext) {
  const occt = await loadOcct();
  const bytes = new Uint8Array(buf);
  const res = ext === "iges" || ext === "igs" ? occt.ReadIgesFile(bytes, null) : occt.ReadStepFile(bytes, null);
  if (!res.success || !res.meshes.length) throw new Error("The STEP/IGES file could not be tessellated.");
  const out = [];
  for (const m of res.meshes) {
    const p = m.attributes.position.array, idx = m.index.array;
    for (const i of idx) out.push(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
  }
  return Float32Array.from(out);
}

/** Read a File (browser) into a triangle soup. */
export async function readCadFile(file) {
  const ext = file.name.split(".").pop().toLowerCase();
  const buf = await file.arrayBuffer();
  if (ext === "stl") return parseSTL(buf);
  if (ext === "obj") return parseOBJ(new TextDecoder().decode(buf));
  if (["step", "stp", "iges", "igs"].includes(ext)) return parseBrep(buf, ext);
  throw new Error(`Unsupported file type .${ext}. Use STL, OBJ, STEP or IGES.`);
}

// ───────────────────────── mesh utilities ─────────────────────────
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => { const m = norm(a) || 1; return [a[0] / m, a[1] / m, a[2] / m]; };

/** Merge coincident vertices. Returns { V: Float64Array xyz, F: Uint32Array (3 per triangle), diag }. */
export function weld(tri) {
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < tri.length; i += 3) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], tri[i + k]); hi[k] = Math.max(hi[k], tri[i + k]); }
  const diag = norm(sub(hi, lo)) || 1, q = diag * 1e-7;
  const map = new Map(), V = [], F = [];
  const nt = tri.length / 9;
  for (let t = 0; t < nt; t++) {
    const ids = [];
    for (let c = 0; c < 3; c++) {
      const o = 9 * t + 3 * c, x = tri[o], y = tri[o + 1], z = tri[o + 2];
      const key = `${Math.round(x / q)},${Math.round(y / q)},${Math.round(z / q)}`;
      let id = map.get(key);
      if (id === undefined) { id = V.length / 3; map.set(key, id); V.push(x, y, z); }
      ids.push(id);
    }
    if (ids[0] !== ids[1] && ids[1] !== ids[2] && ids[0] !== ids[2]) F.push(...ids);
  }
  return { V: Float64Array.from(V), F: Uint32Array.from(F), diag };
}

const vtx = (M, i) => [M.V[3 * i], M.V[3 * i + 1], M.V[3 * i + 2]];

/** Planar loops of a polyline graph given as edge-key adjacency. Returns arrays of keys. */
function chainLoops(adj) {
  const seen = new Set(), loops = [];
  for (const start of adj.keys()) {
    if (seen.has(start)) continue;
    const loop = [start];
    seen.add(start);
    let prev = null, cur = start, closed = false;
    for (;;) {
      const nb = (adj.get(cur) ?? []).find((k) => k !== prev && (!seen.has(k) || (k === start && loop.length > 2)));
      if (nb === undefined) break;
      if (nb === start) { closed = true; break; }
      loop.push(nb); seen.add(nb); prev = cur; cur = nb;
    }
    if (closed && loop.length >= 3) loops.push(loop);
  }
  return loops;
}

function planeBasis(n) {
  const a = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  const u = unit(cross(n, a));
  return [u, cross(n, u)];
}

/** Polygon properties of a closed 3D loop lying in the plane (c, n): area, centroid, perimeter, 2D points. */
function loopProps(pts, c, n) {
  const [u, w] = planeBasis(n);
  const p2 = pts.map((p) => { const d = sub(p, c); return [dot(d, u), dot(d, w)]; });
  let a = 0, cx = 0, cy = 0, per = 0, ixx = 0, iyy = 0;
  for (let i = 0; i < p2.length; i++) {
    const [x0, y0] = p2[i], [x1, y1] = p2[(i + 1) % p2.length], cr = x0 * y1 - x1 * y0;
    a += cr; cx += (x0 + x1) * cr; cy += (y0 + y1) * cr; per += Math.hypot(x1 - x0, y1 - y0);
    ixx += cr * (y0 * y0 + y0 * y1 + y1 * y1); iyy += cr * (x0 * x0 + x0 * x1 + x1 * x1);
  }
  a /= 2;
  if (Math.abs(a) < 1e-30) return null;
  cx /= 6 * a; cy /= 6 * a;
  // polar second moment about the centroid (sign of a cancels)
  const ip = (ixx + iyy) / 12 / a * Math.abs(a) - Math.abs(a) * (cx * cx + cy * cy);
  return { area: Math.abs(a), ip, centroid: add(c, add(mul(u, cx), mul(w, cy))), perimeter: per, p2, u, w, c };
}

function contains(L, q) {
  const d = sub(q, L.c), x = dot(d, L.u), y = dot(d, L.w), P = L.p2;
  let inside = false;
  for (let i = 0, j = P.length - 1; i < P.length; j = i++) {
    const [xi, yi] = P[i], [xj, yj] = P[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

// per-triangle centroid and bounding radius, so a slice only touches triangles near its plane and point
function triBounds(M) {
  if (M.TB) return M.TB;
  const F = M.F, nf = F.length / 3, TB = new Float64Array(4 * nf), V = M.V;
  for (let t = 0; t < nf; t++) {
    const a = 3 * F[3 * t], b = 3 * F[3 * t + 1], c = 3 * F[3 * t + 2];
    const cx = (V[a] + V[b] + V[c]) / 3, cy = (V[a + 1] + V[b + 1] + V[c + 1]) / 3, cz = (V[a + 2] + V[b + 2] + V[c + 2]) / 3;
    let r = 0;
    for (const o of [a, b, c]) r = Math.max(r, Math.hypot(V[o] - cx, V[o + 1] - cy, V[o + 2] - cz));
    TB.set([cx, cy, cz, r], 4 * t);
  }
  return (M.TB = TB);
}

/** Closed cross-section loops of the mesh with the plane through c, normal n, within distance `cut` of c. */
export function slice(M, c, n, cut = Infinity) {
  const nv = M.V.length / 3, V = M.V, F = M.F, TB = triBounds(M), nf = F.length / 3;
  const adj = new Map(), pt = new Map();
  const dist = (i) => (V[3 * i] - c[0]) * n[0] + (V[3 * i + 1] - c[1]) * n[1] + (V[3 * i + 2] - c[2]) * n[2];
  const link = (k0, k1) => { const a = adj.get(k0); if (a) a.push(k1); else adj.set(k0, [k1]); };
  const v = [0, 0, 0], d = [0, 0, 0], ks = [];
  for (let t = 0; t < nf; t++) {
    const ox = TB[4 * t] - c[0], oy = TB[4 * t + 1] - c[1], oz = TB[4 * t + 2] - c[2], tr = TB[4 * t + 3];
    if (Math.abs(ox * n[0] + oy * n[1] + oz * n[2]) > tr) continue;
    if (cut < Infinity && ox * ox + oy * oy + oz * oz > (cut + tr) ** 2) continue;
    for (let e = 0; e < 3; e++) { v[e] = F[3 * t + e]; d[e] = dist(v[e]); }
    ks.length = 0;
    for (let e = 0; e < 3; e++) {
      const f = (e + 1) % 3;
      if ((d[e] > 0) === (d[f] > 0)) continue;
      const a = v[e], b = v[f], k = a < b ? a * nv + b : b * nv + a;
      if (!pt.has(k)) {
        const s = d[e] / (d[e] - d[f]);
        pt.set(k, [V[3 * a] + (V[3 * b] - V[3 * a]) * s, V[3 * a + 1] + (V[3 * b + 1] - V[3 * a + 1]) * s, V[3 * a + 2] + (V[3 * b + 2] - V[3 * a + 2]) * s]);
      }
      ks.push(k);
    }
    if (ks.length !== 2) continue;
    link(ks[0], ks[1]); link(ks[1], ks[0]);
  }
  return chainLoops(adj).map((keys) => loopProps(keys.map((k) => pt.get(k)), c, n)).filter(Boolean);
}

// slice near the point first (fast, and it ignores other parts of a looping vessel); widen if nothing closes
function sliceNear(M, c, n, r) {
  const loops = slice(M, c, n, 4 * r);
  return loops.length ? loops : slice(M, c, n);
}

/** Open boundary rims (edges used by one triangle), as loops of vertex ids. */
function boundaryLoops(M) {
  const count = new Map(), nv = M.V.length / 3, F = M.F;
  for (let t = 0; t < F.length; t += 3) for (let e = 0; e < 3; e++) {
    const a = F[t + e], b = F[t + ((e + 1) % 3)], k = a < b ? a * nv + b : b * nv + a;
    count.set(k, (count.get(k) ?? 0) + 1);
  }
  const adj = new Map();
  for (const [k, c] of count) if (c === 1) {
    const a = Math.floor(k / nv), b = k % nv;
    adj.set(a, [...(adj.get(a) ?? []), b]);
    adj.set(b, [...(adj.get(b) ?? []), a]);
  }
  return chainLoops(adj);
}

function newell(pts) {
  let n = [0, 0, 0];
  for (let i = 0; i < pts.length; i++) n = add(n, cross(pts[i], pts[(i + 1) % pts.length]));
  return n;
}

function pca(M) {
  const nv = M.V.length / 3;
  let m = [0, 0, 0];
  for (let i = 0; i < nv; i++) m = add(m, vtx(M, i));
  m = mul(m, 1 / nv);
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < nv; i++) { const d = sub(vtx(M, i), m); for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += d[a] * d[b]; }
  let v = [1, 0.7, 0.3];
  for (let k = 0; k < 100; k++) v = unit([dot(C[0], v), dot(C[1], v), dot(C[2], v)]);
  return { mean: m, axis: v };
}

/**
 * The two ends as { point, inward normal, area }. Open surfaces: the boundary rims (concentric inner/outer
 * rims of a hollow wall are grouped, keeping the inner one). Closed solids: the largest flat face cluster
 * near each extreme of the principal axis (the end cap).
 */
function findEnds(M) {
  const { mean, axis } = pca(M);
  const rims = boundaryLoops(M).map((ids) => {
    const P = ids.map((i) => vtx(M, i)), c = mul(P.reduce(add, [0, 0, 0]), 1 / P.length), nw = newell(P.map((p) => sub(p, c)));
    return { point: c, normal: unit(nw), area: norm(nw) / 2 };
  }).filter((r) => r.area > 0);
  if (rims.length >= 2) {
    const groups = [];
    for (const r of rims.sort((a, b) => a.area - b.area)) {
      const g = groups.find((x) => norm(sub(x.point, r.point)) < Math.sqrt(r.area / Math.PI));
      if (!g) groups.push({ ...r, total: r.area }); else g.total += r.area;   // inner rim already holds the group's place
    }
    groups.sort((a, b) => b.total - a.total);
    if (groups.length >= 2) {
      const ends = groups.slice(0, 2);
      const mid = mul(add(ends[0].point, ends[1].point), 0.5);
      for (const e of ends) if (dot(e.normal, sub(mid, e.point)) < 0) e.normal = mul(e.normal, -1);
      return { ends, open: true };
    }
  }
  // closed: the end caps are the two largest flat faces (same normal and plane offset), wherever they are,
  // so a U-bend or siphon works too
  const F = M.F, clusters = new Map(), qa = 1e-3 * M.diag;
  for (let t = 0; t < F.length; t += 3) {
    const A = vtx(M, F[t]), B = vtx(M, F[t + 1]), C = vtx(M, F[t + 2]);
    const nn = cross(sub(B, A), sub(C, A)), ar = norm(nn) / 2;
    if (!(ar > 0)) continue;
    const nu = mul(nn, 1 / (2 * ar)), cen = mul(add(add(A, B), C), 1 / 3);
    const key = `${Math.round(nu[0] * 100)},${Math.round(nu[1] * 100)},${Math.round(nu[2] * 100)},${Math.round(dot(cen, nu) / qa)}`;
    let cl = clusters.get(key);
    if (!cl) clusters.set(key, (cl = { n: nu, tris: [] }));
    cl.tris.push(t / 3, cen, ar);
  }
  // split each coplanar group into connected pieces (the two caps of a U-bend share one plane)
  const pieces = [];
  for (const cl of clusters.values()) {
    if (cl.tris.length < 6) continue;
    const parent = new Map(), find = (x) => { while (parent.get(x) !== x) { parent.set(x, parent.get(parent.get(x))); x = parent.get(x); } return x; };
    for (let i = 0; i < cl.tris.length; i += 3) {
      const f = 3 * cl.tris[i], ids = [F[f], F[f + 1], F[f + 2]];
      for (const v of ids) if (!parent.has(v)) parent.set(v, v);
      parent.set(find(ids[1]), find(ids[0])); parent.set(find(ids[2]), find(ids[0]));
    }
    const byRoot = new Map();
    for (let i = 0; i < cl.tris.length; i += 3) {
      const root = find(F[3 * cl.tris[i]]);
      let p = byRoot.get(root);
      if (!p) byRoot.set(root, (p = { n: cl.n, area: 0, sum: [0, 0, 0], tris: 0, verts: new Set() }));
      p.area += cl.tris[i + 2]; p.sum = add(p.sum, mul(cl.tris[i + 1], cl.tris[i + 2])); p.tris++;
      for (let e = 0; e < 3; e++) p.verts.add(F[3 * cl.tris[i] + e]);
    }
    pieces.push(...byRoot.values());
  }
  // a cap is compact (disk, annulus, square); a flat strip along a straight or conical wall is long and thin
  const round = (k) => {
    const [u, w] = planeBasis(k.n), P = [...k.verts].map((i) => vtx(M, i));
    const c = mul(P.reduce(add, [0, 0, 0]), 1 / P.length);
    let a = 0, b = 0, d = 0;
    for (const q of P) {
      const p = sub(q, c), x = dot(p, u), y = dot(p, w);
      a += x * x; b += x * y; d += y * y;
    }
    const tr = a + d, det = a * d - b * b, disc = Math.sqrt(Math.max(tr * tr / 4 - det, 0));
    return tr / 2 - disc > 0.3 * (tr / 2 + disc);
  };
  const flats = pieces.filter((k) => k.tris >= 2 && k.verts.size >= 4 && round(k)).sort((a, b) => b.area - a.area)
    .map((k) => ({ point: mul(k.sum, 1 / k.area), normal: mul(k.n, -1), area: k.area, cap: true }));
  if (flats.length >= 2) {
    const a = flats[0], b = flats.slice(1).find((f) => norm(sub(f.point, a.point)) > 2 * Math.sqrt(a.area / Math.PI));
    if (b && b.area > 0.2 * a.area) {
      // the cap normals came from the winding: flip them all if the solid is wound inside-out
      let vol = 0;
      for (let t = 0; t < F.length; t += 3) vol += dot(vtx(M, F[t]), cross(vtx(M, F[t + 1]), vtx(M, F[t + 2])));
      if (vol < 0) for (const e of [a, b]) e.normal = mul(e.normal, -1);
      return { ends: [a, b], open: false };
    }
  }
  // rounded or irregular ends: start just inside the extremes along the principal axis
  let tmin = Infinity, tmax = -Infinity;
  for (let i = 0; i < M.V.length / 3; i++) { const t = dot(sub(vtx(M, i), mean), axis); tmin = Math.min(tmin, t); tmax = Math.max(tmax, t); }
  const ends = [-1, 1].map((side) => ({ point: add(mean, mul(axis, side < 0 ? tmin : tmax)), normal: mul(axis, -side), area: 0, cap: false }));
  return { ends, open: false };
}

function pickLoop(loops, q, rPrev) {
  const inside = loops.filter((L) => contains(L, q)).sort((a, b) => a.area - b.area);
  if (inside.length) return { lumen: inside[0], outer: inside.find((L) => L.area > inside[0].area * 1.02 && contains(L, inside[0].centroid)) ?? null };
  let best = null, bd = Infinity;
  for (const L of loops) { const dd = norm(sub(L.centroid, q)); if (dd < bd) { bd = dd; best = L; } }
  return best && bd < 1.5 * rPrev ? { lumen: best, outer: null } : null;
}

function cumLength(P) {
  const s = [0];
  for (let i = 1; i < P.length; i++) s.push(s[i - 1] + norm(sub(P[i], P[i - 1])));
  return s;
}

function resample(P, s, m) {
  const L = s[s.length - 1], out = [];
  let j = 0;
  for (let i = 0; i < m; i++) {
    const si = (L * i) / (m - 1);
    while (j < s.length - 2 && s[j + 1] < si) j++;
    const t = Math.min(Math.max((si - s[j]) / (s[j + 1] - s[j] || 1), 0), 1);
    out.push(add(P[j], mul(sub(P[j + 1], P[j]), t)));
  }
  return out;
}

function gaussSmooth(P, sigmaPts) {
  if (sigmaPts < 0.5) return P;
  const n = P.length, w = Math.ceil(3 * sigmaPts), out = [P[0]];
  for (let i = 1; i < n - 1; i++) {
    let acc = [0, 0, 0], ws = 0;
    for (let k = -w; k <= w; k++) {
      // mirror about the end point so the ends stay fixed and straight
      let j = i + k, p;
      if (j < 0) p = sub(mul(P[0], 2), P[-j]);
      else if (j >= n) p = sub(mul(P[n - 1], 2), P[2 * (n - 1) - j]);
      else p = P[j];
      const wk = Math.exp(-(k * k) / (2 * sigmaPts * sigmaPts));
      acc = add(acc, mul(p, wk)); ws += wk;
    }
    out.push(mul(acc, 1 / ws));
  }
  out.push(P[n - 1]);
  return out;
}

/**
 * Radius of the circle with the same Poiseuille resistance as this section. Saint-Venant's approximation
 * R' ∝ I_p / A⁴ is exact for circles and ellipses and within ~8% for a square; the equal-area circle would
 * miss that a flattened or crescent-shaped lumen resists flow more than a round one of the same area.
 */
const hydraulicRadius = (L) => (L.ip > 0 ? (L.area ** 4 / (2 * Math.PI ** 3 * L.ip)) ** 0.25 : Math.sqrt(L.area / Math.PI));

const median3 = (a) => a.map((v, i) => (i === 0 || i === a.length - 1 ? v : [a[i - 1], v, a[i + 1]].sort((x, y) => x - y)[1]));

/** Guess units from the end-section size: arteries and mimics are 0.1–30 mm across. */
export function guessUnits(M) {
  const { ends } = findEnds(M);
  const d = 2 * Math.sqrt(Math.max(...ends.map((e) => e.area), 0) / Math.PI) || M.diag / 10;
  return d >= 0.5 ? "mm" : d >= 0.05 ? "cm" : "m";
}

/**
 * Centerline and lumen of a welded mesh (coordinates in the file's units).
 * opts: { scale (file unit → m), reverse (swap inlet/outlet), m (output stations) }.
 * Returns arrays in metres in a canonical frame (inlet at the origin, inlet→outlet chord along +z, main
 * bend in the x–z plane), plus the 4×3 transform for drawing the original mesh in that frame.
 */
export function extractLumen(M, { scale = 1e-3, reverse = false, m = 1001 } = {}) {
  const warnings = [];
  const { ends, open } = findEnds(M);
  // inlet = the larger end (arteries taper downstream) unless reversed
  ends.sort((a, b) => b.area - a.area);
  if (reverse) ends.reverse();
  const [inlet, outlet] = ends;
  const chord0 = norm(sub(outlet.point, inlet.point));
  if (!(chord0 > 0)) throw new Error("Could not find two separate vessel ends in this file.");

  // pass 1: march along the lumen
  let r0 = inlet.area > 0 ? Math.sqrt(inlet.area / Math.PI) : chord0 / 20;
  let t = inlet.normal, c = add(inlet.point, mul(t, Math.max(0.05 * r0, chord0 * 1e-4)));
  const path = [inlet.point], rad = [];
  let rPrev = r0, firstOk = false;
  for (let k = 0; k < 6000; k++) {
    const loops = sliceNear(M, c, t, rPrev);
    const pick = loops.length ? pickLoop(loops, c, rPrev) : null;
    if (!pick) {
      if (!firstOk) throw new Error("Could not find the lumen next to the inlet end. Is this a closed vessel surface?");
      break;
    }
    firstOk = true;
    const cc = pick.lumen.centroid, rr = Math.sqrt(pick.lumen.area / Math.PI);
    path.push(cc); rad.push(rr);
    rPrev = rr;
    const toEnd = sub(outlet.point, cc);
    if (dot(toEnd, outlet.normal) > -0.6 * rr && norm(toEnd) < 4 * Math.max(rr, Math.sqrt(outlet.area / Math.PI) || rr)) break;   // reached the outlet face
    // new direction: toward recent centroids, blended with the old one for stability
    const back = path[Math.max(0, path.length - 4)];
    const dn = sub(cc, back);
    if (norm(dn) > 1e-12) t = unit(add(mul(t, 0.5), mul(unit(dn), 0.5)));
    const h = Math.min(Math.max(0.25 * rr, chord0 / 3000), chord0 / 40);
    c = add(cc, mul(t, h));
  }
  path.push(outlet.point);
  if (path.length < 6) throw new Error("The lumen could not be followed from end to end.");
  if (cumLength(path).at(-1) > 3 * M.diag) throw new Error("The path through this shape is longer than the part itself, so it is probably not a single vessel lumen. Check the file contains one vessel (lumen or wall) with two open or flat ends.");

  // smooth the centroid path (σ ≈ the mean radius) and resample uniformly
  const s1 = cumLength(path), L1 = s1[s1.length - 1];
  const rMean = rad.reduce((a, b) => a + b, 0) / rad.length;
  let P = resample(path, s1, m);
  P = gaussSmooth(P, 0.5 * (rMean / L1) * (m - 1));
  // pass 2: slice normal to the smoothed tangent at ~400 stations (radius interpolated between them);
  // pass 3 rebuilds the centerline from those perpendicular cuts' centroids, which are unbiased where the
  // lumen narrows quickly (oblique cuts of a narrowing shift the centroid and fake curvature), and cuts again
  const every = Math.max(1, Math.ceil((m - 1) / 400)), idx = [];
  for (let i = 0; i < m - 1; i += every) idx.push(i);
  idx.push(m - 1);
  const ns = idx.length;
  const measure = (P) => {
    const s = cumLength(P), L = s[s.length - 1];
    const rr = new Array(ns), ra = new Array(ns), wall = new Array(ns).fill(null), circ = new Array(ns).fill(1), cen = new Array(ns).fill(null);
    let misses = 0;
    for (let k = 0; k < ns; k++) {
      const i = idx[k], a = P[Math.max(i - 1, 0)], b = P[Math.min(i + 1, m - 1)];
      // near an end, turn the cut toward the end face so it cannot clip the cap
      const wIn = Math.max(0, 1 - s[i] / (1.5 * rMean)), wOut = Math.max(0, 1 - (L - s[i]) / (1.5 * rMean));
      const tan = unit(add(add(mul(unit(sub(b, a)), Math.max(0, 1 - wIn - wOut)), mul(inlet.normal, wIn)), mul(outlet.normal, -wOut)));
      // stay just inside the end faces
      const q = i === 0 ? add(P[0], mul(tan, Math.min(0.02 * L, 0.1 * rMean))) : i === m - 1 ? sub(P[m - 1], mul(tan, Math.min(0.02 * L, 0.1 * rMean))) : P[i];
      const rg = k ? rr[k - 1] ?? rMean : rMean;
      const pick = pickLoop(sliceNear(M, q, tan, rg), q, rg);
      if (pick) {
        rr[k] = hydraulicRadius(pick.lumen);
        ra[k] = Math.sqrt(pick.lumen.area / Math.PI);
        circ[k] = (4 * Math.PI * pick.lumen.area) / pick.lumen.perimeter ** 2;
        cen[k] = pick.lumen.centroid;
        if (pick.outer) wall[k] = Math.sqrt(pick.outer.area / Math.PI) - Math.sqrt(pick.lumen.area / Math.PI);
      } else { rr[k] = null; misses++; }
    }
    return { s, L, rr, ra, wall, circ, cen, misses };
  };
  const sigma = 0.5 * (rMean / L1) * (m - 1);
  let pass = measure(P);
  const refined = [P[0], ...pass.cen.slice(1, -1).filter(Boolean), P[m - 1]];
  if (refined.length > 0.9 * ns) {
    P = gaussSmooth(resample(refined, cumLength(refined), m), sigma);
    pass = measure(P);
  }
  const { s, L, rr, ra, wall, circ } = pass;
  let misses = pass.misses;
  for (let k = 0; k < ns; k++) if (rr[k] == null) {     // fill the rare miss from neighbours
    let a = k - 1, b = k + 1;
    while (a >= 0 && rr[a] == null) a--;
    while (b < ns && rr[b] == null) b++;
    if (a < 0 && b >= ns) throw new Error("No cross-section of the lumen could be measured.");
    rr[k] = a < 0 ? rr[b] : b >= ns ? rr[a] : rr[a] + ((rr[b] - rr[a]) * (k - a)) / (b - a);
  }
  if (misses > 0.05 * ns) warnings.push(`${misses} of ${ns} cross-sections could not be cut cleanly; their radii were interpolated.`);
  const rk = median3(rr), rs = new Array(m);
  for (let k = 0; k < ns - 1; k++) for (let i = idx[k]; i <= idx[k + 1]; i++) rs[i] = rk[k] + ((rk[k + 1] - rk[k]) * (i - idx[k])) / (idx[k + 1] - idx[k] || 1);

  // canonical frame: inlet at origin, chord along +z, main bend along +x
  const ez = unit(sub(P[m - 1], P[0]));
  let bend = [0, 0, 0];
  for (const p of P) { const d = sub(p, P[0]), perp = sub(d, mul(ez, dot(d, ez))); if (norm(perp) > norm(bend)) bend = perp; }
  const ex = norm(bend) > 1e-9 * L ? unit(bend) : planeBasis(ez)[0], ey = cross(ez, ex);
  const o = P[0];
  const toCanon = (p) => { const d = sub(p, o); return [dot(d, ex) * scale, dot(d, ey) * scale, dot(d, ez) * scale]; };
  const C = P.map(toCanon);

  const wallVals = wall.filter((w) => w != null && w > 0);
  const wallMm = wallVals.length > 0.5 * ns ? (wallVals.sort((a, b) => a - b)[Math.floor(wallVals.length / 2)] * scale * 1e3) : null;
  const minCirc = Math.min(...circ);
  if (minCirc < 0.75) warnings.push(`Some cross-sections are far from circular (roundness ${minCirc.toFixed(2)}). They are replaced by the circle of equal Poiseuille resistance (Saint-Venant), which is approximate for strongly non-convex shapes; secondary flows there need 3D CFD.`);
  const areaVals = ra.filter((v) => v != null);
  const chord = norm(sub(P[m - 1], P[0]));
  return {
    s: s.map((v) => v * scale), x: C.map((p) => p[0]), y: C.map((p) => p[1]), z: C.map((p) => p[2]),
    r: rs.map((v) => v * scale), wallMm, circularity: circ, open,
    minAreaDiameterMm: 2e3 * Math.min(...areaVals) * scale,
    caps: ends.every((e) => e.cap || open),
    tortuosity: L / chord, warnings,
    frame: { origin: o, ex, ey, ez, scale },
  };
}

/** Map a triangle soup (file units) into the canonical frame of extractLumen (metres). */
export function toCanonical(tri, frame) {
  const { origin: o, ex, ey, ez, scale } = frame, out = new Float32Array(tri.length);
  for (let i = 0; i < tri.length; i += 3) {
    const d = [tri[i] - o[0], tri[i + 1] - o[1], tri[i + 2] - o[2]];
    out[i] = dot(d, ex) * scale; out[i + 1] = dot(d, ey) * scale; out[i + 2] = dot(d, ez) * scale;
  }
  return out;
}

// ───────────────────────── test shapes (self-test and the example download) ─────────────────────────
/**
 * Triangle soup of a capped tube along a centerline: path(u) → [x,y,z] for u in [0,1], radius(u) → r.
 * wall > 0 makes a hollow vessel wall (inner and outer surfaces, annular end faces).
 */
export function tubeMesh(path, radius, { nu = 160, nt = 40, wall = 0 } = {}) {
  // rotation-minimizing ring frames, so rings never twist against each other
  const frames = [];
  for (let i = 0; i <= nu; i++) {
    const u = i / nu, e = 1e-4, T = unit(sub(path(Math.min(u + e, 1)), path(Math.max(u - e, 0))));
    let N = i ? frames[i - 1][0] : planeBasis(T)[0];
    N = unit(sub(N, mul(T, dot(N, T))));
    frames.push([N, cross(T, N)]);
  }
  const ring = (u, rr) => {
    const p = path(u), [N, B] = frames[Math.round(u * nu)];
    return Array.from({ length: nt }, (_, j) => { const th = (2 * Math.PI * j) / nt; return add(p, add(mul(N, rr * Math.cos(th)), mul(B, rr * Math.sin(th)))); });
  };
  const out = [], tri = (a, b, c) => out.push(...a, ...b, ...c);
  const surface = (rf, flip) => {
    const R = Array.from({ length: nu + 1 }, (_, i) => ring(i / nu, rf(i / nu)));
    for (let i = 0; i < nu; i++) for (let j = 0; j < nt; j++) {
      const a = R[i][j], b = R[i][(j + 1) % nt], c = R[i + 1][j], d = R[i + 1][(j + 1) % nt];
      if (flip) { tri(a, d, b); tri(a, c, d); } else { tri(a, b, d); tri(a, d, c); }
    }
    return R;
  };
  const inner = surface(radius, wall > 0);
  if (wall > 0) {
    const outer = surface((u) => radius(u) + wall, false);
    for (const [i, flip] of [[0, true], [nu, false]]) for (let j = 0; j < nt; j++) {
      const a = inner[i][j], b = inner[i][(j + 1) % nt], c = outer[i][j], d = outer[i][(j + 1) % nt];
      if (flip) { tri(a, b, d); tri(a, d, c); } else { tri(a, d, b); tri(a, c, d); }
    }
  } else {
    for (const [i, flip] of [[0, true], [nu, false]]) {
      const cen = path(i / nu);
      for (let j = 0; j < nt; j++) flip ? tri(cen, inner[i][(j + 1) % nt], inner[i][j]) : tri(cen, inner[i][j], inner[i][(j + 1) % nt]);
    }
  }
  return Float32Array.from(out);
}

/** Binary STL bytes of a triangle soup. */
export function toBinarySTL(tri) {
  const n = tri.length / 9, buf = new ArrayBuffer(84 + 50 * n), dv = new DataView(buf);
  dv.setUint32(80, n, true);
  for (let i = 0; i < n; i++) for (let k = 0; k < 9; k++) dv.setFloat32(84 + 50 * i + 12 + 4 * k, tri[9 * i + k], true);
  return buf;
}
