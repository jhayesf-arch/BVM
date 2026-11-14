// Three.js viewport: lumen surface colored by a field, animated flow particles, probe ring.
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const N_THETA = 64;
const N_PARTICLES = 1400;
const SLOW_MOTION = 40; // particles move 40× slower than real blood

// Perceptually-uniform colormaps sampled at 9 stops.
const MAPS = {
  viridis: ["#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c", "#28ae80", "#5ec962", "#addc30", "#fde725"],
  inferno: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#fcffa4"],
  coolwarm: ["#3b4cc0", "#5977e3", "#7b9ff9", "#9ebeff", "#c0d4f5", "#dddcdc", "#f2cbb7", "#f7ac8e", "#ee8468", "#d65244", "#b40426"],
};
const mapCache = {};
function rgbStops(name) {
  if (!mapCache[name]) mapCache[name] = MAPS[name].map((h) => new THREE.Color(h));
  return mapCache[name];
}
export function sampleMap(name, t, out = new THREE.Color()) {
  const stops = rgbStops(name);
  const x = Math.min(Math.max(t, 0), 1) * (stops.length - 1);
  const i = Math.min(Math.floor(x), stops.length - 2);
  return out.copy(stops[i]).lerp(stops[i + 1], x - i);
}
export function cssGradient(name) {
  return `linear-gradient(to top, ${MAPS[name].join(", ")})`;
}

export class Viewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(35, 1, 0.1, 1000);
    this.camera.position.set(4, 14, 30);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x334455, 1.6));
    const key = new THREE.DirectionalLight(0xffffff, 1.4);
    key.position.set(10, 20, 15);
    this.scene.add(key);
    const grid = new THREE.GridHelper(40, 20, 0x3a4658, 0x2a3442);
    grid.position.y = -6;
    this.scene.add(grid);

    this.group = new THREE.Group();
    this.scene.add(this.group);

    // lumen surface (wall), built once and updated in place
    this.wallGeom = new THREE.BufferGeometry();
    this.wallMat = new THREE.MeshStandardMaterial({
      vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.45,
      roughness: 0.55, metalness: 0.0, depthWrite: false,
    });
    this.wall = new THREE.Mesh(this.wallGeom, this.wallMat);
    this.wall.renderOrder = 2;
    this.group.add(this.wall);

    // particles
    this.pGeom = new THREE.BufferGeometry();
    this.pPos = new Float32Array(N_PARTICLES * 3);
    this.pCol = new Float32Array(N_PARTICLES * 3);
    this.pGeom.setAttribute("position", new THREE.BufferAttribute(this.pPos, 3));
    this.pGeom.setAttribute("color", new THREE.BufferAttribute(this.pCol, 3));
    this.particles = new THREE.Points(
      this.pGeom,
      new THREE.PointsMaterial({ size: 0.13, vertexColors: true, sizeAttenuation: true }),
    );
    this.particles.renderOrder = 1;
    this.group.add(this.particles);
    this.seeds = Array.from({ length: N_PARTICLES }, () => this._seed(Math.random()));

    // probe ring
    this.probe = new THREE.Mesh(
      new THREE.TorusGeometry(1, 0.035, 8, 64),
      new THREE.MeshBasicMaterial({ color: 0xffffff }),
    );
    this.probe.renderOrder = 3;
    this.group.add(this.probe);

    this.state = null;
    this.showParticles = true;
    this.paused = false;
    this.clock = new THREE.Clock();
    new ResizeObserver(() => this._resize()).observe(container);
    this._resize();
    this.renderer.setAnimationLoop(() => this._tick());
  }

  _seed(u) {
    // uniform over the cross-section area: radial fraction = sqrt(rand)
    return { u, rho: Math.sqrt(Math.random()) * 0.94, theta: Math.random() * 2 * Math.PI };
  }

  _resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / Math.max(h, 1);
    this.camera.updateProjectionMatrix();
  }

  // Map model coordinates (m; flow along +z, bends in x) to scene (mm; flow along +X).
  _toScene(x, y, z) {
    return [(z - this.zMid) * 1e3, x * 1e3, y * 1e3];
  }

  update(res, field, cmap, range, { reframe = false } = {}) {
    const { geom } = res, n = geom.s.length;
    const { x, y, z } = geom.points;
    this.zMid = (z[0] + z[n - 1]) / 2;
    this.state = { res, field, cmap, range };

    // ring vertices + colors (+ 2 cap centers)
    const nv = n * N_THETA + 2;
    if (!this.wallGeom.getAttribute("position") || this.wallGeom.getAttribute("position").count !== nv) {
      this.wallGeom.setAttribute("position", new THREE.BufferAttribute(new Float32Array(nv * 3), 3));
      this.wallGeom.setAttribute("color", new THREE.BufferAttribute(new Float32Array(nv * 3), 3));
      const idx = [];
      for (let i = 0; i < n - 1; i++)
        for (let j = 0; j < N_THETA; j++) {
          const a = i * N_THETA + j, b = i * N_THETA + ((j + 1) % N_THETA);
          const c = (i + 1) * N_THETA + j, d = (i + 1) * N_THETA + ((j + 1) % N_THETA);
          idx.push(a, b, d, a, d, c);
        }
      const cin = n * N_THETA, cout = cin + 1, last = (n - 1) * N_THETA;
      for (let j = 0; j < N_THETA; j++) {
        idx.push(cin, (j + 1) % N_THETA, j);
        idx.push(cout, last + j, last + ((j + 1) % N_THETA));
      }
      this.wallGeom.setIndex(idx);
    }
    const pos = this.wallGeom.getAttribute("position").array;
    const col = this.wallGeom.getAttribute("color").array;
    const values = this.fieldValues(res, field);
    const color = new THREE.Color();
    for (let i = 0; i < n; i++) {
      const N = geom.N[i], B = geom.B[i], r = res.radius[i];
      sampleMap(cmap, (values[i] - range[0]) / (range[1] - range[0] || 1), color);
      for (let j = 0; j < N_THETA; j++) {
        const th = (2 * Math.PI * j) / N_THETA, c = Math.cos(th), s = Math.sin(th);
        const [sx, sy, sz] = this._toScene(
          x[i] + r * (c * N[0] + s * B[0]), y[i] + r * (c * N[1] + s * B[1]), z[i] + r * (c * N[2] + s * B[2]));
        const k = 3 * (i * N_THETA + j);
        pos[k] = sx; pos[k + 1] = sy; pos[k + 2] = sz;
        col[k] = color.r; col[k + 1] = color.g; col[k + 2] = color.b;
      }
    }
    const ends = [[0, n * N_THETA], [n - 1, n * N_THETA + 1]];
    for (const [i, v] of ends) {
      const [sx, sy, sz] = this._toScene(x[i], y[i], z[i]);
      pos.set([sx, sy, sz], 3 * v);
      sampleMap(cmap, (values[i] - range[0]) / (range[1] - range[0] || 1), color);
      col.set([color.r, color.g, color.b], 3 * v);
    }
    this.wallGeom.getAttribute("position").needsUpdate = true;
    this.wallGeom.getAttribute("color").needsUpdate = true;
    this.wallGeom.computeVertexNormals();
    this.wallGeom.computeBoundingSphere();

    this.vMax = Math.max(...res.velocity) * 2; // centerline peak of the parabolic profile
    // particle size follows the vessel so a 0.4 mm perforator and a 6 mm carotid both read clearly
    this.particles.material.size = Math.min(0.2, Math.max(0.012, 0.13 * (Math.max(...res.radius) * 1e3) / 1.5));
    if (reframe) this.frame();
  }

  fieldValues(res, field) {
    switch (field) {
      case "wss": return res.wss;
      case "pressure": return Array.from(res.pressure, (p) => p / 133.322);
      case "velocity": return res.velocity;
      case "reynolds": return res.reynolds;
      case "strain": return Array.from(res.strain, (e) => e * 100);
      default: return res.wss;
    }
  }

  frame() {
    const sphere = this.wallGeom.boundingSphere;
    if (!sphere) return;
    const d = sphere.radius / Math.sin((this.camera.fov * Math.PI) / 360) * 1.05;
    this.controls.target.copy(sphere.center);
    this.camera.position.copy(sphere.center).add(new THREE.Vector3(0.15, 0.45, 1).normalize().multiplyScalar(d));
    this.controls.update();
  }

  setOpacity(o) { this.wallMat.opacity = o; this.wallMat.depthWrite = o >= 0.99; }
  setParticles(on) { this.showParticles = on; this.particles.visible = on; }

  setProbe(frac) {
    this.probeFrac = frac;
    this._placeProbe();
  }

  _stationAt(u) {
    const n = this.state.res.geom.s.length;
    const f = Math.min(Math.max(u, 0), 1) * (n - 1);
    const i = Math.min(Math.floor(f), n - 2);
    return [i, f - i];
  }

  _placeProbe() {
    if (!this.state || this.probeFrac == null) return;
    const { res } = this.state, g = res.geom;
    const [i] = this._stationAt(this.probeFrac);
    const [sx, sy, sz] = this._toScene(g.points.x[i], g.points.y[i], g.points.z[i]);
    this.probe.position.set(sx, sy, sz);
    const t = g.T[i];
    const dir = new THREE.Vector3(...this._dirToScene(t));
    this.probe.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), dir);
    const r = res.radius[i] * 1e3 * 1.03;
    this.probe.scale.set(r, r, r);
  }

  _dirToScene(v) { return [v[2], v[0], v[1]]; }

  _tick() {
    const dt = Math.min(this.clock.getDelta(), 0.05);
    this.controls.update();
    if (this.state && this.showParticles) this._advanceParticles(this.paused ? 0 : dt);
    this._placeProbe();
    this.renderer.render(this.scene, this.camera);
  }

  _advanceParticles(dt) {
    const { res } = this.state, g = res.geom, L = g.s[g.s.length - 1];
    const color = new THREE.Color();
    const vMax = this.vMax || 1;
    for (let k = 0; k < N_PARTICLES; k++) {
      const p = this.seeds[k];
      let [i, f] = this._stationAt(p.u);
      // Poiseuille: u(r) = 2 V (1 - (r/R)^2)
      const V = res.velocity[i] * (1 - f) + res.velocity[i + 1] * f;
      const speed = 2 * V * (1 - p.rho * p.rho);
      p.u += (speed * dt) / SLOW_MOTION / L;
      if (p.u >= 1) {
        // re-inject at the inlet weighted by Poiseuille flux, pdf ∝ ρ(1-ρ²), which keeps the
        // particle density uniform in volume (the initial seeding) instead of piling up at the wall
        p.u -= 1;
        p.rho = Math.sqrt(1 - Math.sqrt(1 - Math.random())) * 0.94;
        p.theta = Math.random() * 2 * Math.PI;
        [i, f] = this._stationAt(p.u);
      }
      const lerp = (a) => a[i] * (1 - f) + a[i + 1] * f;
      const r = lerp(res.radius) * p.rho, c = Math.cos(p.theta), s = Math.sin(p.theta);
      const N = g.N[i], B = g.B[i];
      const [sx, sy, sz] = this._toScene(
        lerp(g.points.x) + r * (c * N[0] + s * B[0]),
        lerp(g.points.y) + r * (c * N[1] + s * B[1]),
        lerp(g.points.z) + r * (c * N[2] + s * B[2]));
      this.pPos[3 * k] = sx; this.pPos[3 * k + 1] = sy; this.pPos[3 * k + 2] = sz;
      sampleMap("inferno", 0.25 + 0.75 * Math.min(speed / vMax, 1), color);
      this.pCol[3 * k] = color.r; this.pCol[3 * k + 1] = color.g; this.pCol[3 * k + 2] = color.b;
    }
    this.pGeom.getAttribute("position").needsUpdate = true;
    this.pGeom.getAttribute("color").needsUpdate = true;
  }

  snapshot() { return this.renderer.domElement.toDataURL("image/png"); }

  // ───────────── lesion body + device overlays ─────────────
  _frameAtMm(sMm) {
    const { res } = this.state, g = res.geom, n = g.s.length, L = g.s[n - 1];
    const f = Math.min(Math.max(sMm * 1e-3 / L, 0), 1) * (n - 1);
    const i = Math.min(Math.floor(f), n - 2), t = f - i;
    const lerp = (a) => a[i] * (1 - t) + a[i + 1] * t;
    return { p: [lerp(g.points.x), lerp(g.points.y), lerp(g.points.z)], N: g.N[i], B: g.B[i], r: lerp(res.radius), r0: lerp(g.radius),
      rb: lerp(g.base ?? g.radius), i };
  }

  _pt(fr, rad, th, ds = 0) {
    // point at angle th, radius rad (m), optional axial offset ds (m) along the local tangent
    const g = this.state.res.geom, T = g.T[fr.i], c = Math.cos(th), s = Math.sin(th);
    return this._toScene(
      fr.p[0] + rad * (c * fr.N[0] + s * fr.B[0]) + ds * T[0],
      fr.p[1] + rad * (c * fr.N[1] + s * fr.B[1]) + ds * T[1],
      fr.p[2] + rad * (c * fr.N[2] + s * fr.B[2]) + ds * T[2]);
  }

  /** Ring-swept tube between two arc lengths (mm). radiusFn(sMm, frame) -> radius in m. */
  _tube(aMm, bMm, radiusFn, material, { nTheta = 40, dMm = 0.15 } = {}) {
    const steps = Math.max(2, Math.ceil((bMm - aMm) / dMm) + 1), pos = [], idx = [];
    for (let k = 0; k < steps; k++) {
      const sMm = aMm + ((bMm - aMm) * k) / (steps - 1), fr = this._frameAtMm(sMm), rad = radiusFn(sMm, fr);
      for (let j = 0; j < nTheta; j++) pos.push(...this._pt(fr, rad, (2 * Math.PI * j) / nTheta));
    }
    for (let k = 0; k < steps - 1; k++) for (let j = 0; j < nTheta; j++) {
      const a = k * nTheta + j, b = k * nTheta + ((j + 1) % nTheta), c = a + nTheta, d = b + nTheta;
      idx.push(a, b, d, a, d, c);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    geo.setIndex(idx);
    geo.computeVertexNormals();
    return new THREE.Mesh(geo, material);
  }

  _lines(segments, color, opacity = 1) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.Float32BufferAttribute(segments.flat(), 3));
    return new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }));
  }

  /**
   * extras: { lesion: "thrombus" | "icas" | null, overlay: {type, ...} | null }
   * Must be called after update() so that the geometry matches the shown result.
   */
  setExtras({ lesion = null, overlay = null } = {}) {
    if (this.extras) { this.group.remove(this.extras); this.extras.traverse((o) => { o.geometry?.dispose(); o.material?.dispose?.(); }); }
    this.extras = new THREE.Group();
    this.group.add(this.extras);
    if (!this.state) return;
    const { res } = this.state, g = res.geom, n = g.s.length;
    const rNormLoaded = Math.max(res.radius[0], res.radius[n - 1]);

    // lesion body: annulus between the narrowed lumen and the normal wall
    if (lesion) {
      // narrowing is measured against the local healthy radius, so taper and bulbs are not drawn as plaque
      const base = g.base ?? g.radius;
      const inLesion = Array.from(g.radius, (r, i) => r < base[i] * 0.985);
      const i0 = inLesion.indexOf(true), i1 = inLesion.lastIndexOf(true);
      if (i0 >= 0) {
        const aMm = g.s[i0] * 1e3, bMm = g.s[i1] * 1e3;
        const mat = new THREE.MeshStandardMaterial({
          color: lesion === "thrombus" ? 0x8a1020 : 0xd8b45a, roughness: lesion === "thrombus" ? 0.35 : 0.8,
          transparent: true, opacity: 0.9, side: THREE.DoubleSide, depthWrite: true,
        });
        // loaded healthy radius ≈ base radius × local wall stretch (r / r0)
        const healthy = (fr) => fr.rb * (fr.r / fr.r0);
        this.extras.add(this._tube(aMm, bMm, (s, fr) => healthy(fr) * 0.995, mat));
        this.extras.add(this._tube(aMm, bMm, (s, fr) => Math.min(fr.r * 1.015, healthy(fr) * 0.99), mat));
      }
    }

    if (!overlay) return;
    const o = overlay;
    if (o.type === "stent") {
      const segs = [], crowns = 8, pitch = 1.1, amp = 0.38e-3;
      const rings = [];
      for (let s = o.startMm + 0.2; s <= o.endMm - 0.2 + 1e-9; s += pitch) rings.push(s);
      rings.forEach((sMm, k) => {
        const fr = this._frameAtMm(sMm), rad = fr.r * 1.01, pts = [];
        for (let j = 0; j <= 2 * crowns; j++) pts.push(this._pt(fr, rad, (Math.PI * j) / crowns, (j % 2 ? amp : -amp)));
        for (let j = 0; j < 2 * crowns; j++) segs.push(pts[j], pts[j + 1]);
        if (k < rings.length - 1) {
          const fr2 = this._frameAtMm(rings[k + 1]), rad2 = fr2.r * 1.01;
          for (let j = (k % 2) * 2; j < 2 * crowns; j += 4)
            segs.push(this._pt(fr, rad, (Math.PI * (j + 1)) / crowns, amp), this._pt(fr2, rad2, (Math.PI * (j + 1)) / crowns, -amp));
        }
      });
      this.extras.add(this._lines(segs, 0xeef3f8));
      // radiopaque markers at both ends
      for (const sMm of [o.startMm, o.endMm]) {
        const fr = this._frameAtMm(sMm);
        for (let j = 0; j < 4; j++) {
          const m = new THREE.Mesh(new THREE.SphereGeometry(0.09, 10, 8), new THREE.MeshStandardMaterial({ color: 0xffd34d, emissive: 0x332200 }));
          m.position.set(...this._pt(fr, fr.r * 1.02, (Math.PI * j) / 2));
          this.extras.add(m);
        }
      }
    } else if (o.type === "retriever") {
      const a = o.centerMm - o.lengthMm / 2, b = o.centerMm + o.lengthMm / 2, segs = [];
      const rad = Math.min(o.diameterMm / 2 * 1e-3, rNormLoaded * 0.97), strands = 6, pitchMm = 9, dMm = 0.25;
      for (const dir of [1, -1]) for (let k = 0; k < strands; k++) {
        let prev = null;
        for (let s = Math.max(a, 0.3); s <= Math.min(b, (g.s[n - 1] * 1e3) - 0.3); s += dMm) {
          const fr = this._frameAtMm(s), th = (2 * Math.PI * k) / strands + dir * (2 * Math.PI * (s - a)) / pitchMm;
          const p = this._pt(fr, rad, th);
          if (prev) segs.push(prev, p);
          prev = p;
        }
      }
      this.extras.add(this._lines(segs, 0xffb347));
      // push wire along the centerline from the inlet to the device
      const wire = [];
      for (let s = 0; s < Math.max(a, 0.3); s += 0.5) wire.push(this._pt(this._frameAtMm(s), 0, 0), this._pt(this._frameAtMm(Math.min(s + 0.5, a)), 0, 0));
      this.extras.add(this._lines(wire, 0xcfd6df));
    } else if (o.type === "balloon") {
      const a = Math.max(o.centerMm - o.lengthMm / 2, 0.2), b = Math.min(o.centerMm + o.lengthMm / 2, g.s[n - 1] * 1e3 - 0.2);
      const rb = o.diameterMm / 2 * 1e-3, taper = 1.0;
      const mat = new THREE.MeshStandardMaterial({ color: 0x5fd3ff, transparent: true, opacity: 0.35, roughness: 0.2, depthWrite: false, side: THREE.DoubleSide });
      this.extras.add(this._tube(a, b, (s) => rb * Math.max(0.15, Math.min(1, (s - a) / taper, (b - s) / taper)), mat));
    } else if (o.type === "catheter") {
      const ro = o.odMm / 2 * 1e-3, ri = o.idMm / 2 * 1e-3;
      const mat = new THREE.MeshStandardMaterial({ color: 0xc9d3df, roughness: 0.5, side: THREE.DoubleSide });
      this.extras.add(this._tube(0, o.tipMm, () => ro, mat, { nTheta: 32 }));
      this.extras.add(this._tube(0, o.tipMm, () => ri, new THREE.MeshStandardMaterial({ color: 0x3a4250, side: THREE.DoubleSide }), { nTheta: 32 }));
      const fr = this._frameAtMm(o.tipMm), segs = [];
      for (let j = 0; j < 32; j++) segs.push(this._pt(fr, ri, (2 * Math.PI * j) / 32), this._pt(fr, ro, (2 * Math.PI * j) / 32));
      this.extras.add(this._lines(segs, 0x9aa6b5));
    }
  }
}

export { SLOW_MOTION };
