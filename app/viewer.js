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
}

export { SLOW_MOTION };
