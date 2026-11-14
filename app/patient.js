// Patient-specific geometry viewer: prepared CFD surfaces of a CMHA healthy-control MCA bifurcation
// (Song et al. 2024, CC BY) with imposed lesions and prescribed post-treatment lumens.
// Surfaces come from app/data/patient/<id>/<variant>.stl (binary, faces grouped by patch) plus a JSON
// with patch counts, each M1 wall face's centerline station, and reduced-order profiles along M1.
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { STLLoader } from "three/addons/loaders/STLLoader.js";
import { sampleMap } from "./viewer.js";

const PATCH_COLORS = {
  m1: 0xd04a3b, branch: 0xb8665c, extension: 0x8a97a8, inlet: 0x4c8dff, outlet_1: 0x3ecf8e, outlet_2: 0xf5b041,
};

export class PatientViewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.domElement.className = "patient-canvas";
    this.renderer.domElement.hidden = true;
    container.appendChild(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(35, 1, 0.1, 2000);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x334455, 1.5));
    const key = new THREE.DirectionalLight(0xffffff, 1.3);
    key.position.set(30, 60, 40);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0xffffff, 0.5);
    fill.position.set(-40, -20, -30);
    this.scene.add(fill);
    this.material = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.5, metalness: 0, side: THREE.DoubleSide });
    this.mesh = null;
    this.offset = null;
    this.cache = new Map();
    this.loader = new STLLoader();
    this.active = false;
    new ResizeObserver(() => this._resize()).observe(container);
  }

  async manifest(sid) {
    if (!this._manifest) this._manifest = await (await fetch(`data/patient/${sid}/manifest.json`)).json();
    return this._manifest;
  }

  async load(sid, key) {
    if (!this.cache.has(key)) {
      const [buf, meta] = await Promise.all([
        fetch(`data/patient/${sid}/${key}.stl`).then((r) => { if (!r.ok) throw new Error(`${key}.stl: ${r.status}`); return r.arrayBuffer(); }),
        fetch(`data/patient/${sid}/${key}.json`).then((r) => r.json()),
      ]);
      const geom = this.loader.parse(buf);
      if (!this.offset) {  // centre on the first surface loaded so every variant shares one frame
        geom.computeBoundingBox();
        this.offset = geom.boundingBox.getCenter(new THREE.Vector3()).negate();
      }
      geom.translate(this.offset.x, this.offset.y, this.offset.z);
      geom.computeVertexNormals();
      geom.setAttribute("color", new THREE.BufferAttribute(new Float32Array(geom.attributes.position.count * 3), 3));
      this.cache.set(key, { geom, meta });
    }
    const entry = this.cache.get(key);
    const first = !this.mesh;
    if (this.mesh) this.scene.remove(this.mesh);
    this.mesh = new THREE.Mesh(entry.geom, this.material);
    this.scene.add(this.mesh);
    this.current = entry;
    if (first) this.frame();
    return entry.meta;
  }

  /** mode: "patch" | "diameter" | "wss". Returns the [lo, hi] range used (null for patches). */
  color(mode, range = null) {
    const { geom, meta } = this.current;
    const col = geom.attributes.color.array, c = new THREE.Color();
    const { wall, inlet, outlet_1: o1 } = meta.counts;
    const st = meta.face_station, prof = meta.profile;
    const field = mode === "diameter" ? prof.D_mm : mode === "wss" ? prof.wss_Pa : null;
    let lo = 0, hi = 1;
    if (field) {
      const vals = st.filter((i) => i >= 0).map((i) => field[i]);
      lo = mode === "diameter" ? 0 : 0;
      hi = Math.max(...vals, mode === "wss" ? 5 : 0.5);
      if (range) [lo, hi] = range;
    }
    const nFaces = col.length / 9;
    for (let f = 0; f < nFaces; f++) {
      let hex;
      if (f >= wall) hex = f < wall + inlet ? PATCH_COLORS.inlet : f < wall + inlet + o1 ? PATCH_COLORS.outlet_1 : PATCH_COLORS.outlet_2;
      else if (f >= meta.n_original_wall) hex = PATCH_COLORS.extension;
      else if (st[f] >= 0) {
        if (field) { sampleMap(mode === "wss" ? "inferno" : "viridis", (field[st[f]] - lo) / (hi - lo || 1), c); hex = null; }
        else hex = PATCH_COLORS.m1;
      } else hex = PATCH_COLORS.branch;
      if (hex != null) c.setHex(hex);
      for (let k = 0; k < 3; k++) { col[9 * f + 3 * k] = c.r; col[9 * f + 3 * k + 1] = c.g; col[9 * f + 3 * k + 2] = c.b; }
    }
    geom.attributes.color.needsUpdate = true;
    return field ? [lo, hi] : null;
  }

  frame() {
    const g = this.mesh.geometry;
    g.computeBoundingSphere();
    const sp = g.boundingSphere;
    const d = (sp.radius / Math.sin((this.camera.fov * Math.PI) / 360)) * 1.05;
    this.controls.target.copy(sp.center);
    this.camera.position.copy(sp.center).add(new THREE.Vector3(0.2, 0.55, 1).normalize().multiplyScalar(d));
    this.camera.near = d / 100;
    this.camera.far = d * 10;
    this.camera.updateProjectionMatrix();
    this.controls.update();
  }

  setActive(on) {
    this.active = on;
    this.renderer.domElement.hidden = !on;
    if (on) {
      this._resize();
      this.renderer.setAnimationLoop(() => { this.controls.update(); this.renderer.render(this.scene, this.camera); });
    } else this.renderer.setAnimationLoop(null);
  }

  snapshot() { return this.renderer.domElement.toDataURL("image/png"); }

  _resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / Math.max(h, 1);
    this.camera.updateProjectionMatrix();
  }
}
