/**
 * Fitting for user-defined incompressible (pure) fluids.
 *
 * Mirrors the offline Python pipeline (dev/incompressible_liquids/CPIncomp/
 * ChebyshevFits.py): density and heat capacity are fitted as Chebyshev
 * series in T on the data range, with the degree chosen by generalized
 * cross-validation so rounded near-linear tables get low order instead of
 * ringing between the points. Conductivity is a centered polynomial and
 * viscosity an exp-polynomial (a polynomial fit of ln(mu)), matching the
 * fit types the C++ backend supports for those properties.
 *
 * The output of buildFluidJson is accepted verbatim by
 * add_fluids_as_JSON("INCOMP", ...) — the same schema as the fluid files
 * under dev/incompressible_liquids/json/.
 */

export interface ParsedTable {
  T: number[];
  columns: Record<PropertyKey_, (number | null)[]>;
  warnings: string[];
}

export type PropertyKey_ = "density" | "specific_heat" | "conductivity" | "viscosity";

export interface FitReportEntry {
  points: number;
  degree: number;
  rms: number; // absolute RMS in the property's SI unit
  nrms: number; // RMS / (max - min) of the data
}

export interface FluidFitResult {
  fluidJson: Record<string, unknown>;
  report: Partial<Record<PropertyKey_, FitReportEntry>>;
  warnings: string[];
}

const HEADER_ALIASES: Record<string, PropertyKey_ | "T"> = {
  t: "T", temp: "T", temperature: "T",
  rho: "density", d: "density", density: "density",
  cp: "specific_heat", c: "specific_heat", heat: "specific_heat", specific_heat: "specific_heat",
  k: "conductivity", lambda: "conductivity", l: "conductivity", conductivity: "conductivity",
  mu: "viscosity", v: "viscosity", eta: "viscosity", visc: "viscosity", viscosity: "viscosity",
};

const DEFAULT_COLUMN_ORDER: (PropertyKey_ | "T")[] = ["T", "density", "specific_heat", "conductivity", "viscosity"];

/** Parse a pasted table (TSV/CSV/whitespace). First row may be a header
 * naming the columns (T, rho, cp, k, mu ..., optionally with units in
 * parentheses); without one, columns are assumed in the order
 * T, density, cp, conductivity, viscosity. All values SI: K, kg/m3,
 * J/kg/K, W/m/K, Pa.s. Blank cells = no data.
 *
 * Structural problems throw rather than warn: a header cell mapping to two
 * properties or a row longer than the header would otherwise assign values
 * to the wrong property and fit a silently wrong fluid. */
export function parseTable(text: string): ParsedTable {
  const warnings: string[] = [];
  // One delimiter for the whole table: a spreadsheet paste is tab-separated
  // (so "T (K)" stays one cell); otherwise fall back to ; then , then runs
  // of spaces. Never mix -- mixing is how header cells with units get split
  // into stray one-letter tokens that alias to the wrong property.
  const delimiter = text.includes("\t") ? /\t/ : text.includes(";") ? /;/ : text.includes(",") ? /,/ : / +/;
  const rows = text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0)
    .map((line) => line.split(delimiter).map((cell) => cell.trim()));
  if (rows.length === 0) return { T: [], columns: {} as ParsedTable["columns"], warnings: ["no data"] };

  let mapping: (PropertyKey_ | "T" | null)[];
  let dataRows = rows;
  const firstRowNumeric = rows[0].every((c) => c === "" || Number.isFinite(parseFloat(c)));
  if (!firstRowNumeric) {
    mapping = rows[0].map((h) => {
      // strip parenthesized units ("rho (kg/m3)") before alias lookup
      const bare = h.replace(/\(.*?\)/g, "").trim().toLowerCase().replace(/[^a-z_]/g, "");
      return HEADER_ALIASES[bare] ?? null;
    });
    rows[0].forEach((h, i) => {
      if (mapping[i] === null && h !== "") warnings.push(`column "${h}" not recognised, ignored`);
    });
    if (!mapping.includes("T")) throw new Error("header row found but no temperature column recognised");
    const seen = new Set<string>();
    for (const key of mapping) {
      if (!key) continue;
      if (seen.has(key)) throw new Error(`two header columns map to the same property (${key})`);
      seen.add(key);
    }
    dataRows = rows.slice(1);
  } else {
    mapping = DEFAULT_COLUMN_ORDER.slice(0, rows[0].length);
  }

  const T: number[] = [];
  const columns: Record<string, (number | null)[]> = {};
  for (const key of mapping) if (key && key !== "T") columns[key] = [];

  for (const row of dataRows) {
    if (row.length > mapping.length) {
      throw new Error(`row "${row.join(" ")}" has ${row.length} cells but the header has ${mapping.length} columns`);
    }
    const tIdx = mapping.indexOf("T");
    const tVal = parseFloat(row[tIdx] ?? "");
    if (!Number.isFinite(tVal)) {
      warnings.push(`skipped row without temperature: "${row.join(" ")}"`);
      continue;
    }
    T.push(tVal);
    mapping.forEach((key, i) => {
      if (!key || key === "T") return;
      const v = parseFloat(row[i] ?? "");
      columns[key].push(Number.isFinite(v) ? v : null);
    });
  }

  // sort by temperature (mirrors the offline loaders, DATA_AUDIT.md)
  const order = T.map((_, i) => i).sort((a, b) => T[a] - T[b]);
  const sortedT = order.map((i) => T[i]);
  for (const key of Object.keys(columns)) columns[key] = order.map((i) => columns[key][i]);
  for (let i = 1; i < sortedT.length; i++) {
    if (sortedT[i] === sortedT[i - 1]) warnings.push(`duplicate temperature ${sortedT[i]} K`);
  }
  if (sortedT.length > 0 && sortedT[0] <= 0) warnings.push("temperatures must be in Kelvin (found T <= 0)");

  return { T: sortedT, columns: columns as ParsedTable["columns"], warnings };
}

/** Solve the linear least-squares system design*x ~= y via normal equations
 * with Gaussian elimination (partial pivoting). The design matrices here are
 * small (<= 9 columns) and Chebyshev-based, hence well-conditioned. */
function lstsq(design: number[][], y: number[]): number[] {
  const n = design[0].length;
  const A: number[][] = Array.from({ length: n }, () => new Array<number>(n + 1).fill(0));
  for (let r = 0; r < design.length; r++) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) A[i][j] += design[r][i] * design[r][j];
      A[i][n] += design[r][i] * y[r];
    }
  }
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(A[r][col]) > Math.abs(A[pivot][col])) pivot = r;
    [A[col], A[pivot]] = [A[pivot], A[col]];
    if (Math.abs(A[col][col]) < 1e-300) throw new Error("singular fit matrix");
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = A[r][col] / A[col][col];
      for (let c = col; c <= n; c++) A[r][c] -= f * A[col][c];
    }
  }
  return A.map((row, i) => row[n] / A[i][i]);
}

function chebRow(u: number, degree: number): number[] {
  const row = new Array<number>(degree + 1);
  row[0] = 1;
  if (degree >= 1) row[1] = u;
  for (let k = 2; k <= degree; k++) row[k] = 2 * u * row[k - 1] - row[k - 2];
  return row;
}

/** Clenshaw evaluation of a Chebyshev series on [Tmin, Tmax]. */
export function chebEval(coeffs: number[], T: number, Trange: [number, number]): number {
  const u = (2 * T - (Trange[1] + Trange[0])) / (Trange[1] - Trange[0]);
  let ukp1 = 0;
  let ukp2 = 0;
  for (let k = coeffs.length - 1; k > 0; k--) {
    const uk = 2 * u * ukp1 - ukp2 + coeffs[k];
    ukp2 = ukp1;
    ukp1 = uk;
  }
  return coeffs[0] + u * ukp1 - ukp2;
}

interface Fit1D {
  coeffs: number[];
  degree: number;
  rms: number;
  nrms: number;
}

/** Chebyshev fit in T with GCV degree selection, mirroring
 * ChebyshevFits.fit_from_data (1D pure-fluid case). */
export function chebFitGCV(T: number[], y: (number | null)[], Trange: [number, number], maxDegree = 8): Fit1D | null {
  const pts: { u: number; v: number }[] = [];
  for (let i = 0; i < T.length; i++) {
    const v = y[i];
    if (v !== null && Number.isFinite(v)) {
      pts.push({ u: (2 * T[i] - (Trange[1] + Trange[0])) / (Trange[1] - Trange[0]), v });
    }
  }
  if (pts.length < 3) return null;
  const values = pts.map((p) => p.v);
  const spread = Math.max(...values) - Math.min(...values);

  let best: { gcv: number; fit: Fit1D } | null = null;
  const cap = Math.min(maxDegree, pts.length - 1);
  for (let degree = 1; degree <= cap; degree++) {
    if (degree + 1 >= pts.length) break; // need leftover DOF for validation
    const design = pts.map((p) => chebRow(p.u, degree));
    let coeffs: number[];
    try {
      coeffs = lstsq(design, values);
    } catch {
      continue;
    }
    let rss = 0;
    for (let i = 0; i < pts.length; i++) {
      let pred = 0;
      for (let k = 0; k <= degree; k++) pred += design[i][k] * coeffs[k];
      rss += (pred - values[i]) ** 2;
    }
    const gcv = (pts.length * rss) / (pts.length - (degree + 1)) ** 2;
    if (best === null || gcv < best.gcv) {
      const rms = Math.sqrt(rss / pts.length);
      best = { gcv, fit: { coeffs, degree, rms, nrms: spread > 0 ? rms / spread : rms } };
    }
  }
  if (best === null) {
    // very short tables: exactly-determined lowest order
    const degree = Math.max(1, Math.min(cap, pts.length - 1));
    const design = pts.map((p) => chebRow(p.u, degree));
    const coeffs = lstsq(design, values);
    return { coeffs, degree, rms: 0, nrms: 0 };
  }
  return best.fit;
}

/** Centered monomial fit (for conductivity / ln-viscosity), GCV-selected. */
function polyFitGCV(T: number[], y: (number | null)[], Tbase: number, maxDegree = 4): Fit1D | null {
  const pts: { dT: number; v: number }[] = [];
  for (let i = 0; i < T.length; i++) {
    const v = y[i];
    if (v !== null && Number.isFinite(v)) pts.push({ dT: T[i] - Tbase, v });
  }
  if (pts.length < 3) return null;
  const values = pts.map((p) => p.v);
  const spread = Math.max(...values) - Math.min(...values);
  // scale (T - Tbase) to O(1) for conditioning, then rescale coefficients
  const scale = Math.max(...pts.map((p) => Math.abs(p.dT))) || 1;

  let best: { gcv: number; fit: Fit1D } | null = null;
  for (let degree = 1; degree <= Math.min(maxDegree, pts.length - 2); degree++) {
    const design = pts.map((p) => {
      const row = new Array<number>(degree + 1);
      row[0] = 1;
      for (let k = 1; k <= degree; k++) row[k] = row[k - 1] * (p.dT / scale);
      return row;
    });
    let scaled: number[];
    try {
      scaled = lstsq(design, values);
    } catch {
      continue;
    }
    let rss = 0;
    for (let i = 0; i < pts.length; i++) {
      let pred = 0;
      for (let k = 0; k <= degree; k++) pred += design[i][k] * scaled[k];
      rss += (pred - values[i]) ** 2;
    }
    const gcv = (pts.length * rss) / (pts.length - (degree + 1)) ** 2;
    if (best === null || gcv < best.gcv) {
      const coeffs = scaled.map((c, k) => c / scale ** k);
      const rms = Math.sqrt(rss / pts.length);
      best = { gcv, fit: { coeffs, degree, rms, nrms: spread > 0 ? rms / spread : rms } };
    }
  }
  return best?.fit ?? null;
}

export interface FluidDefinition {
  name: string;
  description: string;
  reference: string;
  table: ParsedTable;
}

/** Fit all supplied properties and assemble the INCOMP JSON definition.
 * Density and heat capacity are required; conductivity/viscosity optional. */
export function fitFluid(def: FluidDefinition): FluidFitResult {
  const { T, columns } = def.table;
  const warnings = [...def.table.warnings];
  if (!/^[A-Za-z][A-Za-z0-9_-]*$/.test(def.name)) {
    throw new Error("fluid name must start with a letter and contain only letters, digits, '-' or '_'");
  }
  if (T.length < 3) throw new Error("need at least 3 temperature points");
  const Trange: [number, number] = [T[0], T[T.length - 1]];
  if (!(Trange[0] > 0 && Trange[1] > Trange[0])) throw new Error("temperatures must be positive Kelvin, increasing");
  const Tbase = 0.5 * (Trange[0] + Trange[1]);

  const report: Partial<Record<PropertyKey_, FitReportEntry>> = {};
  const fluidJson: Record<string, unknown> = {
    name: def.name,
    description: def.description || def.name,
    reference: def.reference || "user-defined (CoolProp GUI)",
    Tmin: Trange[0],
    Tmax: Trange[1],
    TminPsat: Trange[1], // no vapour-pressure data: psat unavailable rather than wrong
    Tbase,
    xbase: 0.0,
    xid: "pure",
  };

  const caloric: PropertyKey_[] = ["density", "specific_heat"];
  for (const prop of caloric) {
    const fit = chebFitGCV(T, columns[prop] ?? [], Trange);
    if (fit === null) throw new Error(`${prop} needs at least 3 numeric values`);
    if (fit.coeffs.some((c) => !Number.isFinite(c))) throw new Error(`${prop} fit failed (non-finite coefficients)`);
    // guard against a fit swinging non-positive anywhere in range
    for (let i = 0; i <= 40; i++) {
      const Tq = Trange[0] + ((Trange[1] - Trange[0]) * i) / 40;
      if (chebEval(fit.coeffs, Tq, Trange) <= 0) throw new Error(`${prop} fit is not positive over the whole range`);
    }
    fluidJson[prop + "_cheb"] = {
      type: "chebyshev",
      Trange: [Trange[0], Trange[1]],
      xbase: 0.0,
      coeffs: fit.coeffs.map((c) => [c]),
      NRMS: fit.nrms,
      fit_source: "tabular_data",
    };
    const n = (columns[prop] ?? []).filter((v) => v !== null).length;
    report[prop] = { points: n, degree: fit.degree, rms: fit.rms, nrms: fit.nrms };
  }

  const cond = polyFitGCV(T, columns.conductivity ?? [], Tbase);
  if (cond !== null) {
    fluidJson.conductivity = { type: "polynomial", coeffs: cond.coeffs.map((c) => [c]), NRMS: cond.nrms };
    report.conductivity = {
      points: (columns.conductivity ?? []).filter((v) => v !== null).length,
      degree: cond.degree,
      rms: cond.rms,
      nrms: cond.nrms,
    };
  } else if ((columns.conductivity ?? []).some((v) => v !== null)) {
    warnings.push("conductivity: fewer than 3 values, not fitted");
  }

  const muColumn = columns.viscosity ?? [];
  if (muColumn.some((v) => v !== null && v <= 0)) {
    warnings.push("viscosity: non-positive values ignored");
  }
  const lnMu = muColumn.map((v) => (v !== null && v > 0 ? Math.log(v) : null));
  const visc = polyFitGCV(T, lnMu, Tbase);
  if (visc !== null) {
    const rmsAbs = visc.rms; // RMS of ln(mu): relative error scale for mu itself
    fluidJson.viscosity = { type: "exppolynomial", coeffs: visc.coeffs.map((c) => [c]), NRMS: visc.nrms };
    report.viscosity = {
      points: lnMu.filter((v) => v !== null).length,
      degree: visc.degree,
      rms: rmsAbs,
      nrms: visc.nrms,
    };
  } else if (muColumn.some((v) => v !== null)) {
    warnings.push("viscosity: fewer than 3 usable values, not fitted");
  }

  return { fluidJson, report, warnings };
}
