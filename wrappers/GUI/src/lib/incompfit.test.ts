import { describe, expect, it } from "vitest";
import { chebEval, chebFitGCV, fitFluid, parseTable } from "./incompfit";

// The ExamplePure template fluid's tabulated data (Therminol D12,
// dev/incompressible_liquids/CPIncomp/ExampleObjects.py) -- the same table
// the offline Python pipeline fits, so the results can be cross-checked
// against the committed json/ExamplePure.json entries.
const EXAMPLE_PURE = [
  "T\trho\tcp\tk\tmu",
  ...[
    [323.15, 740, 2235, 0.105, 0.804],
    [333.15, 733, 2280, 0.104, 0.704],
    [343.15, 726, 2326, 0.102, 0.623],
    [353.15, 717, 2361, 0.1, 0.556],
    [363.15, 710, 2406, 0.098, 0.498],
    [373.15, 702, 2445, 0.096, 0.451],
    [383.15, 695, 2485, 0.095, 0.41],
    [393.15, 687, 2528, 0.093, 0.374],
    [403.15, 679, 2571, 0.091, 0.346],
    [413.15, 670, 2607, 0.089, 0.317],
    [423.15, 662, 2645, 0.087, 0.289],
  ].map((row) => row.join("\t")),
].join("\n");

describe("parseTable", () => {
  it("recognises a header row and maps aliases", () => {
    const t = parseTable("temperature,density,CP\n300,998,4184\n310,994,4187");
    expect(t.T).toEqual([300, 310]);
    expect(t.columns.density).toEqual([998, 994]);
    expect(t.columns.specific_heat).toEqual([4184, 4187]);
  });

  it("assumes T,rho,cp,k,mu order without a header", () => {
    const t = parseTable("300 998 4184\n310 994 4187");
    expect(t.T).toEqual([300, 310]);
    expect(t.columns.density).toEqual([998, 994]);
  });

  it("sorts by temperature and keeps blank cells as null", () => {
    const t = parseTable("T\trho\tcp\n310\t994\t\n300\t998\t4184");
    expect(t.T).toEqual([300, 310]);
    expect(t.columns.density).toEqual([998, 994]);
    expect(t.columns.specific_heat).toEqual([4184, null]);
  });
});

describe("chebFitGCV", () => {
  it("reproduces exactly Chebyshev data", () => {
    // rho = 1000 - 5*u on [280, 360]: coefficients [1000, -5]
    const T = [280, 300, 320, 340, 360];
    const rho = T.map((t) => 1000 - 5 * ((2 * t - 640) / 80));
    const fit = chebFitGCV(T, rho, [280, 360]);
    expect(fit).not.toBeNull();
    expect(fit!.degree).toBe(1);
    expect(fit!.coeffs[0]).toBeCloseTo(1000, 9);
    expect(fit!.coeffs[1]).toBeCloseTo(-5, 9);
  });

  it("keeps the order low for noisy near-linear data (GCV)", () => {
    const T = Array.from({ length: 11 }, (_, i) => 300 + 10 * i);
    // linear + rounding-level noise
    const y = T.map((t, i) => 1000 - 0.7 * (t - 300) + (i % 2 === 0 ? 0.3 : -0.3));
    const fit = chebFitGCV(T, y, [300, 400]);
    expect(fit).not.toBeNull();
    expect(fit!.degree).toBeLessThanOrEqual(3);
  });
});

describe("fitFluid", () => {
  it("fits the ExamplePure table consistently with the offline pipeline", () => {
    const table = parseTable(EXAMPLE_PURE);
    const { fluidJson, report } = fitFluid({
      name: "ExamplePureTest",
      description: "",
      reference: "",
      table,
    });
    expect(fluidJson.Tmin).toBe(323.15);
    expect(fluidJson.Tmax).toBe(423.15);
    expect(fluidJson.xid).toBe("pure");

    // Cross-check: the offline Python fit (GCV-selected Chebyshev) of the
    // same data evaluates to 736.46 kg/m3 at 328.15 K, matching the
    // committed ExamplePure golden value used in the C++ self-tests.
    const dens = fluidJson.density_cheb as { coeffs: number[][]; Trange: [number, number] };
    const rho = chebEval(dens.coeffs.map((r) => r[0]), 328.15, dens.Trange);
    expect(rho).toBeCloseTo(736.46, 1);

    // all four properties fitted, with sane relative errors
    for (const key of ["density", "specific_heat", "conductivity", "viscosity"] as const) {
      expect(report[key], key).toBeDefined();
      expect(report[key]!.nrms).toBeLessThan(0.05);
    }
    // viscosity ships as exppolynomial (positive by construction)
    expect((fluidJson.viscosity as { type: string }).type).toBe("exppolynomial");
    expect((fluidJson.conductivity as { type: string }).type).toBe("polynomial");
  });

  it("rejects tables without enough caloric data", () => {
    const table = parseTable("T\trho\n300\t998\n310\t994\n320\t990");
    expect(() => fitFluid({ name: "X", description: "", reference: "", table })).toThrow(/specific_heat/);
  });

  it("rejects invalid fluid names", () => {
    const table = parseTable(EXAMPLE_PURE);
    expect(() => fitFluid({ name: "bad name!", description: "", reference: "", table })).toThrow(/name/);
  });

  it("rejects fits that swing non-positive", () => {
    // density data crossing zero: every fit of it crosses zero too
    const rows = ["T\trho\tcp", "300\t100\t4000", "310\t50\t4000", "320\t1\t4000", "330\t-50\t4000", "340\t-100\t4000"];
    expect(() => fitFluid({ name: "X", description: "", reference: "", table: parseTable(rows.join("\n")) })).toThrow(/positive/);
  });
});
