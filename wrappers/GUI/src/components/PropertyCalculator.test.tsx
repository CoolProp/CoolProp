import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { invoke } from "@tauri-apps/api/core";
import PropertyCalculator from "./PropertyCalculator";

const mockInvoke = vi.mocked(invoke);

function makeReadyState(id = 1) {
  mockInvoke.mockImplementation(async (cmd: string, args?: unknown) => {
    if (cmd === "create_state") return id as never;
    if (cmd === "free_state") return undefined as never;
    if (cmd === "update_state") return undefined as never;
    if (cmd === "get_property") {
      const p = (args as { param: string }).param;
      const fixtures: Record<string, number> = {
        T: 300, P: 101325, Dmass: 996.0, Dmolar: 55340,
        Hmass: 112650, Hmolar: 2030, Smass: 393.5, Smolar: 7.09,
        Umass: 112540, Umolar: 2028, Cvmass: 4130, Cvmolar: 74.4,
        Cpmass: 4180, Cpmolar: 75.3, speed_of_sound: 1500,
        viscosity: 0.000853, conductivity: 0.609, Q: -1, Phase: 0,
      };
      return (fixtures[p] ?? 0) as never;
    }
    throw new Error(`unexpected invoke: ${cmd}`);
  });
}

// The form's <label> elements aren't htmlFor-linked to their inputs, so
// getByLabelText doesn't resolve. Find inputs by their preceding label text
// instead — this also documents the expected DOM structure.
function inputByLabel(label: string): HTMLInputElement {
  const labelEl = screen.getByText(label);
  const input = labelEl.parentElement?.querySelector("input");
  if (!input) throw new Error(`No input found next to label "${label}"`);
  return input as HTMLInputElement;
}
function selectByLabel(label: string): HTMLSelectElement {
  const labelEl = screen.getByText(label);
  const select = labelEl.parentElement?.querySelector("select");
  if (!select) throw new Error(`No select found next to label "${label}"`);
  return select as HTMLSelectElement;
}

beforeEach(() => {
  mockInvoke.mockReset();
  mockInvoke.mockImplementation(async () => undefined as never);
});

describe("PropertyCalculator", () => {
  it("renders the default input pair (P, T) with mass-basis defaults", () => {
    makeReadyState();
    render(<PropertyCalculator fluids={["Water"]} basis="mass" />);
    expect(inputByLabel("P (Pa)").value).toBe("101325");
    expect(inputByLabel("T (K)").value).toBe("300");
  });

  it("switching the input pair updates labels and default values", async () => {
    makeReadyState();
    const user = userEvent.setup();
    render(<PropertyCalculator fluids={["Water"]} basis="mass" />);
    await user.selectOptions(selectByLabel("Input pair"), "h, s");
    expect(inputByLabel("h (J/kg)").value).toBe("100000");
    expect(inputByLabel("s (J/kg·K)").value).toBe("1000");
  });

  it("flipping basis from mass→molar resets to molar defaults", async () => {
    makeReadyState();
    const user = userEvent.setup();
    const { rerender } = render(
      <PropertyCalculator fluids={["Water"]} basis="mass" />
    );
    rerender(<PropertyCalculator fluids={["Water"]} basis="molar" />);
    await user.selectOptions(selectByLabel("Input pair"), "ρ, T");
    expect(inputByLabel("ρ (mol/m³)").value).toBe("55000");
  });

  it("clicking Calculate invokes update_state and renders results", async () => {
    makeReadyState();
    const user = userEvent.setup();
    render(<PropertyCalculator fluids={["Water"]} basis="mass" />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: /calculate/i })).toBeEnabled()
    );
    await user.click(screen.getByRole("button", { name: /calculate/i }));
    await waitFor(() =>
      expect(
        mockInvoke.mock.calls.find((c) => c[0] === "update_state")
      ).toBeDefined()
    );
    await waitFor(() => {
      expect(screen.getAllByText(/^996\./).length).toBeGreaterThan(0);
    });
  });

  it("shows the create_state error when the backend rejects the fluid", async () => {
    mockInvoke.mockImplementation(async (cmd: string) => {
      if (cmd === "create_state") throw new Error("Unknown fluid: Foo");
      return undefined as never;
    });
    render(<PropertyCalculator fluids={["Foo"]} basis="mass" />);
    await waitFor(() => {
      expect(screen.getByText(/Unknown fluid: Foo/)).toBeInTheDocument();
    });
  });
});
