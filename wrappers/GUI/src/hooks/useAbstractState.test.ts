import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor, act } from "@testing-library/react";
import { invoke } from "@tauri-apps/api/core";
import { useAbstractState } from "./useAbstractState";

const mockInvoke = vi.mocked(invoke);

beforeEach(() => {
  // mockReset clobbers the default impl from setup.ts; re-install one that
  // succeeds for everything by default. Tests then override per-call.
  mockInvoke.mockReset();
  mockInvoke.mockImplementation(async () => undefined as never);
});

describe("useAbstractState", () => {
  it("creates a state and reports ready", async () => {
    mockInvoke.mockImplementation(async (cmd: string) => {
      if (cmd === "create_state") return 42 as never;
      return undefined as never;
    });
    const { result } = renderHook(() => useAbstractState("HEOS", "Water"));
    await waitFor(() => expect(result.current.ready).toBe(true));
    expect(result.current.id).toBe(42);
    expect(result.current.error).toBeNull();
  });

  it("propagates create_state error", async () => {
    mockInvoke.mockImplementation(async (cmd: string) => {
      if (cmd === "create_state") throw new Error("backend not available");
      return undefined as never;
    });
    const { result } = renderHook(() => useAbstractState("BAD", "Water"));
    await waitFor(() => expect(result.current.error).not.toBeNull());
    expect(result.current.ready).toBe(false);
    expect(result.current.error).toContain("backend not available");
  });

  it("frees the state on unmount", async () => {
    const calls: Array<[string, unknown]> = [];
    mockInvoke.mockImplementation(async (cmd: string, args?: unknown) => {
      calls.push([cmd, args]);
      if (cmd === "create_state") return 7 as never;
      return undefined as never;
    });
    const { result, unmount } = renderHook(() =>
      useAbstractState("HEOS", "Water")
    );
    await waitFor(() => expect(result.current.ready).toBe(true));
    unmount();
    await waitFor(() => {
      expect(calls.some(([c, a]) => c === "free_state" && (a as { id: number }).id === 7)).toBe(true);
    });
  });

  it("update() and get() forward correct arguments", async () => {
    const seen: Array<{ cmd: string; args: unknown }> = [];
    mockInvoke.mockImplementation(async (cmd: string, args?: unknown) => {
      seen.push({ cmd, args });
      if (cmd === "create_state") return 5 as never;
      if (cmd === "get_property") return 996.0 as never;
      return undefined as never;
    });
    const { result } = renderHook(() => useAbstractState("HEOS", "Water"));
    await waitFor(() => expect(result.current.ready).toBe(true));

    await act(async () => {
      await result.current.update("PT_INPUTS", 101325, 300);
    });
    const updateCall = seen.find((c) => c.cmd === "update_state");
    expect(updateCall?.args).toMatchObject({
      id: 5,
      inputPair: "PT_INPUTS",
      v1: 101325,
      v2: 300,
    });

    let rho = 0;
    await act(async () => {
      rho = await result.current.get("Dmass");
    });
    expect(rho).toBe(996.0);
    const getCall = seen.find((c) => c.cmd === "get_property");
    expect(getCall?.args).toMatchObject({ id: 5, param: "Dmass" });
  });

  it("update()/get() throw when not ready", async () => {
    mockInvoke.mockRejectedValueOnce(new Error("nope"));
    const { result } = renderHook(() => useAbstractState("BAD", "Water"));
    await waitFor(() => expect(result.current.error).not.toBeNull());
    await expect(result.current.update("PT_INPUTS", 1, 1)).rejects.toThrow(
      "AbstractState not ready"
    );
    await expect(result.current.get("Dmass")).rejects.toThrow(
      "AbstractState not ready"
    );
  });
});
