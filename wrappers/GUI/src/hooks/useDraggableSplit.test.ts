import { describe, it, expect } from "vitest";
import { act, renderHook } from "@testing-library/react";
import { useDraggableSplit } from "./useDraggableSplit";

function fakeMouseDown(clientX: number) {
  return {
    clientX,
    preventDefault: () => {},
  } as unknown as React.MouseEvent;
}

describe("useDraggableSplit", () => {
  it("starts at the initial width", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    expect(result.current.width).toBe(240);
  });

  it("widens when the pointer moves right", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    act(() => result.current.startDrag(fakeMouseDown(100)));
    act(() => {
      window.dispatchEvent(new MouseEvent("mousemove", { clientX: 200 }));
    });
    expect(result.current.width).toBe(340);
    act(() => {
      window.dispatchEvent(new MouseEvent("mouseup"));
    });
  });

  it("clamps to the minimum", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    act(() => result.current.startDrag(fakeMouseDown(500)));
    act(() => {
      window.dispatchEvent(new MouseEvent("mousemove", { clientX: 0 }));
    });
    expect(result.current.width).toBe(180);
    act(() => {
      window.dispatchEvent(new MouseEvent("mouseup"));
    });
  });

  it("clamps to the maximum", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    act(() => result.current.startDrag(fakeMouseDown(0)));
    act(() => {
      window.dispatchEvent(new MouseEvent("mousemove", { clientX: 9999 }));
    });
    expect(result.current.width).toBe(600);
    act(() => {
      window.dispatchEvent(new MouseEvent("mouseup"));
    });
  });

  it("stops responding after mouseup", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    act(() => result.current.startDrag(fakeMouseDown(100)));
    act(() => {
      window.dispatchEvent(new MouseEvent("mousemove", { clientX: 150 }));
    });
    expect(result.current.width).toBe(290);
    act(() => {
      window.dispatchEvent(new MouseEvent("mouseup"));
    });
    act(() => {
      window.dispatchEvent(new MouseEvent("mousemove", { clientX: 999 }));
    });
    // No further change after mouseup
    expect(result.current.width).toBe(290);
  });

  it("clears body cursor / userSelect after drag ends", () => {
    const { result } = renderHook(() => useDraggableSplit(240, 180, 600));
    act(() => result.current.startDrag(fakeMouseDown(0)));
    expect(document.body.style.cursor).toBe("col-resize");
    expect(document.body.style.userSelect).toBe("none");
    act(() => {
      window.dispatchEvent(new MouseEvent("mouseup"));
    });
    expect(document.body.style.cursor).toBe("");
    expect(document.body.style.userSelect).toBe("");
  });
});
