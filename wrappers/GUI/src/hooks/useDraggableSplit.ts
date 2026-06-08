import { useRef, useState } from "react";

/** Stateful left-pane width with a mousedown handler that drags it within [min, max]. */
export function useDraggableSplit(initial: number, min: number, max: number) {
  const [width, setWidth] = useState(initial);
  const drag = useRef<{ startX: number; startW: number } | null>(null);

  const startDrag = (e: React.MouseEvent) => {
    e.preventDefault();
    drag.current = { startX: e.clientX, startW: width };
    const onMove = (ev: MouseEvent) => {
      if (!drag.current) return;
      const dx = ev.clientX - drag.current.startX;
      setWidth(Math.max(min, Math.min(max, drag.current.startW + dx)));
    };
    const onUp = () => {
      drag.current = null;
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
      document.body.style.cursor = "";
      document.body.style.userSelect = "";
    };
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
    document.body.style.cursor = "col-resize";
    document.body.style.userSelect = "none";
  };

  return { width, startDrag };
}
