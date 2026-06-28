import { openUrl } from "@tauri-apps/plugin-opener";
import type { MouseEvent, ReactNode } from "react";

interface Props {
  href: string;
  className?: string;
  title?: string;
  children: ReactNode;
  /** Extra handler run before the URL is opened (e.g. to dismiss a modal). */
  onClick?: (e: MouseEvent<HTMLAnchorElement>) => void;
}

/**
 * Anchor for external URLs. In a Tauri webview a bare `<a target="_blank">`
 * does not open the system browser, so we intercept the click and hand the
 * URL to the opener plugin. `href`/`target`/`rel` are kept for accessibility
 * and the right-click "copy link" menu.
 */
export default function ExternalLink({
  href,
  className,
  title,
  children,
  onClick,
}: Props) {
  return (
    <a
      href={href}
      className={className}
      title={title}
      target="_blank"
      rel="noreferrer"
      onClick={(e) => {
        e.preventDefault();
        onClick?.(e);
        void openUrl(href).catch((err) =>
          console.error("Failed to open external URL", href, err),
        );
      }}
    >
      {children}
    </a>
  );
}
