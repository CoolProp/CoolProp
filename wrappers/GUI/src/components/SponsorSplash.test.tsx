import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import SponsorSplash from "./SponsorSplash";

const KEY_V7 = "coolprop.sponsorSplash.seen.major.7";
const KEY_V8 = "coolprop.sponsorSplash.seen.major.8";

beforeEach(() => {
  localStorage.clear();
});

describe("SponsorSplash", () => {
  it("renders when no seen-key is set for the current major version", () => {
    render(<SponsorSplash version="7.2.1" />);
    expect(screen.getByText(/Enjoying CoolProp/i)).toBeInTheDocument();
  });

  it("does not render when the current major-version key is set", () => {
    localStorage.setItem(KEY_V7, "1");
    render(<SponsorSplash version="7.2.1" />);
    expect(screen.queryByText(/Enjoying CoolProp/i)).not.toBeInTheDocument();
  });

  it("renders again when only a previous major-version key is set", () => {
    localStorage.setItem(KEY_V7, "1");
    render(<SponsorSplash version="8.0.0" />);
    expect(screen.getByText(/Enjoying CoolProp/i)).toBeInTheDocument();
    expect(localStorage.getItem(KEY_V8)).toBeNull();
  });

  it("dismissing via 'Maybe later' hides it and sets the current major key", async () => {
    const user = userEvent.setup();
    render(<SponsorSplash version="7.2.1" />);
    await user.click(screen.getByRole("button", { name: /Maybe later/i }));
    expect(screen.queryByText(/Enjoying CoolProp/i)).not.toBeInTheDocument();
    expect(localStorage.getItem(KEY_V7)).toBe("1");
  });

  it("the Sponsor link points at the Sponsors URL and sets the seen key", async () => {
    const user = userEvent.setup();
    render(<SponsorSplash version="7.2.1" />);
    const link = screen.getByRole("link", { name: /Sponsor on GitHub/i });
    expect(link).toHaveAttribute("href", "https://github.com/sponsors/CoolProp");
    await user.click(link);
    expect(localStorage.getItem(KEY_V7)).toBe("1");
    // Clicking the CTA also dismisses the splash for the session.
    expect(screen.queryByText(/Enjoying CoolProp/i)).not.toBeInTheDocument();
  });
});
