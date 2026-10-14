/**
 * Zensical's tooltip markup, reused for content that does not fit in `title`.
 * Callers render the inner panel. The container supplies the coordinate
 * system and decides the lifetime: one inside the page container goes away
 * with it on instant navigation, one on `document.body` outlives it, so its
 * owner has to call `destroy`.
 */

let tooltipSequence = 0;

export function createNativeTooltip(container, className) {
  const tooltip = document.createElement("div");
  tooltip.className = `md-tooltip2 ${className}`;
  tooltip.setAttribute("role", "tooltip");
  tooltip.id = `__native_tooltip_${++tooltipSequence}`;
  tooltip.style.setProperty("--md-tooltip-tail", "0px");

  const inner = document.createElement("div");
  inner.className = "md-tooltip2__inner md-typeset";
  tooltip.append(inner);
  container.append(tooltip);

  const show = (host, render) => {
    render(inner);

    // The theme centers the panel using its measured width.
    const panel = inner.getBoundingClientRect();
    tooltip.style.setProperty("--md-tooltip-width", `${panel.width}px`);

    const base = container.getBoundingClientRect();
    const rect = host.getBoundingClientRect();
    tooltip.style.setProperty("--md-tooltip-host-x", `${rect.left - base.left}px`);
    tooltip.style.setProperty("--md-tooltip-host-y", `${rect.top - base.top}px`);
    tooltip.style.setProperty("--md-tooltip-x", `${rect.width / 2}px`);
    tooltip.style.setProperty("--md-tooltip-y", `${8 + rect.height}px`);
    tooltip.classList.remove("md-tooltip2--top");
    tooltip.classList.add("md-tooltip2--bottom", "md-tooltip2--active");
  };

  const hide = () => tooltip.classList.remove("md-tooltip2--active");

  return { id: tooltip.id, show, hide, destroy: () => tooltip.remove() };
}
