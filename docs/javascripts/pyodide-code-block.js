/**
 * The runnable Python examples built by `docs/scripts/pyodide_code_block.py`.
 *
 * Each block carries a `data-pyodide-url` derived from its wheel build. If the
 * runtime ABI does not match, `micropip.install` rejects the wheel.
 */

import { createNativeTooltip } from "./native-tooltip.js";

// Runtime state. Instant navigation never re-evaluates this module, so one
// interpreter and one editor bundle serve every page.
let pyodidePromise;
let codeEditorPromise;
let installedWheelUrl;
// `setStdout` and `setStderr` are global to that interpreter, so the queue
// has to outlive a navigation too.
let runQueue = Promise.resolve();

// Page state, replaced on every navigation. Each entry of `sessions` holds one
// `session`'s Python globals dict.
const blocks = [];
const sessions = new Map();
let tooltips = [];
let navigationGeneration = 0;

function shortcutKeys(modifierClass, modifierName) {
  const keys = document.createElement("span");
  keys.className = "keys";

  const modifier = document.createElement("kbd");
  modifier.className = modifierClass;
  modifier.textContent = modifierName;
  keys.append(modifier, document.createTextNode(" + "));

  const enter = document.createElement("kbd");
  enter.className = "key-enter";
  enter.textContent = "Enter";
  keys.append(enter);
  return keys;
}

function addRunTooltip(record) {
  const tooltip = createNativeTooltip(document.body, "pyodide-code-block__tooltip");
  record.runButton.setAttribute("aria-describedby", tooltip.id);

  let hovered = false;
  let focused = false;
  const update = () => {
    if (hovered || focused) {
      tooltip.show(record.runButton, (inner) => {
        inner.replaceChildren(
          shortcutKeys("key-control", "Ctrl"),
          document.createTextNode(" or "),
          shortcutKeys("key-command", "Cmd"),
        );
      });
    } else {
      tooltip.hide();
    }
  };

  record.runButton.addEventListener("pointerenter", () => {
    hovered = true;
    update();
  });
  record.runButton.addEventListener("pointerleave", () => {
    hovered = false;
    update();
  });
  record.runButton.addEventListener("focus", () => {
    focused = true;
    update();
  });
  record.runButton.addEventListener("blur", () => {
    focused = false;
    update();
  });

  return tooltip;
}

// Clear a rejected promise so a later run can retry the load.
function loadCodeEditor() {
  if (codeEditorPromise === undefined) {
    codeEditorPromise = import("./pyodide-code-editor.js").catch((error) => {
      codeEditorPromise = undefined;
      throw new Error("The Python editor did not load.", { cause: error });
    });
  }

  return codeEditorPromise;
}

function loadPyodideRuntime(runtimeUrl) {
  if (pyodidePromise === undefined) {
    pyodidePromise = (async () => {
      const { loadPyodide } = await import(runtimeUrl);
      const pyodide = await loadPyodide({ fullStdLib: false });
      await pyodide.loadPackage("micropip");
      return pyodide;
    })().catch((error) => {
      pyodidePromise = undefined;
      throw error;
    });
  }

  return pyodidePromise;
}

function wheelUrl(wheel) {
  return new URL(`../assets/wasm/${wheel}`, import.meta.url).href;
}

async function installWheel(pyodide, wheel) {
  const url = wheelUrl(wheel);
  if (installedWheelUrl === url) return;

  await pyodide.pyimport("micropip").install(url);
  installedWheelUrl = url;
}

function openSession(pyodide, name) {
  let globals = sessions.get(name);
  if (globals === undefined) {
    // `globals.get` returns a new proxy for the `dict` type each time.
    const dictType = pyodide.globals.get("dict");
    try {
      globals = dictType();
    } finally {
      dictType.destroy();
    }
    sessions.set(name, globals);
  }

  return globals;
}

// Queued behind `runQueue` so an in-flight run keeps its WebAssembly globals
// and cannot leak proxies into a later page with the same session name.
function destroySessions() {
  sessions.forEach((globals) => globals.destroy());
  sessions.clear();
}

const reducedMotionMedia = window.matchMedia("(prefers-reduced-motion: reduce)");
const hoverMedia = window.matchMedia("(hover)");

// Make overflowing output keyboard reachable, matching the theme's scrollable
// code blocks. Recalculate this after each render.
function syncOutputScroll(pre) {
  if (hoverMedia.matches && pre.scrollWidth > pre.clientWidth) {
    pre.setAttribute("tabindex", "0");
  } else {
    pre.removeAttribute("tabindex");
  }
}

// Attaching a mouse to a touch device flips `hover` without re-rendering, so
// the decision above has to be revisited rather than baked in at render time.
hoverMedia.addEventListener("change", () => {
  for (const record of blocks) syncOutputScroll(record.outputPre);
});

// Re-adding the class only restarts the keyframes after a forced reflow.
function restartReveal(pre) {
  pre.classList.remove("pyodide-code-block__output--revealed");
  void pre.offsetHeight;
  pre.classList.add("pyodide-code-block__output--revealed");
}

// Replace the output text, animating the panel's height in either direction,
// including a transition to the loading state.
function setOutput(record, text, { reveal = true } = {}) {
  const panel = record.outputPanel;
  const pre = record.outputPre;
  const apply = () => {
    pre.classList.remove("pyodide-code-block__output--stale");
    pre.textContent = text;
    // New text starts at column zero, not the previous run's scroll offset.
    pre.scrollLeft = 0;
    syncOutputScroll(pre);
  };

  // The measure and reflow below are pointless when the panel cannot animate,
  // and would leave an inline height behind.
  if (reducedMotionMedia.matches) {
    if (reveal) pre.classList.remove("pyodide-code-block__output--revealed");
    apply();
    return;
  }

  const startHeight = panel.getBoundingClientRect().height;
  apply();

  // Measure with automatic height so an interrupted transition cannot leave a
  // stale inline target.
  panel.style.height = "auto";
  const endHeight = panel.getBoundingClientRect().height;

  if (Math.abs(startHeight - endHeight) < 1) {
    panel.style.height = "";
    if (reveal) restartReveal(pre);
    return;
  }

  panel.style.height = `${startHeight}px`;
  void panel.offsetHeight;
  panel.style.height = `${endHeight}px`;
  if (reveal) restartReveal(pre);

  const onTransitionEnd = (event) => {
    if (event.target !== panel || event.propertyName !== "height") return;
    panel.removeEventListener("transitionend", onTransitionEnd);
    // Let the panel follow later layout, font, and theme changes.
    panel.style.height = "";
  };
  panel.addEventListener("transitionend", onTransitionEnd);
}

// Skip the loading message when the interpreter and wheel are ready.
// `installedWheelUrl` is set only after the interpreter and wheel both load.
function isEnvironmentPrepared(record) {
  return installedWheelUrl === wheelUrl(record.root.dataset.wheel);
}

// Pyodide's boot holds the main thread, so the message must be painted before
// the load begins. A second frame confirms the first one was.
function nextPaint() {
  return new Promise((resolve) => {
    requestAnimationFrame(() => requestAnimationFrame(resolve));
  });
}

function syncClearButton(record) {
  record.clearButton.disabled = record.outputPre.textContent.trim() === "";
}

function formatElapsed(milliseconds) {
  // Round before the comparison, so a value that rounds up to a whole minute
  // reads as one minute, not sixty seconds.
  const seconds = Number((milliseconds / 1000).toFixed(3));
  if (seconds < 60) return `${seconds.toFixed(3)} s`;

  const minutes = Math.floor(seconds / 60);
  const remainder = seconds - minutes * 60;
  return `${minutes} min ${remainder.toFixed(3)} s`;
}

// `runPythonAsync` blocks intermediate paints, but the final value still uses
// `performance.now()`.
function startTimer(record) {
  const start = performance.now();
  record.timer.textContent = formatElapsed(0);

  // Update once per frame instead of queuing interval writes between paints.
  let handle;
  const tick = () => {
    record.timer.textContent = formatElapsed(performance.now() - start);
    handle = requestAnimationFrame(tick);
  };
  handle = requestAnimationFrame(tick);

  return () => {
    cancelAnimationFrame(handle);
    record.timer.textContent = formatElapsed(performance.now() - start);
  };
}

// Loading the runtime or installing the wheel may outlive the current page.
// Recheck `generation` after each step, not only before the run.
async function run(record, generation) {
  record.outputPre.classList.remove("pyodide-code-block__output--revealed");
  record.runButton.disabled = true;
  record.clearButton.disabled = true;

  // Render all output at once so the panel animates to its final height.
  const lines = [];
  const isPreparing = !isEnvironmentPrepared(record);
  // Mark existing output as stale before replacing it. Preparation replaces the
  // output immediately, so it needs no stale state.
  if (!isPreparing && record.outputPre.textContent.trim() !== "") {
    record.outputPre.classList.add("pyodide-code-block__output--stale");
  }

  const session = record.root.dataset.session;
  try {
    if (isPreparing) {
      // Added before the render so the keyframes start in the painting pass.
      record.outputPre.classList.add("pyodide-code-block__output--loading");
      setOutput(record, "Preparing the Python environment\n", { reveal: false });
      await nextPaint();
      if (generation !== navigationGeneration) return;
    }

    const pyodide = await loadPyodideRuntime(record.root.dataset.pyodideUrl);
    if (generation !== navigationGeneration) return;
    await installWheel(pyodide, record.root.dataset.wheel);
    if (generation !== navigationGeneration) return;
    const globals = openSession(pyodide, session);
    // Drop the message once preparation ends
    if (isPreparing) {
      record.outputPre.classList.remove("pyodide-code-block__output--loading");
      setOutput(record, "", { reveal: false });
    }

    // Started after the boot and wheel install, which the loading message
    // already accounts for, so they stay out of the reported time.
    const stopTimer = startTimer(record);
    pyodide.setStdout({ batched: (value) => lines.push(`${value}\n`) });
    pyodide.setStderr({ batched: (value) => lines.push(`${value}\n`) });
    let result;
    try {
      result = await pyodide.runPythonAsync(record.editor.value(), {
        globals,
      });
    } finally {
      // Stop the timer on both success and failure.
      stopTimer();
      pyodide.setStdout();
      pyodide.setStderr();
    }
    // Python `None` arrives as `undefined`. A non-convertible object arrives as
    // a proxy this code owns. A primitive has no `destroy`, and interpolating
    // either one runs Python `__repr__`, which can raise.
    if (result !== undefined && result !== null) {
      try {
        lines.push(`${result}\n`);
      } finally {
        if (result instanceof pyodide.ffi.PyProxy) result.destroy();
      }
    }
    if (generation === navigationGeneration) setOutput(record, lines.join(""));
  } catch (error) {
    record.outputPre.classList.remove("pyodide-code-block__output--loading");
    // Preserve output produced before the exception. `lines` is empty if the
    // runtime or wheel loading fails.
    lines.push(String(error));
    if (generation === navigationGeneration) setOutput(record, lines.join(""));
  } finally {
    record.isRunning = false;
    // A stale block was detached with the replaced page container. The two
    // writes that do more than set a property are guarded: this one re-enables
    // controls the run queue can still reach, and `setOutput` above would
    // measure and animate a detached panel. Other writes to it are inert.
    if (generation === navigationGeneration) {
      record.runButton.disabled = false;
      syncClearButton(record);
    }
  }
}

// Runs share one Pyodide interpreter, so serialize them instead of interleaving
// their output.
function queueRun(record) {
  // The button is disabled only when `run` starts. `Mod-Enter` can enqueue a
  // run directly, so guard the queue here.
  if (record.isRunning) return;
  record.isRunning = true;

  // Ignore runs queued before navigation. Their editors no longer exist.
  const generation = navigationGeneration;
  runQueue = runQueue
    .catch(() => undefined)
    .then(() => {
      if (generation !== navigationGeneration) {
        record.isRunning = false;
        return undefined;
      }
      return run(record, generation);
    });
}

function createRecord(root) {
  return {
    root,
    source: root.querySelector("[data-source]"),
    outputPanel: root.querySelector(".pyodide-code-block__output"),
    outputPre: root.querySelector("[data-output]"),
    runButton: root.querySelector("[data-run]"),
    clearButton: root.querySelector("[data-clear]"),
    timer: root.querySelector("[data-timer]"),
    editor: null,
    isRunning: false,
  };
}

function initializeBlock(record, codeEditor) {
  record.editor = codeEditor.createPyodideCodeEditor(
    record.source,
    getComputedStyle(record.outputPre),
    () => {
      queueRun(record);
    },
    record.root.dataset.editable !== "false",
  );

  syncClearButton(record);
  // Start with the button disabled in HTML and enable it after the editor is
  // ready.
  record.runButton.disabled = false;
  record.runButton.addEventListener("click", () => {
    queueRun(record);
  });
  // Keep the last run's duration when clearing its output.
  record.clearButton.addEventListener("click", () => {
    setOutput(record, "", { reveal: false });
    syncClearButton(record);
  });
}

async function initializeBlocks(generation) {
  let codeEditor;
  try {
    // Load the CodeMirror bundle only on pages that contain these blocks.
    codeEditor = await loadCodeEditor();
  } catch (error) {
    if (generation === navigationGeneration) {
      for (const record of blocks) {
        record.outputPre.textContent = String(error);
        record.runButton.disabled = true;
      }
    }
    return;
  }

  if (generation !== navigationGeneration) return;
  for (const record of blocks) initializeBlock(record, codeEditor);
}

document$.subscribe(({ body }) => {
  navigationGeneration += 1;
  for (const record of blocks) record.editor?.destroy();
  blocks.length = 0;
  // These hang off `document.body`, which navigation does not replace.
  for (const tooltip of tooltips) tooltip.destroy();
  tooltips = [];
  runQueue = runQueue.catch(() => undefined).then(destroySessions);

  for (const root of body.querySelectorAll(".pyodide-code-block")) {
    blocks.push(createRecord(root));
  }
  if (blocks.length === 0) return;

  tooltips = blocks.map(addRunTooltip);
  initializeBlocks(navigationGeneration);
});
