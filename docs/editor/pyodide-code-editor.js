/* Source for the checked-in bundle `docs/javascripts/pyodide-code-editor.js`.
 *
 * The bundle is committed so that building the site needs no JavaScript
 * toolchain, and no Just recipe rebuilds it. `docs-build` and `docs-serve`
 * must stay free of npm. To rebuild by hand after editing this file:
 *
 *     npm install --no-save \
 *       @codemirror/commands@6.11.1 @codemirror/lang-python@6.2.1 \
 *       @codemirror/language@6.12.4 @codemirror/state@6.7.5 \
 *       @codemirror/view@6.43.12 esbuild@0.25.12
 *     ./node_modules/.bin/esbuild docs/editor/pyodide-code-editor.js \
 *       --bundle --format=esm --minify --legal-comments=inline \
 *       --outfile=docs/javascripts/pyodide-code-editor.js
 *
 * Those pins are what the committed bundle was produced from. Running the
 * command above reproduces it byte for byte, which is the only check that
 * this file and the bundle have not drifted.
 *
 * `--legal-comments=inline` is what carries the CodeMirror notice below into
 * the bundle. CodeMirror is MIT and requires it be retained.
 */

/*! CodeMirror 6 (@codemirror/*, @lezer/*), MIT License
 *
 * Copyright (C) 2018-2025 by Marijn Haverbeke <marijn@haverbeke.nl>,
 * Adrian Heine <mail@adrianheine.de>, and others.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

import { defaultKeymap, history, historyKeymap, indentWithTab } from "@codemirror/commands";
import { pythonLanguage } from "@codemirror/lang-python";
import {
  bracketMatching,
  HighlightStyle,
  indentUnit,
  syntaxHighlighting,
} from "@codemirror/language";
import { EditorState } from "@codemirror/state";
import { EditorView, keymap } from "@codemirror/view";
import { tags } from "@lezer/highlight";

// Map the grammar onto the theme's own highlight variables, so the editor
// follows the site's light and dark palettes with no colours of its own.
const highlightStyle = HighlightStyle.define([
  { tag: tags.comment, color: "var(--md-code-hl-comment-color)" },
  { tag: tags.number, color: "var(--md-code-hl-number-color)" },
  {
    tag: [tags.bool, tags.null, tags.special(tags.variableName)],
    color: "var(--md-code-hl-special-color)",
  },
  {
    tag: [
      tags.keyword,
      tags.controlKeyword,
      tags.definitionKeyword,
      tags.moduleKeyword,
      tags.operatorKeyword,
    ],
    color: "var(--md-code-hl-keyword-color)",
  },
  { tag: tags.string, color: "var(--md-code-hl-string-color)" },
  {
    tag: [
      tags.function(tags.variableName),
      tags.function(tags.definition(tags.variableName)),
      tags.function(tags.propertyName),
    ],
    color: "var(--md-code-hl-function-color)",
  },
  {
    tag: [tags.className, tags.definition(tags.className), tags.namespace, tags.typeName],
    color: "var(--md-code-hl-name-color)",
  },
  { tag: tags.operator, color: "var(--md-code-hl-operator-color)" },
  { tag: tags.punctuation, color: "var(--md-code-hl-punctuation-color)" },
  {
    tag: tags.invalid,
    backgroundColor: "var(--md-code-hl-color)",
    color: "var(--md-code-fg-color)",
  },
]);

// `Mod-Enter` runs the block. It is bound first so it wins over the defaults.
function editorKeymap(onRun) {
  return keymap.of([
    { key: "Mod-Enter", run: () => (onRun(), true) },
    indentWithTab,
    ...historyKeymap,
    ...defaultKeymap,
  ]);
}

export function createPyodideCodeEditor(host, style, onRun, editable = true) {
  // The host holds the Pygments-highlighted static render, whose text is the
  // example's source. Read it before CodeMirror replaces the children, so the
  // block still shows highlighted code when this bundle fails to load.
  let doc = host.textContent.replace(/\n$/, "");
  host.replaceChildren();
  host.style.setProperty("--pyodide-code-block-editor-font-family", style.fontFamily);
  host.style.setProperty("--pyodide-code-block-editor-font-size", style.fontSize);
  // A non-editable block disables typing, pasting, and document-changing
  // commands through `readOnly`. `editable` marks the content as non-editable
  // to the browser. `tabindex` keeps the block keyboard reachable so readers
  // can still focus, select, and copy its code.
  let readOnlyExtensions = editable
    ? []
    : [
        EditorState.readOnly.of(true),
        EditorView.editable.of(false),
        EditorView.contentAttributes.of({ tabindex: "0" }),
      ];
  let state = EditorState.create({
    doc,
    extensions: [
      pythonLanguage,
      history(),
      bracketMatching(),
      indentUnit.of("    "),
      syntaxHighlighting(highlightStyle),
      editorKeymap(onRun),
      EditorView.contentAttributes.of({ "aria-label": host.ariaLabel }),
      ...readOnlyExtensions,
    ],
  });
  let view = new EditorView({ parent: host, state });
  return {
    destroy: () => view.destroy(),
    value: () => view.state.doc.toString(),
  };
}
