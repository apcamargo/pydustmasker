"""Build runnable Python examples for the Zensical documentation.

The block body contains raw Python::

    /// pyodide-code-block
        session: playground

    import pydustmasker
    print(pydustmasker.__version__)
    ///

Keep the body at the block's own indentation and separate it from the options
with a blank line, because PyMdown Blocks treats indented lines after the
opening fence as options. Blocks sharing a `session` share Python globals in
the browser. Blocks are editable by default, and `editable: false` renders one
readers can run but not modify.

The extension discovers the WebAssembly wheel in `docs/assets/wasm`. The wheel
is the Pyodide build published to PyPI, downloaded by hand after each release.
`pyodide-version.txt` beside it holds the Pyodide version whose ABI the wheel
targets (the `PYODIDE_VERSION` of `.github/workflows/pyodide.yml`), because
`micropip.install` rejects a wheel built for a different runtime.
"""

import sys
import xml.etree.ElementTree as etree
from functools import cache
from pathlib import Path
from typing import Any, cast

# Imported through `pymdownx.blocks` rather than `markdown`, because Pyrefly
# resolves `markdown` to its typeshed stubs here but to the installed sources
# inside the untyped `pymdownx`, and would see two different `Markdown` classes.
from pymdownx.blocks import BlocksExtension, BlocksProcessor, Markdown
from pymdownx.blocks.block import Block, type_any
from pymdownx.superfences import SuperFencesBlockPreprocessor

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override

#: Directory containing the Pyodide wheel downloaded from PyPI.
WHEEL_DIR = Path(__file__).resolve().parents[2] / "docs" / "assets" / "wasm"

#: Pyodide version the wheel beside it was built against.
PYODIDE_VERSION_PATH = WHEEL_DIR / "pyodide-version.txt"


@cache
def _wheel_name() -> str:
    """Return the single Pyodide wheel used by the documentation."""
    wheels = sorted(WHEEL_DIR.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError(
            f"pyodide-code-block needs exactly one wheel in {WHEEL_DIR}, "
            f"found {len(wheels)}. Download the pyemscripten wheel from PyPI."
        )
    return wheels[0].name


@cache
def _pyodide_url() -> str:
    """Return the Pyodide URL recorded beside the wheel."""
    if not PYODIDE_VERSION_PATH.is_file():
        raise ValueError(
            f"pyodide-code-block needs {PYODIDE_VERSION_PATH} with the Pyodide "
            f"version the wheel was built against."
        )
    version = PYODIDE_VERSION_PATH.read_text(encoding="utf-8").strip()
    if not version:
        raise ValueError(f"pyodide-code-block found an empty {PYODIDE_VERSION_PATH}.")
    return f"https://cdn.jsdelivr.net/pyodide/v{version}/full/pyodide.mjs"


class PyodideCodeBlock(Block):
    """Render a Python example the browser runs with Pyodide."""

    NAME = "pyodide-code-block"
    ARGUMENT = False
    OPTIONS = {
        "session": ("default", type_any),
        "editable": (True, type_any),
    }

    @override
    def on_validate(self, parent: etree.Element) -> bool:
        """Validate the block options."""
        if self.options["attrs"]:
            raise ValueError("pyodide-code-block does not support attrs")

        session = self.options["session"]
        if not isinstance(session, str):
            raise TypeError("pyodide-code-block session must be str")
        if not session.strip():
            raise ValueError("pyodide-code-block session must be a non-empty string")

        if not isinstance(self.options["editable"], bool):
            raise TypeError("pyodide-code-block editable must be bool")

        self.wheel = _wheel_name()
        self.pyodide_url = _pyodide_url()
        return True

    @override
    def on_create(self, parent: etree.Element) -> etree.Element:
        """Create the editor, controls, and output panel.

        CodeMirror replaces the source ``div``'s children. The output uses a bare
        ``pre`` so ``pymdownx.highlight`` does not rewrite it as a code block.
        """
        block = etree.SubElement(
            parent,
            "div",
            {
                "class": "pyodide-code-block",
                "data-wheel": self.wheel,
                "data-pyodide-url": self.pyodide_url,
                "data-session": self.options["session"],
                "data-editable": ("true" if self.options["editable"] else "false"),
            },
        )

        frame = etree.SubElement(
            block, "div", {"class": "pyodide-code-block__editor-frame"}
        )
        self.body = etree.SubElement(
            frame,
            "div",
            {
                "class": "pyodide-code-block__editor",
                "aria-label": "Python code",
                "data-source": "",
            },
        )

        controls = etree.SubElement(
            block, "div", {"class": "pyodide-code-block__controls"}
        )
        timer = etree.SubElement(
            controls,
            "span",
            {
                "class": "pyodide-code-block__timer",
                "data-timer": "",
            },
        )
        timer.text = ""
        actions = etree.SubElement(
            controls, "div", {"class": "pyodide-code-block__actions"}
        )
        clear = etree.SubElement(
            actions,
            "button",
            {
                "class": "md-button pyodide-code-block__button",
                "type": "button",
                "disabled": "",
                "data-clear": "",
            },
        )
        clear.text = "Clear"
        run = etree.SubElement(
            actions,
            "button",
            {
                "class": "md-button pyodide-code-block__button pyodide-code-block__button--run",
                "type": "button",
                "disabled": "",
                "data-run": "",
            },
        )
        run.text = "Run"

        output = etree.SubElement(
            block, "div", {"class": "highlight pyodide-code-block__output"}
        )
        etree.SubElement(
            output,
            "pre",
            {"aria-live": "polite", "data-output": ""},
        )
        return block

    @override
    def on_add(self, block: etree.Element) -> etree.Element:
        """Return the source element for Markdown processing."""
        return self.body

    @override
    def on_markdown(self) -> str:
        """Keep the Python source out of the Markdown parser."""
        return "raw"

    @override
    def on_end(self, block: etree.Element) -> None:
        """Render the static source with the site's ordinary Python highlighter."""
        source = self.body.text or ""
        fences = cast(
            SuperFencesBlockPreprocessor, self.md.preprocessors["fenced_code_block"]
        )
        fallback = fences.highlight(
            source,
            "python",
            {},
            self.md,
            classes=[],
            id_value="",
            attrs={},
        )
        self.body.text = self.md.htmlStash.store(fallback)


class PyodideCodeBlockExtension(BlocksExtension):
    """Register the `pyodide-code-block` Block with PyMdown."""

    @override
    def extendMarkdownBlocks(self, md: Markdown, block_mgr: BlocksProcessor) -> None:
        """Register the block implementation."""
        block_mgr.register(PyodideCodeBlock, {})


def makeExtension(**kwargs: Any) -> PyodideCodeBlockExtension:
    """Return the Markdown extension entry point."""
    return PyodideCodeBlockExtension(**kwargs)
