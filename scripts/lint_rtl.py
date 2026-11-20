#!/usr/bin/env python3
"""Cheap structural lint of the Verilog sources (no simulator needed).

Checks, per file set given on the command line (default: rtl/*.v and tb/*.v):
  * balanced module/endmodule, function/endfunction, case/endcase, generate/endgenerate, begin/end, (), [], {}
  * every instantiated module is defined, every named port / parameter connection exists, every input port is connected
It does not replace a compiler: run `make -C sim all` with Icarus Verilog, Verilator or Vivado for that.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KW = ("module", "function", "case", "generate")


def strip(src: str) -> str:
    src = re.sub(r"/\*.*?\*/", " ", src, flags=re.S)
    src = re.sub(r"//[^\n]*", " ", src)
    return re.sub(r'"[^"\n]*"', '""', src)


def balanced(text: str, start: int, open_c="(", close_c=")"):
    """Index just past the group starting at text[start] == open_c, and the group content."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == open_c:
            depth += 1
        elif text[i] == close_c:
            depth -= 1
            if depth == 0:
                return i + 1, text[start + 1:i]
    raise ValueError("unbalanced")


def split_top(s: str):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def parse_modules(text: str):
    mods = {}
    for m in re.finditer(r"\bmodule\s+(\w+)\s*", text):
        name, pos = m.group(1), m.end()
        params = []
        if text[pos] == "#":
            pos = text.index("(", pos)
            pos, body = balanced(text, pos)
            for item in split_top(body):
                head = item.split("=")[0]
                params.append(re.findall(r"\w+", head)[-1])
        while text[pos] in " \n\t":
            pos += 1
        ports = []
        if text[pos] == "(":
            pos, body = balanced(text, pos)
            for item in split_top(body):
                ports.append((re.findall(r"\w+", re.sub(r"\[[^\]]*\]", " ", item))[-1],
                              "input" in item.split()))
        mods[name] = {"params": params, "ports": ports}
    return mods


def check_pairs(path: Path, text: str, errs: list):
    toks = re.findall(r"\b(module|endmodule|function|endfunction|case|casez|casex|endcase|generate|endgenerate|begin|end)\b", text)
    c = {k: toks.count(k) for k in set(toks)}
    g = lambda k: c.get(k, 0)  # noqa: E731
    for a, b in (("module", "endmodule"), ("function", "endfunction"), ("generate", "endgenerate")):
        if g(a) != g(b):
            errs.append(f"{path}: {a}/{b} {g(a)} vs {g(b)}")
    if g("case") + g("casez") + g("casex") != g("endcase"):
        errs.append(f"{path}: case/endcase mismatch")
    if g("begin") != g("end"):
        errs.append(f"{path}: begin/end {g('begin')} vs {g('end')}")
    for o, cl in ("()", "[]", "{}"):
        if text.count(o) != text.count(cl):
            errs.append(f"{path}: unbalanced {o}{cl}: {text.count(o)} vs {text.count(cl)}")


def main(argv):
    files = [Path(a) for a in argv] or sorted((ROOT / "rtl").glob("*.v")) + sorted((ROOT / "tb").glob("*.v"))
    texts = {f: strip(f.read_text()) for f in files}
    errs, mods = [], {}
    for f, t in texts.items():
        check_pairs(f, t, errs)
        mods.update(parse_modules(t))
    for f, t in texts.items():
        for m in re.finditer(r"\b(\w+)\s*(#\s*\()?", t):
            name = m.group(1)
            if name not in mods:
                continue
            pos = m.end()
            if m.group(2):
                pos = m.end() - 1
                pos, pbody = balanced(t, pos)
            else:
                pbody = ""
            im = re.match(r"\s*(\w+)\s*\(", t[pos:])
            if not im or im.group(1) in ("input", "output"):
                continue
            inst = im.group(1)
            if t[max(0, m.start() - 7):m.start()].strip().endswith("module"):
                continue
            cstart = pos + im.end() - 1
            _, cbody = balanced(t, cstart)
            conns = set(re.findall(r"\.(\w+)\s*\(", cbody))
            pconns = set(re.findall(r"\.(\w+)\s*\(", pbody))
            spec = mods[name]
            names = {p for p, _ in spec["ports"]}
            for c_ in conns - names:
                errs.append(f"{f.name}: instance {inst} of {name}: unknown port .{c_}")
            for p_ in pconns - set(spec["params"]):
                errs.append(f"{f.name}: instance {inst} of {name}: unknown parameter .{p_}")
            for p, is_in in spec["ports"]:
                if is_in and p not in conns:
                    errs.append(f"{f.name}: instance {inst} of {name}: input {p} not connected")
    for e in errs:
        print("LINT:", e)
    print(f"lint_rtl: {len(files)} files, {len(mods)} modules, {len(errs)} problems")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
