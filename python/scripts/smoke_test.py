#!/usr/bin/env python3
"""Smoke test of the repository (see ``btpsd/smoke.py`` for what it checks).

    python3 python/scripts/smoke_test.py [--quick] [--xsim] [--list] [--only name,name]
"""
import sys

from _common import ROOT  # noqa: F401  (puts python/ on sys.path)

from btpsd.smoke import main

if __name__ == "__main__":
    sys.exit(main())
