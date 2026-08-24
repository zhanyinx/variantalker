# MAFigate Makefile
# Development automation for MAFigate application

.PHONY: help install test test-cov test-fast parity-reference app-load-check lint format clean run build-mac build-win deps-check qa dev build check-vendor sync-vendor check-params

# Default target
help:
	@echo "MAFigate Development Commands"
	@echo "============================"
	@echo ""
	@echo "Setup Commands:"
	@echo "  install     Install dependencies"
	@echo ""
	@echo "Development Commands:"
	@echo "  run         Start the Streamlit application"
	@echo "  test        Run all tests"
	@echo "  test-cov    Run tests with coverage, measured in a tree with no bin/"
	@echo "  parity-reference  Parity sweep over the full GERSOM reference (needs PARITY_REFERENCE_ROOT)"
	@echo "  app-load-check    Load a MAF through both of the running app's load paths"
	@echo "  lint        Run code linting"
	@echo "  format      Format code with black and isort"
	@echo ""
	@echo "Vendored Pipeline Code (see vendor/README.md):"
	@echo "  check-vendor  Verify vendor/ still matches bin/ (gates the installer builds)"
	@echo "  sync-vendor   Re-copy bin/ into vendor/ to repair drift"
	@echo "  check-params  Verify config/pipeline_params.py still matches nextflow.config"
	@echo "  check-deps    Verify this interpreter has what requirements.txt declares"
	@echo ""
	@echo "Build Commands (see build/BUILD_INSTRUCTIONS.md):"
	@echo "  build-mac   Build macOS .dmg installer (macOS only)"
	@echo "  build-win   Build Windows .exe installer (Windows only)"
	@echo "  version     Print APP_VERSION, the one place the version is written"
	@echo "  release-tag Print the tag a release is cut from (mafigate-v<version>)"
	@echo "  release-preflight  Refuse to publish a draft the tree has moved on from"
	@echo ""
	@echo "Maintenance Commands:"
	@echo "  clean       Clean temporary files and caches"
	@echo "  deps-check  Check for dependency updates"

# Installation commands
#
# Through $(PYTHON) -m pip rather than a bare `pip`, and for the reason issue #162 found
# the hard way: `pip` and `streamlit` are separate PATH entries and need not belong to the
# same python. On the machine where #162 was measured they did not, so `make install`
# filled miniconda while `make run` launched Homebrew's streamlit — every declared
# dependency installed, none of it where the app ran, and `plotly` therefore absent
# wherever it looked. $(PYTHON) is defined below, beside check-vendor, and one variable
# for installing and launching is the whole point: it cannot disagree with itself.
# It is `./setup.sh` that *builds* the environment, not this target: make expands $(PYTHON)
# before any recipe runs, so a target that created .venv could not then install into it. So
# this fills whatever $(PYTHON) resolved to — the venv when there is one — and says out loud
# when that is the ambient interpreter, which is the install issue #258 exists to stop being
# the silent default. It warns rather than refuses: a contributor who has deliberately
# arranged their own interpreter is not doing anything wrong.
install:
	@echo "📦 Installing dependencies into: $(PYTHON)"
	@if [ "$(PYTHON)" = "python3" ]; then \
		echo "⚠️  That is whatever python3 is first on PATH, and this fills it directly."; \
		echo "   ./setup.sh builds an isolated environment instead."; \
	fi
	@$(PYTHON) -m pip install -r requirements.txt
	@echo "🔍 Checking that what requirements.txt asked for arrived..."
	@$(PYTHON) check_dependencies.py

# Vendored pipeline code (see vendor/README.md)
#
# vendor/ holds copies of the pipeline's filter code so the app reaches the same
# PASS/NOPASS verdict as bin/filter_variants.py. A copy is only worth anything while
# it is provably still a copy, so check-vendor gates both installer builds: a .dmg or
# .exe must never be cut from a drifted copy, because the drift would be invisible in
# the shipped app (which carries no bin/ to compare against).
#
# The check is deliberately stdlib-only — no pytest, no pandas — so it runs on a bare
# build machine. tests/test_vendor_drift.py asserts it agrees with the pytest guard.

# The one interpreter every recipe below goes through (issue #162), and since issue #258 it
# prefers the virtual environment ./setup.sh builds in this directory — after MAFIGATE_PYTHON,
# which still wins, and before a bare python3, which is what is left when neither exists.
#
# Following the scripts here is not a nicety. A developer who runs ./setup.sh and then
# `make test` would otherwise install into one interpreter and test in another, which is
# issue #162 rebuilt out of the two halves this repo just spent a ticket joining.
#
# Both facts are *asked of mafigate_python.sh* rather than restated here. That file is where
# the two scripts read them from, and a third copy of `.venv/bin/python` in this Makefile
# would be exactly the kind of duplicate that file exists to prevent: it would keep working
# after the environment moved, and quietly point make at the wrong python. It is sourced
# with `sh`, not bash, which is why nothing in it uses a bashism.
#
# Made absolute here because test-cov below runs pytest from a copy of this tree in /tmp,
# where a relative path would silently mean *that* copy's environment.
VENV_DIR := $(shell sh -c '. ./mafigate_python.sh && printf "%s" "$$MAFIGATE_VENV"' 2>/dev/null)
VENV_PYTHON := $(shell sh -c '. ./mafigate_python.sh && mafigate_venv_python' 2>/dev/null)
VENV_PYTHON := $(if $(VENV_PYTHON),$(CURDIR)/$(VENV_PYTHON))
PYTHON ?= $(if $(MAFIGATE_PYTHON),$(MAFIGATE_PYTHON),$(if $(VENV_PYTHON),$(VENV_PYTHON),python3))

check-vendor:
	@echo "🔒 Checking vendored pipeline code against bin/..."
	@$(PYTHON) vendor/_sync.py --check

sync-vendor:
	@echo "♻️  Re-copying bin/ into vendor/..."
	@$(PYTHON) vendor/_sync.py

# The parameter contract (issue #31). Sibling to check-vendor and a different kind of
# drift: that guard covers copied *code*, this one copied *values* — nextflow.config's
# params block, restated in config/pipeline_params.py because the packaged app ships no
# nextflow.config to read.
#
# Deliberately NOT a build gate. A moved pipeline default is not a broken build: it is a
# decision about whether the app should follow, and following it moves the parity
# baselines. So it fails in pytest and in CI, where a human is reading, rather than
# stopping a .dmg from being cut. It also needs pytest, which check-vendor deliberately
# does not.
check-params:
	@echo "🔒 Checking the parameter contract against nextflow.config..."
	@$(PYTHON) -m pytest tests/test_param_contract.py -q

# Build commands (standalone installers; see build/BUILD_INSTRUCTIONS.md)
#
# There is no VERSION variable here any more (issue #260). It was `VERSION ?= <a literal>` —
# a fourth copy of a number, a whole major version behind what the app's own APP_VERSION
# said — and, being an override, a way for anyone to build a DMG labelled with anything. The
# version now has one home, config/constants.py, and one reader, build/version.py.
#
# The digits are not quoted here on purpose: a comment recording a stale copy by reprinting
# it leaves a copy in the file, which is what the guard sweeps for.
#
# Read inside the recipes rather than as a `$(shell …)` variable at the top of the file, so
# `make help` and `make clean` do not need a working python3 to expand something they never
# print.
.PHONY: version release-tag release-preflight
version:
	@$(PYTHON) build/version.py

# What `git tag` should be handed for a release. Its own namespace, because the public
# repository's v1.x tags are the Nextflow pipeline's, on an independent cadence (#229).
release-tag:
	@$(PYTHON) build/version.py --tag

# Run this immediately before pressing Publish. It refuses a draft built from a tree the
# repository has moved on from, which is how mafigate-v1.0.0 came within one click of
# shipping a known wrong-variant bug under the fixed version's number (#328). TAG= to
# check a rehearsal tag instead of the one APP_VERSION derives.
release-preflight:
	@$(PYTHON) build/release_preflight.py $(if $(TAG),--tag $(TAG),)

build-mac: check-vendor
	@echo "🍎 Building macOS .dmg installer..."
	@if [ "$$(uname)" = "Darwin" ]; then \
		build/mac/build_dmg.sh; \
	else \
		echo "❌ macOS builds can only be created on macOS systems"; \
	fi

build-win: check-vendor
	@echo "🪟 Building Windows .exe installer..."
	@if command -v python >/dev/null 2>&1 && python -c "import sys; sys.exit(0 if sys.platform == 'win32' else 1)"; then \
		cd build/windows && cmd /c build_installer.bat; \
	else \
		echo "❌ Windows builds must be created on Windows (see build/BUILD_INSTRUCTIONS.md)"; \
	fi

# The third drift guard, and the youngest (issue #162). Its siblings above check copied
# code and copied values; this one checks that the interpreter about to run the app has
# what requirements.txt declares — because a declaration is not an installation, and
# "installed" is a property of one interpreter rather than of the machine.
#
# Presence only, deliberately not versions: every requirement pins one exact version since
# issue #256, but this check runs on bare machines before anything is installed, where
# absence is the answer that matters and ordering releases is a packaging library's job.
# See check_dependencies.py.
.PHONY: check-deps
check-deps:
	@echo "🔍 Checking this interpreter has what MAFigate needs..."
	@$(PYTHON) check_dependencies.py

# Development commands
#
# run depends on check-deps for the reason run_mafigate.sh refuses to launch: the app
# degrades quietly when a chart library is missing, so the check belongs before the boot
# rather than in the reader's chart. Launched through $(PYTHON) -m streamlit so the
# interpreter that was checked is the interpreter that runs.
run: check-deps
	@echo "🌟 Starting MAFigate application..."
	$(PYTHON) -m streamlit run MAFigate.py

# Through $(PYTHON) -m pytest, so the suite runs in the same interpreter as the app. A
# suite run somewhere else is how issue #162 stayed invisible: six chart tests skip when
# plotly is absent, and "absent" was true of the app's interpreter.
test:
	@echo "🧪 Running tests..."
	$(PYTHON) -m pytest tests/ -v

# Coverage lives here rather than in pytest.ini's addopts, so that plain `make test`
# works on a stock install (requirements.txt does not ship pytest-cov, and an
# unrecognised --cov makes pytest exit 4 before collecting anything).
#
# MEASURED WITH bin/ ABSENT (issue #41), which is why this target builds a throwaway tree
# whose repository root holds nothing but streamlit_app/ and runs there.
#
# It is not a nicety. The parity module calls the app's entry point *in process*, so with
# bin/ present those runs count as coverage of the app and with bin/ absent they skip and
# do not. The same command therefore reported two different numbers depending on what was
# checked out beside the app — measured here at 39.58% with bin/ and 39.15% without, which
# pytest-cov's own summary line rounds to 40% and 39% — and a floor that moves with the
# checkout is not a floor: a CI job without the pipeline would fail a threshold a
# developer's machine had just passed. Pinning the measurement to the smaller,
# reproducible half also pins it to the shape the packaged .dmg/.exe actually ships.
#
# It is a real second check as well as a stable measurement. The suite has to be *green*
# in that tree for the floor to be reached at all, which is how three tests in
# test_param_contract.py were found erroring rather than skipping without nextflow.config
# — a state `make test` on a full checkout cannot see, and the state every packaged build
# is in.
#
# The floor is a ratchet against regression, not an aspiration. Measured on this tree:
# 39% (2,812 statements, 1,711 missed), with vendor/, tests/ and build/ omitted via
# .coveragerc. COV_MIN is set a few points below that so an ordinary refactor does not
# trip it, and it is deliberately no longer the 15 that issue #28 specified: 15 was chosen
# against a 17% measurement, and the parity work that closed since has more than doubled
# real coverage, which would leave the floor unable to catch the loss of half the suite.
# Raise it again as coverage improves. (The 80% that used to sit in pytest.ini never ran
# at all — the file's section header was wrong.)
#
# Needs pytest-cov, which requirements.txt still does not declare (issue #19).
#
# The copy drops the virtual environment along with the caches (issue #258). `--cov=.`
# measures the whole tree, so a copied environment would put every installed package's source
# into the denominator — streamlit and pandas included — and the floor below would stop
# describing this app. $(PYTHON) is absolute, so the run still uses the real environment.
#
# Written as $(if $(VENV_DIR),…) rather than as a bare path: $(VENV_DIR) is empty when
# mafigate_python.sh cannot be read, and an empty one inside the quotes would leave
# `rm -rf "$$tree/streamlit_app/"` — the whole copied tree. The argument disappears instead.
COV_MIN ?= 35

# The premise of the whole measurement, checked rather than assumed: inside the copied
# tree the harness must report no pipeline. Asking the harness is the point — an
# `[ ! -d bin ]` on a directory we just created ourselves can never fail, and a guard
# that cannot fire is decoration. This one can: run it in the real checkout and it exits
# 1, which is how it was verified. If the copy ever grows a repository root, or BIN_DIR
# stops being what decides, this stops the run instead of silently measuring the wrong
# tree and calling the number stable.
define PREMISE
import sys
sys.path.insert(0, ".")
from tests.parity.harness import PIPELINE_AVAILABLE
sys.exit(1 if PIPELINE_AVAILABLE else 0)
endef
export PREMISE

test-cov:
	@echo "🧪 Running tests with coverage (floor $(COV_MIN)%, measured with bin/ absent)..."
	@$(PYTHON) -c "import pytest_cov" 2>/dev/null || { \
		echo "❌ pytest-cov is not installed, and pytest exits 4 on --cov without it —"; \
		echo "   before collecting anything, so the failure says nothing about coverage."; \
		echo "   requirements.txt lists it under Development / tests: pip install pytest-cov"; \
		exit 1; \
	}
	@parent=$$(mktemp -d) && tree="$$parent/pipeline-less" && mkdir -p "$$tree" && \
	cp -R "$(CURDIR)" "$$tree/streamlit_app" && \
	rm -rf "$$tree/streamlit_app/htmlcov" "$$tree/streamlit_app/.pytest_cache" \
		$(if $(VENV_DIR),"$$tree/streamlit_app/$(VENV_DIR)") && \
	if ! ( cd "$$tree/streamlit_app" && $(PYTHON) -c "$$PREMISE" ); then \
		echo "❌ the harness still sees a pipeline in the coverage tree"; \
		rm -rf "$$parent"; exit 1; \
	fi; \
	status=0; \
	( cd "$$tree/streamlit_app" && $(PYTHON) -m pytest tests/ -v \
		--cov=. --cov-report=term-missing --cov-report=html \
		--cov-fail-under=$(COV_MIN) ) || status=$$?; \
	if [ -d "$$tree/streamlit_app/htmlcov" ]; then \
		rm -rf "$(CURDIR)/htmlcov"; \
		cp -R "$$tree/streamlit_app/htmlcov" "$(CURDIR)/htmlcov"; \
	fi; \
	rm -rf "$$parent"; \
	exit $$status

test-fast:
	@echo "⚡ Running fast tests only..."
	$(PYTHON) -m pytest tests/ -v -m "not slow"

# The parity sweep over the full 50-sample GERSOM reference (issue #25, issue #27).
#
# Deliberately NOT a pytest test. The reference is a shared clinical drive, present on
# a handful of machines and in no CI job, so a skip-if-absent gate would be green
# everywhere by skipping — the exact failure issue #24 found in test_parity.py, where a
# module-level skipif silently removed the whole parity suite from any checkout without
# bin/. A target you have to type cannot be mistaken for a check that ran.
#
# It is also not a docs-only script any more, which is how the blind spot survived:
# reference.py was written, measured, and then never run again, while the fixture
# baseline went on reporting "every diverging row attributed" — true of the fixtures,
# false of the reference by 3,687 rows, all of them divergence #6.
#
# --strict asserts the two invariants that need no committed numbers, so neither can go
# stale the way a recorded count does:
#
#   1. every diverging row is explained by some probe (`unattributed` == 0). This is the
#      property that was false by 3,687 rows while the fixture baseline read as
#      complete. Run it whenever PROBES or the fixture cells change.
#   2. the four-cell decomposition reproduces bin/filter_variants.py's PASS set exactly
#      (`decomposition_mismatches` == 0), compared per sample as a set of variant keys
#      in a fresh subprocess. Strictly stronger than any PASS total, which two runs can
#      agree on by cancelling errors.
#
# The everywhere-cheap half of invariant 1 lives in
# tests/parity/test_baseline_integrity.py, which reads baseline.json and needs
# neither bin/ nor the reference.
#
# PARITY_REFERENCE_ROOT must be set: the reference lives on a clinical shared drive whose
# path names the institute, the unit and the project, so issue #247 took it out of the
# tree rather than strip this module from the public export. reference.py says so, with
# the variable's name, if it is missing.
parity-reference:
	@echo "🔬 Parity sweep over the full reference (~3 min, needs PARITY_REFERENCE_ROOT)..."
	$(PYTHON) tests/parity/reference.py --strict --report --jobs $${JOBS:-8}

# Issue #32's runtime acceptance criterion: "loading and browsing a MAF in the running
# app still works, including the upload path". Not a pytest test — it boots the whole
# Streamlit app twice through streamlit.testing.v1.AppTest, once per load call site
# ("Open With" and upload), and asserts the page raised nothing, loaded every row,
# filtered, rendered, and got the pipeline's dtypes.
#
# A target rather than a docs-only script, for the reason parity-reference gives above: a
# script that is written, measured once and never run again is how a blind spot survives.
#
# IN tests/ RATHER THAN UNDER docs/, WHICH IS WHERE IT WAS WRITTEN (issue #345). It ran
# from `../docs/wayfinder/issue-32/` until the day someone cloned the public repository and
# found that `.publicignore` strips that whole tree: this target — the one the README calls
# the quickest evidence a change has not broken the app — stopped with a missing-file
# traceback in every public clone, while the page advertising it worked everywhere else.
#
# The move costs nothing and the two reasons it looked expensive were both false. Neither
# installer ships it: both copy lists are allow-lists, deny-by-default, and the .dmg also
# `rm -rf`s tests/ — all of which tests/test_vendor_drift.py's
# `test_neither_installer_ships_the_test_suite_or_fixtures` already held before this file
# moved, so no new exclusion and no new guard were needed. And unlike parity-reference
# above, nothing here is unable to travel: the three fixtures it loads are already in the
# public tree, and the export gate's own six scan rules find nothing in either file. It was
# stripped for its address, not for its content.
#
# What #345 did need was a guard against the *class*, which is
# tests/test_makefile_targets_are_runnable.py: no recipe may name a path the export strips,
# and every script a recipe runs must be present. The two halves fail in different trees —
# the first here, the second in a public clone — because neither tree can see this bug alone.
#
# Not named test_*.py, so pytest does not collect it: it is a script that boots the app six
# times over ~90s, and `make test` is not where that belongs. The name is the interface.
#
# Declared .PHONY here rather than appended to the list at the top of this file. Make
# accepts any number of .PHONY lines, and keeping this one next to its target means a
# branch adding a target does not touch a line every other such branch also touches —
# which is a merge conflict for no reason. This ticket's own branch hit exactly that.
.PHONY: app-load-check
app-load-check:
	@echo "🖥  Loading a MAF through both of the running app's load paths..."
	$(PYTHON) tests/run_app_check.py .

lint:
	@echo "🔍 Running code linting..."
	@if command -v flake8 >/dev/null 2>&1; then \
		flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics; \
		flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics; \
	else \
		echo "⚠️  flake8 not installed, skipping lint check"; \
	fi

format:
	@echo "✨ Formatting code..."
	@if command -v black >/dev/null 2>&1; then \
		black . --line-length 100; \
	else \
		echo "⚠️  black not installed, skipping formatting"; \
	fi
	@if command -v isort >/dev/null 2>&1; then \
		isort . --profile black; \
	else \
		echo "⚠️  isort not installed, skipping import sorting"; \
	fi

# Maintenance commands
clean:
	@echo "🧹 Cleaning temporary files..."
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "htmlcov" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name ".coverage" -delete 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf dist/ .tox/ .coverage.*
	@# The build stamp a local installer build left behind (issue #263). Removed because it
	@# outlives the build that wrote it: run the app from this tree afterwards and About says
	@# it is a .dmg, which is false of what is running. Its absence is the app's source-checkout
	@# answer, so deleting it restores the truth rather than losing information. Nothing in a
	@# build depends on this — `build-mac` and `build-win` write it themselves, and neither
	@# runs `clean`.
	rm -f config/build_stamp.py
	@echo "✅ Cleanup complete"

deps-check:
	@echo "🔍 Checking for dependency updates..."
	@if command -v pip-check >/dev/null 2>&1; then \
		pip-check; \
	else \
		echo "⚠️  pip-check not installed, install with: pip install pip-check"; \
	fi

# Quality assurance
qa: lint test
	@echo "✅ Quality assurance complete"

# Full development cycle
dev: clean format lint test
	@echo "🎉 Development cycle complete"

# Production build preparation
build: clean format lint test
	@echo "🚀 Build preparation complete"
