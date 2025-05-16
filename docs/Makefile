# docs/Makefile — Minimal Sphinx build with autosummary

# You can override these from the environment or command line
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
AUTOSUMM_DIR  = $(SOURCEDIR)/_autosummary

.PHONY: help clean-autosummary html clean

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

# Remove old autosummary stubs so we always regenerate
clean-autosummary:
	@echo "Cleaning old autosummary stubs…"
	@rm -rf "$(AUTOSUMM_DIR)"

# Build HTML — first clean autosummary, then let Sphinx generate everything
html: clean-autosummary
	@echo "Building HTML with autosummary stubs…"
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)"

# Generic passthrough for other Sphinx targets (e.g. latexpdf, linkcheck)
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)"
