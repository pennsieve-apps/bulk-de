.PHONY: help test fixture run validate clean

SERVICE_NAME ?= bulk-de
PYTHON ?= python3
INPUT ?= ./data/fixture/rejoin

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test      - pytest: the synthetic fixtures through the processor (no Docker)"
	@echo "make fixture   - write the synthetic fixtures to data/fixture/{rejoin,sparc,genesets}/"
	@echo "make run       - build the image and run it in Docker over INPUT=$(INPUT) (dev.env), writing to data/output/"
	@echo "make validate  - app.yml against the manifest schema"
	@echo "make clean     - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture

run: clean
	@if [ "$(INPUT)" = "./data/fixture/rejoin" ]; then $(MAKE) fixture PYTHON=$(PYTHON); fi
	mkdir -p data/output
	INPUT=$(INPUT) docker compose -f docker-compose.yml down --remove-orphans
	INPUT=$(INPUT) docker compose -f docker-compose.yml build
	INPUT=$(INPUT) docker compose -f docker-compose.yml up --exit-code-from processor
	INPUT=$(INPUT) docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

validate:
	$(PYTHON) scripts/validate.py

clean:
	rm -rf data/output
