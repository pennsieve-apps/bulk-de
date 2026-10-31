#!/usr/bin/env python3
"""app.yml against the Pennsieve app manifest schema (SCHEMA, default: the sibling workflow-bulk-rnaseq repo)."""
import json
import os
import sys
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parent.parent
schema = json.loads(Path(os.environ.get("SCHEMA", ROOT.parent / "workflow-bulk-rnaseq/docs/app-manifest.v1.schema.json")).read_text())
errors = list(jsonschema.validators.validator_for(schema)(schema).iter_errors(yaml.safe_load((ROOT / "app.yml").read_text())))
for e in errors:
    print(f"app.yml: {e.message} at {list(e.path)}")
print(f"app.yml: {'valid' if not errors else f'{len(errors)} error(s)'}")
sys.exit(1 if errors else 0)
