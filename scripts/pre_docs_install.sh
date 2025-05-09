#!/usr/bin/env bash
sed 's/"defs":/"$defs":/; s|#\/defs|#\/$defs|g' nextflow_schema.json > tmp && mv tmp nextflow_schema.json
nf-core pipelines schema docs --format markdown --output docs/params_doc.md