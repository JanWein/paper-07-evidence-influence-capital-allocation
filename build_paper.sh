#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p build
for pass in 1 2 3; do
  pdflatex -no-shell-escape -halt-on-error -interaction=nonstopmode -output-directory=build manuscript.tex > "build/pass_${pass}.log"
done
echo "Built build/manuscript.pdf"
