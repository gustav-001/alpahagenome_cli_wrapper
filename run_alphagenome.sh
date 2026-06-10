#!/usr/bin/env bash
set -e

if ! command -v uv &>/dev/null; then
    echo "uv not found! installing now..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
fi

uv run "$(dirname "$0")/alphagenome_pipeline.py" "$@"
