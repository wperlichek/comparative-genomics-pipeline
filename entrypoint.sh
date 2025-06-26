#!/bin/bash
set -e

# Run the pipeline if no arguments are given, else run the provided command
if [ $# -eq 0 ]; then
    comparative-genomics-pipeline || true
    echo "Pipeline finished. Container will stay alive. Press Ctrl+C to exit."
    tail -f /dev/null
else
    exec "$@"
fi
