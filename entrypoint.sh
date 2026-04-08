#!/bin/bash
set -e

if [ "$#" -gt 0 ]; then
    exec "$@"
fi

echo "CanIonSynth container is ready."
echo "Run 'julia --project=. tests/run_tests.jl' to execute the test suite."
echo "Run generation scripts from /workspace using the checked-in environments."

exec /bin/bash
