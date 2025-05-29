#!/bin/bash
echo "Auto commit script is running..."
script_path=$(realpath "$0")
script_dir=$(dirname "$script_path")

cd "$script_dir"

git add .
git commit -m "Auto commit on $(date)"
git push origin main
