#!/bin/bash -e

declare -A overleaf

overleaf[vub-article]=5efb2eb292c0b10001045e23
overleaf[vub-thesis]=5efb2eb392c0b10001045e38
overleaf[bruface-thesis]=5efb2eb492c0b10001045e4e
overleaf[vub-presentation]=5efb2eb692c0b10001045e63

TARGET=$1
TARGET_DIR="$TARGET-overleaf"
OVERLEAF_URL="https://git.overleaf.com/${overleaf[$TARGET]}"

git config --global credential.helper store
echo "\
protocol=https
host=git.overleaf.com
username=$OVERLEAF_USERNAME
password=$OVERLEAF_PASSWORD" | git credential fill | git credential approve

git clone "$OVERLEAF_URL" "$TARGET_DIR"

# Clean up everything
rm -rf "$TARGET_DIR"/*
unzip -d "$TARGET_DIR" "template-$TARGET.zip"
(
    cd "$TARGET_DIR"
    ls -al
    git add .
    git commit -m "Update from Gitlab." || echo "Nothing to do."
    git push
)

rm -rf "$TARGET_DIR"
