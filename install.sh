#!/usr/bin/env bash
# xrd-skill installer — Cross-platform auto-detection
set -euo pipefail

SKILL_NAME="xrd-skill"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

detect_platform() {
    if [ -d "$HOME/.claude/skills" ]; then echo "claude:$HOME/.claude/skills"
    elif [ -d "$HOME/.agents/skills" ]; then echo "universal:$HOME/.agents/skills"
    elif [ -d ".cursor/rules" ]; then echo "cursor:.cursor/rules"
    elif [ -d ".github/skills" ]; then echo "copilot:.github/skills"
    elif [ -d "$HOME/.gemini/skills" ]; then echo "gemini:$HOME/.gemini/skills"
    else echo "unknown:$HOME/.claude/skills"; fi
}

install_skill() {
    local target_dir="$1"
    mkdir -p "$target_dir"
    cp -R "$SCRIPT_DIR" "$target_dir/$SKILL_NAME"
    echo "Installed $SKILL_NAME to $target_dir/$SKILL_NAME"
}

PLATFORM=$(detect_platform)
IFS=':' read -r name path <<< "$PLATFORM"
install_skill "$path"
echo "Done. Use: /xrd-skill <task>"
