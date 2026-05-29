# xrd-skill

XRD data analysis skill: parse, background-subtract, peak-find, phase-match.

## Install

```bash
# Claude Code
git clone <repo-url> ~/.claude/skills/xrd-skill

# Universal (.agents)
git clone <repo-url> ~/.agents/skills/xrd-skill

# Or run installer
./install.sh
```

## Usage

```
/xrd-skill Analyze this XRD file
/xrd-skill 分析这个XRD数据，预期物相是Cu(OH)₂
/xrd-skill Find peaks and match phases
```

Triggers on: XRD, xrd, 衍射, 物相, 寻峰, phase matching

## Requires

- Python 3.9+ with numpy
- nature-figure skill (for plotting)
- Notion API key (optional, for result upload)
