cd /tmp/claude-0/-home-claude/d69497fa-df66-5b49-b3d2-73e86eed1d01/scratchpad
while IFS= read -r line; do printf '%s\0' "$line"; done < specs.txt | xargs -0 -P 2 -I{} sh -c 'python3 cmfit.py spec "$1" >> log_specs.txt 2>&1' _ {}
echo ALLDONE >> log_specs.txt
