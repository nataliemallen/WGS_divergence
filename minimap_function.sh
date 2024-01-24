run_minimap() {
    if [ "$#" -ne 3 ]; then
        echo "Usage: run_minimap <file1> <file2> <name>"
        return 1
    fi

    file1="$1"
    file2="$2"
    name="$3"

    minimap_command() {
        minimap2 -cx asm"$1" -t8 --cs "$2" "$3" > "${4}_${1}.paf"
        sort -k6,6 -k8,8n "${4}_${1}.paf" > "${4}_${1}.srt.paf"
        paftools.js call "${4}_${1}.paf" > "${4}_${1}.var.txt"
    }

    minimap_command 5 "$file1" "$file2" "$name"
    minimap_command 10 "$file1" "$file2" "$name"
    minimap_command 20 "$file1" "$file2" "$name"
}

# Example usage:
# run_minimap "chicken" "turkey" chicken_turkey
