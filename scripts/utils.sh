#!/bin/sh

function doc2odt {
    if [ -z $1 ]; then
        echo "Usage: doc2odt folder"
    else
        owd=`pwd`
        cd $1
        libreoffice --headless --convert-to odt *.doc*
        cd $owd
   fi
}   

function Rdoc {
    if [ -f $1 ]; then 
        Rscript -e "cat(tools::Rd2txt('$1'),file=stdout());"
    elif [ -z $2 ]; then
        Rscript -e "help($1)"
    else
        Rscript -e "library($1); help($2);"    
    fi
}

function rgb2hex {
    if [ -z $3 ]; then
        echo "Usage: rgb2hex R G B"
    else 
        hex=`printf "#%02x%02x%02x\n" $1 $2 $3`
        printf "$1 $2 $3 == $hex "
        printf "  \033[48;2;%d;%d;%dm           \033[0m " $1 $2 $3
        printf "  \033[38;2;%d;%d;%dmHello World\033[0m\n" $1 $2 $3
    fi
}

function colors-16 { 
    printf "\n   "
    for c in {0..15};
    do
        printf "\033[48;5;%sm %3d \033[0m" "$c" $((c+1)); #)
        if (( c == 7 )); then
            printf "\n   ";
        fi;
    done
    printf "\n\n";
    
    printf "\n   "
    for c in {0..15};
    do
        printf "\033[38;5;%sm %3d \033[0m" "$c" "$c";
        if (( c == 7 )); then
            printf "\n   ";
        fi;
    done
    printf "\n\n";
}

function curl-cheat {
    curl cheat.sh/$1
}
function colors-show {
    for c in {0..255}; do     
        printf "\033[48;5;%sm%3d\033[0m " "$c" "$c";     
        if (( c == 15 )) || (( c > 15 )) && (( (c-15) % 16 == 0 )); then
            printf "\n";
        fi; 
    done
}

function curl-weather {
    if [ -z $1 ]; then
        curl wttr.in/Potsdam
    else
        curl wttr.in/$1
    fi
}

function substring () {
    if [[ -z $2 ]]; then
        echo "Usage: substring string start [end]"
    elif [[ -z $3 ]]; then
        echo ${1:$2}
    else
        echo ${1:$2:$3}
    fi
}

function csvcols (){
    local search_cols="$1"   # Expecting space-separated string: "alias barcode sample"
    local input_file="$2"
    local format="$3"        # "csv" or "tab"
    local delim

    # Set delimiter based on format
    if [[ "$format" == "tab" ]]; then
        delim=$'\t'
    else
        delim=','
    fi

    # Run awk:
    # 1. Read only the first line (NR==1).
    # 2. Map every column name to its index number.
    # 3. Look up the requested columns in that map.
    # 4. Print indices in order (0 if not found).
    awk -F"$delim" -v cols="$search_cols" '
    NR==1 {
        # Build map: header_name -> column_index
        for(i=1; i<=NF; i++) {
            map[$i] = i
        }
        
        # Split the requested columns string into an array
        n = split(cols, targets, " ")
        
        # Iterate through requested targets and print their index
        for(j=1; j<=n; j++) {
            val = (targets[j] in map) ? map[targets[j]] : 0
            printf "%s ", val
        }
        print "" # Newline
        exit     # Stop reading file immediately after header
    }' "$input_file"
}

# Usage: echo "line" | csvextract "2 5 1" [delimiter]
# Expects input to be a string of indices corresponding to targets cols of a csv file
function csvextract() {
    local indices="$1"
    local delim="${2:-,}"  # Default to comma if not provided

    awk -F"$delim" -v cols="$indices" '
    BEGIN {
        # Split the input string "2 5 1" into an array
        n = split(cols, inds, " ")
    }
    {
        # Loop through the desired indices for every row
        for(i=1; i<=n; i++) {
            col_idx = inds[i]

            # If index is 0 (missing), use empty string; otherwise use column value
            val = (col_idx == 0) ? "" : $(col_idx)

            # Print value + comma (unless it is the last item)
            printf "%s%s", val, (i==n ? "" : ",")
        }
        # Print newline at the end of the row
        print ""
    }'
}


# Reheader bam files
# Usage: echo "line" | csvextract bam_file ["index"]
function bamrehead(){
	local bam="$1"
	local index="$2"
	local temp_bam="${bam}.rehead.tmp.bam"

	samtools reheader <( \
		samtools view -H "$bam" | \
		sed  's/SN:\([0-9XY]\)/SN:chr\1/g' | \
		sed  's/SN:MT/SN:chrM/g' \
	) "$bam" > "$temp_bam"

	# Check if it worked before overwritting
	if [[ $? -eq 0 ]]; then
		mv "$temp_bam" "$bam"
		
		echo "  ✓ Fixed header for: $(basename "$bam")"
		if [[ $index == "index" ]]; then
			samtools index $bam
			echo "  ✓ Complete indexing: ${bam}"
		fi
	else
		echo "  X Error reheadering: $(basename "$bam")"
			rm -f "$temp_bam"
		return 1
	fi
}
