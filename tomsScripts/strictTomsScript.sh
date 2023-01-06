#!/bin/bash

if (( $# != 3 )); then
        echo "Usage: $0 fileToRead targetDoc filtFastq" >&2
fi

fileToRead=$1
targetDoc=$2
filtFastq=$3

discard=0
write="noWrite"
k=0

target1=GCTGGAGGCATAAACCCCAT
target2=ATGGGGTTTATGCCTCCAGC

lSide1=AGGCATG
lSide2=CATGCCT

rSide1=GGGTGC
rSide2=GCACCC

# Levenshtein distance to calculate string distance
# Taken from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Bash
function levenshtein {
    if (( $# != 2 )); then
        echo "Usage: $0 word1 word2" >&2
    elif (( ${#1} < ${#2} )); then
        levenshtein "$2" "$1"
    else
        local str1len=${#1}
        local str2len=${#2}
        local d

        for (( i = 0; i <= (str1len+1)*(str2len+1); i++ )); do
            d[i]=0
        done

        for (( i = 0; i <= str1len; i++ )); do
            d[i+0*str1len]=$i
        done

        for (( j = 0; j <= str2len; j++ )); do
            d[0+j*(str1len+1)]=$j
        done

        for (( j = 1; j <= str2len; j++ )); do
            for (( i = 1; i <= str1len; i++ )); do
                [ "${1:i-1:1}" = "${2:j-1:1}" ] && local cost=0 || local cost=1
                del=$(( d[(i-1)+str1len*j]+1 ))
                ins=$(( d[i+str1len*(j-1)]+1 ))
                alt=$(( d[(i-1)+str1len*(j-1)]+cost ))
                d[i+str1len*j]=$( echo -e "$del\n$ins\n$alt" | sort -n | head -1 )
            done
        done
        echo ${d[str1len+str1len*(str2len)]}
    fi
}

# Read file line by line
while read -r line; do
	# If line is the one containing genomic reads
    if [ $k -eq 0 ]; then
        header=${line}
	elif [ $k -eq 1 ]; then
        seq=${line}
        if [[ ${line} == *${lSide1}* && ${line} == *${rSide1}* ]]; then
            trimLine=${seq#*${lSide1}}
            trimLine=${trimLine%${rSide1}*}
            levenshtein $trimLine $target1 >> $targetDoc
            write="write"
        elif [[ ${line} == *${lSide2}* && ${line} == *${rSide2}* ]]; then
            trimLine=${seq%${lSide2}*}
            trimLine=${trimLine#*${rSide2}}
            levenshtein $trimLine $target2 >> $targetDoc
            write="write"
        else
            ((discard++))
        fi
    elif [[ $k -eq 3 && $write == "noWrite" ]]; then
        k=-1
    elif [[ $k -eq 3 && $write == "write" ]]; then
        cat >>${filtFastq} <<EOF
${header}
${seq}
+
${line}
EOF
        write="noWrite"
	    k=-1
	fi
	((k++))
done < ${fileToRead}

if [[ $write == "write" ]]; then
        cat >>${filtFastq} <<EOF
${header}
${seq}
+
${line}
EOF
fi

echo "$isSubstring were perfect matches (perfect WT sequences)"
echo "${discard} were found in the non ARID1A amplified region, and where therefore discarded"

# # Read file line by line
# cat ${fileToRead} | while read line; do
# 	# If line is the one containing genomic reads
# 	if [ $i -eq 4 ]; then
#         if [[ ${lSide1} != *${line}*  && ${lSide2} != *${line}* && ${rSide1} != *${line}* && ${rSide2} != *${line}* ]]; then
#             if [[ $line == *$target1* ]]; then
#                 levenshtein $line $target1 >> $perfMatch
#                 ((isSubstring++))
#             elif [[ $line == *$target2* ]]; then
#                 levenshtein $line $target2 >> $perfMatch
#                 ((isSubstring++))
#             else
#                 poss1=$(levenshtein $line $target1)
#                 poss2=$(levenshtein $line $target2)
#                 if [ $poss1 -gt $poss2 ]; then
#                     levenshtein $line $target2 >> $targetDoc
#                 else
#                     levenshtein $line $target1 >> $targetDoc
#                 fi			
#             fi
#         else
#             echo Discarded!
#         fi
# 	i=0
# 	fi
# 	((i++))
# done

# echo "$isSubstring were perfect matches (perfect WT sequences)"
# echo "${discard} were found in the non ARID1A amplified region, and where therefore discarded"