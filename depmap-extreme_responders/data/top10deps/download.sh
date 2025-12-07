cat ../../results/n1s_tbl.tsv |awk '{ print $2 }' |grep "^ACH" |xargs -I@ -P4 bash -c 'wget "https://depmap.org/portal/cell_line/prefdep/crispr/@" --quiet -O @.json'
