#!/bin/bash

for dir in */; do

    species=${dir%/}
    runs_file="${species}/${species}.runs.tsv"
    output_file="${species}/${species}_SI_gene.tsv"

    echo "Processing $species"

    if [ ! -f "$runs_file" ]; then
        echo "  Skipping (no runs file)"
        continue
    fi

    awk -F'\t' '
    NR>1 {

        gsub("gene-","",$3)
        g=$3

        if(!(g in max_obs) || $15>max_obs[g])
            max_obs[g]=$15

        sum_exp[g]+=$17
        count[g]++
        sum_risky[g]+=$16
        sum_entropy[g]+=$11
        gene_len[g]=$6
    }

    END{
        print "gene_id","obs_max_nt_homopolymer",
              "mean_expected_max_nt_homopolymer",
              "SI","total_risky",
              "mean_codon_entropy","gene_length"

        for(g in max_obs){

            mean_exp=sum_exp[g]/count[g]

            if(mean_exp==0)
                SI=0
            else
                SI=max_obs[g]/mean_exp

            mean_entropy=sum_entropy[g]/count[g]

            print g,
                  max_obs[g],
                  mean_exp,
                  SI,
                  sum_risky[g],
                  mean_entropy,
                  gene_len[g]
        }
    }
    ' OFS='\t' "$runs_file" > "$output_file"

done

echo "DONE generating SI_gene for all species."
