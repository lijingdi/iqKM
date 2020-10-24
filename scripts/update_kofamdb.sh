#!/bin/bash

script_dir=`dirname $0`

## step1 Download KoFamKOALA profile files ##
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -P "$script_dir"/../db/ && tar -xvzf "$script_dir"/../db/profiles.tar.gz &&
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz -P "$script_dir"/../db/ && gunzip "$script_dir"/../db/ko_list.gz &&

## step2 Add GA threshold to each profile ##
python3 "$script_dir"/update_kofam.py "$script_dir"/../db/profiles "$script_dir"/../db/ko_list "$script_dir"/../db/ &&

## step3 Concatenate all the adjusted hmm files to Kofam.db ##
cat "$script_dir"/../db/profiles/*.hmm > "$script_dir"/../db/kofam.hmm &&

## step4 Press hmm.db for hmmsearch ##
hmmpress "$script_dir"/../db/kofam.hmm


### if Kofam db is updated, the minimum distance within KM need to be updated too
