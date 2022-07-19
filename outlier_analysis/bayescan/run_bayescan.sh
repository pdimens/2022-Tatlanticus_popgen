#! /usr/bin/env bash

BayeScan2.1 $1 \
	-o bft_bscan \
	-n 15000 \
	-pilot 15000 \
	-burn 100000 \
	-pr_odds 100 \
	-threads $2 \
