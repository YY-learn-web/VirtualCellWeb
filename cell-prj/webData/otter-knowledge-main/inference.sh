#!/usr/bin/env bash

CUDA_VISIBLE_DEVICES=0 python inference.py --input_path '/Users/yd/codebuddy/cell-prj/webData/otter-knowledge-main/DB_drugsRenew.csv' \
--sequence_column 'smiles' --input_type Drug --relation_name smiles \
--model_path '/Users/yd/codebuddy/cell-prj/webData/otter-knowledge-main/models/' \
--output_path '/Users/yd/codebuddy/cell-prj/webData/otter-knowledge-main/DB_drugsRenew2vec.json'

