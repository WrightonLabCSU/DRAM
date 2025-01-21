#!/usr/bin/env bash

#set -o errexit  # exit  on error
#set -o pipefail  # enable pipe fail to prevent things like `error here | true` always succeeding
#[[ "${DEBUG}" == 'true' ]] && set -o xtrace  # print every line in debug
set -o xtrace


kegg_pep_root_dir=$1
output_dir="${2:-./DRAM_output}"

mkdir -p "${output_dir}"

cat "${kegg_pep_root_dir}"/*/*pep > "${output_dir}"/kegg.pep