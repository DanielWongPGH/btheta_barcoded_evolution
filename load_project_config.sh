#!/usr/bin/env bash

# Source this file from Bash: source ./load_project_config.sh
_btheta_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PYTHON:-}" ]]; then
    _btheta_config_python="${PYTHON}"
elif [[ -x "${_btheta_config_dir}/.venv/bin/python" ]]; then
    _btheta_config_python="${_btheta_config_dir}/.venv/bin/python"
else
    _btheta_config_python="python3"
fi
eval "$("${_btheta_config_python}" "${_btheta_config_dir}/project_config.py" --shell)"
unset _btheta_config_python
unset _btheta_config_dir
