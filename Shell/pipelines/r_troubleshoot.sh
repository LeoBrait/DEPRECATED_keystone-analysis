source Shell/settings.sh
mkdir -p "logs"
conda deactivate
Rscript R/pipelines/r_troubleshoot.R > logs/r_libs_report.out