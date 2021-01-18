#! /bin/bash
#
# Nom du job
#SBATCH -J sc_revcov
#
# Partition visee
#SBATCH --partition=fast
#
#SBATCH --mem=90G
#SBATCH --cpus=50
# 
# Temps de presence du job
#SBATCH --time=0-5:00:00
#
# Adresse mel de l'utilisateur
#SBATCH --mail-user=elise.amblard@curie.fr
#
# Envoi des mails
#SBATCH --mail-type=abort,end
#
# Nom du fichier de log de Slurm pour le job
#SBATCH -o /home/externe/cnrs/aregimbeau/dev/elise/revcov_log.out

module load r/4.0.2

Rscript 2_rev_cov_hubs.R
