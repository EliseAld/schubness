#! /bin/bash
#
# Nom du job
#SBATCH -J sc_2
#
# Partition visee
#SBATCH --partition=fast
#
#SBATCH --mem=50G
#SBATCH --cpus=30
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
#SBATCH -o /home/externe/cnrs/aregimbeau/dev/elise/2_log.out

module load r/4.0.2

Rscript 2_hubness_param.R
