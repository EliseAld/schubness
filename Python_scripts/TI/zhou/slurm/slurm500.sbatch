#! /bin/bash
#
# Nom du job
#SBATCH -J stab500
#
# Partition visee
#SBATCH --partition=long
#
#SBATCH --mem=50G
#SBATCH --cpus=50
# 
# Temps de presence du job
#SBATCH --time=3-6:00:00
#
# Adresse mel de l'utilisateur
#SBATCH --mail-user=elise.amblard@curie.fr
#
# Envoi des mails
#SBATCH --mail-type=abort,end
#
# Nom du fichier de log de Slurm pour le job
#SBATCH -o /home/externe/curie/eamblard/dev/TI/stab500.out

module load r/4.0.2

Rscript stab_500_SLURM_IHopeThisOneIsCorrect.R
