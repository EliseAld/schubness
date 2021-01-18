#! /bin/bash
#
# Nom du job
#SBATCH -J stab100_k
#
# Partition visee
#SBATCH --partition=long
#
#SBATCH --mem=80G
#SBATCH --cpus=45
# 
# Temps de presence du job
#SBATCH --time=2-23:50:00
#
# Adresse mel de l'utilisateur
#SBATCH --mail-user=elise.amblard@curie.fr
#
# Envoi des mails
#SBATCH --mail-type=abort,end
#
# Nom du fichier de log de Slurm pour le job
#SBATCH -o /home/externe/curie/eamblard/dev/TI/stab100_ktuned.out

module load r/4.0.2

Rscript stab_100_SLURM_IHopeThisOneIsCorrect_ktuned.R
