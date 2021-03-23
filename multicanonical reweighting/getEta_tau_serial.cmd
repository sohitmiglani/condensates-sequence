#! /bin/bash
#
#SBATCH --nodes=1 # node count
#SBATCH --ntasks-per-node=6
#SBATCH --time=00:05:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=bweiner@princeton.edu

declare -a arr=(
"L24_b1 0.8869204103436666 0 0.8913772968278058 0.8940729942980511 -0.011"
"L24_b2 0.918411448369 0 0.9230265812753768 0.9258179923074595 -0.011"
"L24_b3 0.9029335888339999 0 0.9074709435517586 0.9102153113245967 -0.011"
"L24_b4 0.8804533285380001 0 0.8848777171236183 0.8875537586068549 -0.011"
"L24_b6 0.8391164779299999 0 0.8433331436482411 0.8458835463004031 -0.011"
"L24_b12 0.7631303322280001 0 0.7669651580180904 0.7692846091008065 -0.011"
)

module load anaconda3/5.0.1

for i in "${arr[@]}"
do
python multicanonical_getEta_tau.py $i &
done
wait