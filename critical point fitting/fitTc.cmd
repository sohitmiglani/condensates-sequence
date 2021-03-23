#! /bin/bash
#
#SBATCH --nodes=1 # node count
#SBATCH --ntasks-per-node=21
#SBATCH --time=3:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=bweiner@princeton.edu

declare -a arr=(
"L8_b1 1.1332188333 0"
"L8_b1 1.1332188333 1"
"L8_b1 1.1332188333 2"
"L8_b2 1.14977861975 0"
"L8_b2 1.14977861975 1"
"L8_b2 1.14977861975 2"
"L8_b4 1.09911374365 0"
"L8_b4 1.09911374365 1"
"L8_b4 1.09911374365 2"
"L12_b1 1.0208149553499999 0"
"L12_b1 1.0208149553499999 1"
"L12_b1 1.0208149553499999 2"
"L12_b2 1.04614831116 0"
"L12_b2 1.04614831116 1"
"L12_b2 1.04614831116 2"
"L12_b3 1.02424688029 0"
"L12_b3 1.02424688029 1"
"L12_b3 1.02424688029 2"
"L12_b6 0.9593788350110001 0"
"L12_b6 0.9593788350110001 1"
"L12_b6 0.9593788350110001 2"
)

module load anaconda3/5.0.1

for i in "${arr[@]}"
do
python fitTc.py $i &
done
wait