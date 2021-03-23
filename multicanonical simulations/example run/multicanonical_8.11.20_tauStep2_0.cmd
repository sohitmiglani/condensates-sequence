#! /bin/bash
#
#SBATCH --nodes=1 # node count
#SBATCH --ntasks-per-node=18
#SBATCH --time=72:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=bweiner@princeton.edu

declare -a arr=(
"./multicanonical multicanonical_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b12.txt annealing_schedule_e0.7716181316764409.txt 100 300000000 4000 0.05 preweighting_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409.txt 6920"
"./multicanonical multicanonical_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b12.txt annealing_schedule_e0.7716181316764409.txt 100 300000000 4000 0.05 preweighting_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409.txt 123"
"./multicanonical multicanonical_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b12.txt annealing_schedule_e0.7716181316764409.txt 100 300000000 4000 0.05 preweighting_L24_b12_beta0_0.7692846091008065_beta1_0.7716181316764409.txt 6611"
"./multicanonical multicanonical_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b1.txt annealing_schedule_e0.8967850458479947.txt 100 300000000 4000 0.05 preweighting_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947.txt 8915"
"./multicanonical multicanonical_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b1.txt annealing_schedule_e0.8967850458479947.txt 100 300000000 4000 0.05 preweighting_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947.txt 8465"
"./multicanonical multicanonical_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b1.txt annealing_schedule_e0.8967850458479947.txt 100 300000000 4000 0.05 preweighting_L24_b1_beta0_0.8940729942980511_beta1_0.8967850458479947.txt 9722"
"./multicanonical multicanonical_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b2.txt annealing_schedule_e0.9286263380879677.txt 100 300000000 4000 0.05 preweighting_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677.txt 2351"
"./multicanonical multicanonical_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b2.txt annealing_schedule_e0.9286263380879677.txt 100 300000000 4000 0.05 preweighting_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677.txt 7094"
"./multicanonical multicanonical_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b2.txt annealing_schedule_e0.9286263380879677.txt 100 300000000 4000 0.05 preweighting_L24_b2_beta0_0.9258179923074595_beta1_0.9286263380879677.txt 5002"
"./multicanonical multicanonical_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b3.txt annealing_schedule_e0.912976328446916.txt 100 300000000 4000 0.05 preweighting_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916.txt 2222"
"./multicanonical multicanonical_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b3.txt annealing_schedule_e0.912976328446916.txt 100 300000000 4000 0.05 preweighting_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916.txt 5463"
"./multicanonical multicanonical_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b3.txt annealing_schedule_e0.912976328446916.txt 100 300000000 4000 0.05 preweighting_L24_b3_beta0_0.9102153113245967_beta1_0.912976328446916.txt 2907"
"./multicanonical multicanonical_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b4.txt annealing_schedule_e0.8902460349221437.txt 100 300000000 4000 0.05 preweighting_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437.txt 8667"
"./multicanonical multicanonical_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b4.txt annealing_schedule_e0.8902460349221437.txt 100 300000000 4000 0.05 preweighting_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437.txt 2907"
"./multicanonical multicanonical_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b4.txt annealing_schedule_e0.8902460349221437.txt 100 300000000 4000 0.05 preweighting_L24_b4_beta0_0.8875537586068549_beta1_0.8902460349221437.txt 7858"
"./multicanonical multicanonical_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395_j0.05_rep0 lattice_fcc.txt polySpecs_L24_b6.txt annealing_schedule_e0.8484494215672395.txt 100 300000000 4000 0.05 preweighting_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395.txt 3678"
"./multicanonical multicanonical_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395_j0.05_rep1 lattice_fcc.txt polySpecs_L24_b6.txt annealing_schedule_e0.8484494215672395.txt 100 300000000 4000 0.05 preweighting_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395.txt 2259"
"./multicanonical multicanonical_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395_j0.05_rep2 lattice_fcc.txt polySpecs_L24_b6.txt annealing_schedule_e0.8484494215672395.txt 100 300000000 4000 0.05 preweighting_L24_b6_beta0_0.8458835463004031_beta1_0.8484494215672395.txt 8441"
)

for i in "${arr[@]}"
do
$i &
done
wait
