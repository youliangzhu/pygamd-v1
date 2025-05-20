Description: 
Ubiquitin in water.

MD in gromacs 2024.4:
gmx grompp -f test_NPT.mdp -c eq_NPT.gro -p topol.top -o test_NPT.tpr -maxwarn 1
gmx mdrun -deffnm test_NPT -gpu_id 3
gmx energy -f test_NPT.edr

Results:
Volume: 221.554596      Potential: -283311.031250

Computational efficiency:
GeForce RTX 4080 SUPER: 580.643 ns/day 

MD in pygamd 1.4.7:
python gro_to_xml.py --gro=eq_NPT.gro --top=topol.top 
python gala_NPT.py --gpu=3

Results:
volume: 220.7103882     total_potential: -279253.625

Computational efficiency:
GeForce RTX 4080 SUPER: 117.043 ns/day

***NOTE***
1. All parameters should be located at ffnonbonded.itp and ffbonded.itp.
2. DO NOT use any ifdef or ifndef statement in any itp files.
3. All include statements should be moved to topol.top.
4. Creating a different .itp file for each different molecule contains bonding information.
5. For structures other than the example, please use fully equilibrated structures.