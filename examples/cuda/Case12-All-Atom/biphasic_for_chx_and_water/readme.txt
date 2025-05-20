Description: 
A biphasic system consisting of cyclohexane and water.

MD in gromacs 2024.4:
gmx grompp -f test_NPT.mdp -c eq_NPT.gro -p topol.top -o test_NPT.tpr -maxwarn 1
gmx mdrun -deffnm test_NPT -gpu_id 0
gmx energy -f test_NPT.edr

Results:
Volume: 194.091614      Potential:-49548.703125

Computational efficiency:
GeForce RTX 4080 SUPER: 356.423 ns/day 

MD in pygamd 1.4.7:
python gro_to_xml.py --gro=eq_NPT.gro --top=topol.top 
python gala_NPT.py --gpu=0

Results:
volume: 192.6091766     total_potential: -47726.875

Computational efficiency:
GeForce RTX 4080 SUPER: 132.624 ns/day

***NOTE***
1. All parameters should be located at ffnonbonded.itp and ffbonded.itp.
2. DO NOT use any ifdef or ifndef statement in any itp files.
3. All include statements should be moved to topol.top.
4. Creating a different .itp file for each different molecule contains bonding information.
5. For structures other than the example, please use fully equilibrated structures.