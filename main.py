from vina import Vina

v = Vina(sf_name='vina')

# inputs
receptor_list = ['a.pdbqt','b.pdbqt','c.pdbqt','d.pdbqt','e.pdbqt']
ligand_list = ['a_ligand.pdbqt','b_ligand.pdbqt','c_ligand.pdbqt','d_ligand.pdbqt','e_ligand.pdbqt']
filename_list = ['1AX8','1B9G','1D2S','3N9Y','3V6O']
center_cords_list = [[55.40,-33.70,5.30],[1.2,0.7,1.8],[26.8,16.9,36.1],[32.4,-1.1,-8.2],[55.6,51.9,-2.6]]
box_size_cords_list = [[50,50,50],[30,30,30],[40,40,40],[120,80,100],[90,100,150]]
exhaust=200
poses=10

for i in range(len(receptor_list)):

    v.set_receptor(receptor_list[i])
    v.set_ligand_from_file(ligand_list[i])
    filename = filename_list[i]
    print(filename)

    v.compute_vina_maps(center=center_cords_list[i], box_size=box_size_cords_list[i])

    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(f'200/{filename}_ligand_minimized.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=exhaust, n_poses=poses)
    v.write_poses(f'{exhaustiveness}/{filename}_ligand_vina_out.pdbqt', n_poses=30, overwrite=True)
