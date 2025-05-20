'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its
	derivatives are free from error. In no event shall the copyright
	holder or contributors be liable for any indirect, incidental,
	special, exemplary, or consequential loss or damage that results
	from its use. We also have no responsibility for providing the
	service of functional extension of this program to general users.
USER OBLIGATION
	If any results obtained with PYGAMD are published in the scientific
	literature, the users have an obligation to distribute this program
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu,
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''
# read a snapshot file with XML format
import xml.etree.ElementTree as ET

class read_xml:
    def __init__(self, filename):
        self.num_particles = 0
        self.timestep = 0
        self.dimension = 3
        self.xml_read = False
        self.CG_xml = False
        self.invariant_data = False
        self.variant_data = False
        self.read_indicator = {'box': False, 'position': False, 'type': False,
                               'image': False, 'mass': False, 'velocity': False,
                               'charge': False, 'body': False, 'diameter': False,
                               'rotangle': False, 'force': False, 'virial': False,
                               'molecule': False, 'init': False, 'cris': False,
                               'orientation': False, 'quaternion': False,
                               'rotation': False, 'inert': False, 'asphere': False,
                               'patch': False, 'bond': False, 'angle': False,
                               'dihedral': False, 'vsite': False}
        self.init_data()
        self.read_file(filename)

    def read_file(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()

        self.reset_params()
        for element in root.iter():
            if element.tag == 'galamost_xml' and element.attrib.get('version') == '1.3':
                print('info : read xml file with version 1.3')
                self.xml_read = True
                continue
            if element.tag == 'polymer_xml' and element.attrib.get('version') == '1.3':
                print('info : read polymer_xml file with version 1.3')
                self.xml_read = True
                self.CG_xml = True
                continue

            if self.xml_read:
                if element.tag == 'configuration':
                    self.reset_params()
                    self.num_particles = element.attrib.get('natoms')
                    if self.num_particles:
                        self.num_particles_read = True
                        self.num_particles = element.attrib.get('natoms')
                        print("info : number of particles", self.num_particles)
                    else:
                        self.num_particles = len(self.position)
                    continue

                if element.tag == 'configuration':
                    self.reset_params()
                    self.timestep_read = True
                    if self.num_particles_read:
                        self.timestep = element.attrib.get('time_step')
                        print("info : timestep", self.timestep)
                    continue

                if element.tag == 'configuration':
                    self.reset_params()
                    self.dimension_read = True
                    if self.dimension_read:
                        self.dimension = element.attrib.get('dimensions')
                        print("info : dimension", self.dimension)
                    continue

                if element.tag == 'box':
                    self.reset_params()
                    self.box_read = True
                    if self.invariant_data:
                        self.read_indicator['box'] = True
                    if self.box_read:
                        self.box.append(float(element.attrib.get('lx')))
                        self.box.append(float(element.attrib.get('ly')))
                        self.box.append(float(element.attrib.get('lz')))
                        print('info : box ', float(element.attrib.get('lx')), float(element.attrib.get('ly')), float(element.attrib.get('lz')))
                    continue

                if element.tag == 'position':
                    self.reset_params()
                    self.position_read = True
                    if self.invariant_data:
                        self.read_indicator['position'] = True
                        continue
                    if self.position_read:
                        position_content1 = "".join(element.text)
                        position_content2 = position_content1.strip().split()
                        start = 0
                        end = len(position_content2)
                        step = 3
                        for i in range(start, end, step):
                            x = i
                            self.position.append(position_content2[x:x + step])
                    continue

                if element.tag == 'type':
                    self.reset_params()
                    self.type_read = True
                    if self.invariant_data:
                        self.read_indicator['type'] = True
                        continue
                    if self.type_read:
                        type_content1 = "".join(element.text)
                        type_content2 = type_content1.strip().split()
                        for i in type_content2:
                            self.type.append(i)
                    continue

                if element.tag == 'image':
                    self.reset_params()
                    self.image_read = True
                    if self.invariant_data:
                        self.read_indicator['image'] = True
                        continue
                    if self.image_read:
                        image_content1 = "".join(element.text)
                        image_content2 = image_content1.strip().split()
                        start = 0
                        end = len(image_content2)
                        step = 3
                        for i in range(start, end, step):
                            x = i
                            self.image.append(image_content2[x:x + step])
                    continue

                if element.tag == 'mass':
                    self.reset_params()
                    self.mass_read = True
                    if self.invariant_data:
                        self.read_indicator['mass'] = True
                        continue
                    if self.mass_read:
                        mass_content1 = "".join(element.text)
                        mass_content2 = mass_content1.strip().split()
                        for i in mass_content2:
                            self.mass.append(float(i))
                    continue

                if element.tag == "velocity":
                    self.reset_params()
                    self.velocity_read = True
                    if self.invariant_data:
                        self.read_indicator['velocity'] = True
                        continue
                    if self.velocity_read:
                        velocity_content1 = "".join(element.text)
                        velocity_content2 = velocity_content1.strip().split()
                        start = 0
                        end = len(velocity_content2)
                        step = 3
                        for i in range(start, end, step):
                            x = i
                            output_list = [int(j) if j.isdigit() else j for j in velocity_content2[x:x + step]]
                            self.velocity.append(output_list)
                    continue

                if element.tag == "charge":
                    self.reset_params()
                    self.charge_read = True
                    if self.invariant_data:
                        self.read_indicator['charge'] = True
                    if self.charge_read:
                        charge_content1 = "".join(element.text)
                        charge_content2 = charge_content1.strip().split()
                        for i in charge_content2:
                            self.charge.append(i)
                    continue

                if element.tag == "body":#-----------------------------------------
                    self.reset_params()
                    self.body_read = True
                    if self.invariant_data:
                        self.read_indicator['body'] = True
                    continue

                if element.tag == "diameter":#-----------------------------------------
                    self.reset_params()
                    self.diameter_read = True
                    if self.invariant_data:
                        self.read_indicator['diameter'] = True
                    continue

                if element.tag == "rotangle":#-----------------------------------------
                    self.reset_params()
                    self.rotangle_read = True
                    if self.invariant_data:
                        self.read_indicator['rotangle'] = True
                    continue

                if element.tag == "force":#-----------------------------------------
                    self.reset_params()
                    self.force_read = True
                    if self.invariant_data:
                        self.read_indicator['force'] = True
                    continue

                if element.tag == "virial":#-----------------------------------------
                    self.reset_params()
                    self.virial_read = True
                    if self.invariant_data:
                        self.read_indicator['virial'] = True
                    continue

                if element.tag == "molecule":#-----------------------------------------
                    self.reset_params()
                    self.molecule_read = True
                    if self.invariant_data:
                        self.read_indicator['molecule'] = True
                    continue

                if element.tag == "init":#-----------------------------------------
                    self.reset_params()
                    self.init_read = True
                    if self.invariant_data:
                        self.read_indicator['init'] = True
                    continue

                if element.tag == "cris":#-----------------------------------------
                    self.reset_params()
                    self.cris_read = True
                    if self.invariant_data:
                        self.read_indicator['cris'] = True
                    continue

                if element.tag == "orientation":#-----------------------------------------
                    self.reset_params()
                    self.orientation_read = True
                    if self.invariant_data:
                        self.read_indicator['orientation'] = True
                    continue

                if element.tag == "quaternion":#-----------------------------------------
                    self.reset_params()
                    self.quaternion_read = True
                    if self.invariant_data:
                        self.read_indicator['quaternion'] = True
                    continue

                if element.tag == "rotation":#-----------------------------------------
                    self.reset_params()
                    self.rotation_read = True
                    if self.invariant_data:
                        self.read_indicator['rotation'] = True
                    continue

                if element.tag == "inert":#-----------------------------------------
                    self.reset_params()
                    self.inert_read = True
                    if self.invariant_data:
                        self.read_indicator['inert'] = True
                    continue

                if element.tag == "asphere":#-----------------------------------------
                    self.reset_params()
                    self.asphere_read = True
                    if self.invariant_data:
                        self.read_indicator['asphere'] = True
                    continue

                if element.tag == "patch":#-----------------------------------------
                    self.reset_params()
                    self.patch_read = True
                    if self.invariant_data:
                        self.read_indicator['patch'] = True
                    continue

                if element.tag == "bond":
                    self.reset_params()
                    self.bond_read = True
                    if self.invariant_data:
                        self.read_indicator['bond'] = True
                    if self.bond_read:
                        bond_content1 = "".join(element.text)
                        bond_content2 = bond_content1.strip().split()
                        start = 0
                        end = len(bond_content2)
                        step = 3
                        for i in range(start, end, step):
                            x = i
                            output_list = [int(j) if j.isdigit() else j for j in bond_content2[x:x + step]]
                            self.bond.append(output_list)
                    continue

                if element.tag == "angle":
                    self.reset_params()
                    self.angle_read = True
                    if self.invariant_data:
                        self.read_indicator['angle'] = True
                    if self.angle_read:
                        angle_content1 = "".join(element.text)
                        angle_content2 = angle_content1.strip().split()
                        start = 0
                        end = len(angle_content2)
                        step = 4
                        for i in range(start, end, step):
                            x = i
                            output_list = [int(j) if j.isdigit() else j for j in angle_content2[x:x + step]]
                            self.angle.append(output_list)
                    continue
					
                if element.tag == "dihedral":
                    self.reset_params()
                    self.dihedral_read = True
                    if self.invariant_data:
                        self.read_indicator['dihedral'] = True
                    if self.dihedral_read:
                        dihedral_content1 = "".join(element.text)
                        dihedral_content2 = dihedral_content1.strip().split()
                        start = 0
                        end = len(dihedral_content2)
                        step = 5
                        for i in range(start, end, step):
                            x = i
                            output_list = [int(j) if j.isdigit() else j for j in dihedral_content2[x:x + step]]
                            self.dihedral.append(output_list)
                    continue

                if element.tag == "vsite":
                    self.reset_params()
                    self.vsite_read = True
                    if self.invariant_data:
                        self.read_indicator['vsite'] = True
                    if self.vsite_read:
                        vsite_content1 = "".join(element.text)
                        vsite_content2 = vsite_content1.strip().split()
                        start = 0
                        end = len(vsite_content2)
                        step = 5
                        for i in range(start, end, step):
                            x = i
                            output_list = [int(j) if j.isdigit() else j for j in vsite_content2[x:x + step]]
                            self.vsite.append(output_list)
                    continue

        # check data
        if not self.xml_read:
            raise RuntimeError('Error! the input file is not a xml file with version 1.3!')

        if not self.CG_xml:
            if len(self.position) != int(self.num_particles):
                raise RuntimeError('Error! number of position ', len(self.position),
                                ' is not equal to the number of particles ', int(self.num_particles))
            print("info :", len(self.position), "positions")

            if len(self.image) > 0:
                if len(self.image) != int(self.num_particles):
                    raise RuntimeError('Error! number of image ', len(self.image),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.image), "images")

            if len(self.mass) > 0:
                if len(self.mass) != int(self.num_particles):
                    raise RuntimeError('Error! number of mass ', len(self.mass),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.mass), "masses")

            if len(self.type) != int(self.num_particles):
                raise RuntimeError('Error! number of type ', len(self.type),
                                ' is not equal to the number of particles ', int(self.num_particles))
            print("info :", len(self.type), "types")

            if len(self.velocity) > 0:
                if len(self.velocity) != int(self.num_particles):
                    raise RuntimeError('Error! number of velocity ', len(self.velocity),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.velocity), "velocities")

            if len(self.charge) > 0:
                if len(self.charge) != int(self.num_particles):
                    raise RuntimeError('Error! number of charge ', len(self.charge),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.charge), "charges")

            if len(self.body) > 0:
                if len(self.body) != int(self.num_particles):
                    raise RuntimeError('Error! number of body ', len(self.body),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.body), "bodies")

            if len(self.diameter) > 0:
                if len(self.diameter) != int(self.num_particles):
                    raise RuntimeError('Error! number of diameter ', len(self.diameter),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.diameter), "diameters")

            if len(self.rotangle) > 0:
                if len(self.rotangle) != int(self.num_particles):
                    raise RuntimeError('Error! number of rotangle ', len(self.rotangle),
                                    ' is not equal to the number of particles ', int(self.num_particles))
                print("info :", len(self.rotangle), "rotangles")                      

            if len(self.force) > 0:
                if len(self.force) != int(self.num_particles):
                    raise RuntimeError('Error! number of force ', len(self.force),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.force), "forces")

            if len(self.virial) > 0:
                if len(self.virial) != self.num_particles:
                    raise RuntimeError('Error! number of virial ', len(self.virial),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.virial), "virials")

            if len(self.molecule) > 0:
                if len(self.molecule) != self.num_particles:
                    raise RuntimeError('Error! number of molecule ', len(self.molecule),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.molecule), "molecules")

            if len(self.init) > 0:
                if len(self.init) != self.num_particles:
                    raise RuntimeError('Error! number of init ', len(self.init),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.init), "inits")

            if len(self.cris) > 0:
                if len(self.cris) != self.num_particles:
                    raise RuntimeError('Error! number of cris ', len(self.cris),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.cris), "crises")

            if len(self.orientation) > 0:
                if len(self.orientation) != self.num_particles:
                    raise RuntimeError('Error! number of orientation ', len(self.orientation),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.orientation), "orientations")

            if len(self.quaternion) > 0:
                if len(self.quaternion) != self.num_particles:
                    raise RuntimeError('Error! number of quaternion ', len(self.quaternion),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.quaternion), "quaternions")

            if len(self.rotation) > 0:
                if len(self.rotation) != self.num_particles:
                    raise RuntimeError('Error! number of rotation ', len(self.rotation),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.rotation), "rotations")

            if len(self.inert) > 0:
                if len(self.inert) != self.num_particles:
                    raise RuntimeError('Error! number of inert ', len(self.inert),
                                    ' is not equal to the number of particles ', self.num_particles)
                print("info :", len(self.inert), "inerts")

        if len(self.asphere) > 0:
            print("info :", len(self.asphere), "aspheres")

        if len(self.patch) > 0:
            print("info :", len(self.patch), "patches")

        if len(self.bond) > 0:
            print("info :", len(self.bond), "bonds")

        if len(self.angle) > 0:
            print("info :", len(self.angle), "angles")

        if len(self.dihedral) > 0:
            print("info :", len(self.dihedral), "dihedrals")

        if len(self.vsite) > 0:
            print("info :", len(self.vsite), "vsites")

    def init_data(self):
        # data
        if not self.read_indicator['box']:
            self.box = []
        if not self.read_indicator['position']:
            self.position = []
        if not self.read_indicator['type']:
            self.type = []
        if not self.read_indicator['image']:
            self.image = []
        if not self.read_indicator['mass']:
            self.mass = []
        if not self.read_indicator['velocity']:
            self.velocity = []
        if not self.read_indicator['charge']:
            self.charge = []
        if not self.read_indicator['body']:
            self.body = []
        if not self.read_indicator['diameter']:
            self.diameter = []
        if not self.read_indicator['rotangle']:
            self.rotangle = []
        if not self.read_indicator['force']:
            self.force = []
        if not self.read_indicator['virial']:
            self.virial = []
        if not self.read_indicator['molecule']:
            self.molecule = []
        if not self.read_indicator['init']:
            self.init = []
        if not self.read_indicator['cris']:
            self.cris = []
        if not self.read_indicator['orientation']:
            self.orientation = []
        if not self.read_indicator['quaternion']:
            self.quaternion = []
        if not self.read_indicator['rotation']:
            self.rotation = []
        if not self.read_indicator['inert']:
            self.inert = []
        if not self.read_indicator['asphere']:
            self.asphere = []
        if not self.read_indicator['patch']:
            self.patch = []
        if not self.read_indicator['bond']:
            self.bond = []
        if not self.read_indicator['angle']:
            self.angle = []
        if not self.read_indicator['dihedral']:
            self.dihedral = []
        if not self.read_indicator['vsite']:
            self.vsite = []

    # reset parameters
    def reset_params(self):
        # indicators
        self.num_particles_read = False
        self.timestep_read = False
        self.dimension_read = False
        self.box_read = False
        self.position_read = False
        self.type_read = False
        self.image_read = False
        self.mass_read = False
        self.velocity_read = False
        self.charge_read = False
        self.body_read = False
        self.diameter_read = False
        self.rotangle_read = False
        self.force_read = False
        self.virial_read = False
        self.molecule_read = False
        self.init_read = False
        self.cris_read = False
        self.orientation_read = False
        self.quaternion_read = False
        self.rotation_read = False
        self.inert_read = False
        self.asphere_read = False
        self.patch_read = False
        self.bond_read = False
        self.angle_read = False
        self.dihedral_read = False
        self.vsite_read = False