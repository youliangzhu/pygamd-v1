"""
code: XML
作用：读写XML文件
"""

from io import StringIO
import xml.etree.cElementTree as ET


class XML:
    def __init__(self, file_name=None):
        self.all_beads_point = None
        self.tree = ET.parse(file_name)
        self.filename = file_name
        self.read_indicator = {'box': False, 'position': False, 'type': False, 'image': False,
                               'velocity': False, 'charge': False, 'body': False, 'mass': False,
                               'diameter': False, 'rotangle': False, 'force': False, 'virial': False,
                               'init': False, 'cris': False, 'pair': False, 'molecule': False,
                               'orientation': False, 'quaternion': False, 'rotation': False, 
                               'asphere': False, 'patch': False, 'bond': False, 'inert': False,
                               'angle': False, 'dihedral': False, 'vsite': False, 'constraint': False}
        self.init_data()

    @staticmethod
    def __get_config(root):
        gala_node = root.getroot()
        configuration = list(gala_node)[0]
        return configuration

    def __record_all_info(self):
        for element in self.__get_config(self.tree):
            if element.tag == "box":
                self.box = [float(element.attrib["lx"]), float(element.attrib["ly"]), float(element.attrib["lz"])]
            elif element.tag == "position":
                self.position = extract_3d(self.position, element.text)
            elif element.tag == "velocity":
                self.velocity = extract_3d(self.velocity, element.text)
            elif element.tag == "body":
                self.body = extract_1d_int(self.body, element.text)
            elif element.tag == "h_cris":
                self.cris = extract_1d_int(self.cris, element.text)
            elif element.tag == "h_init":
                self.init = extract_1d_int(self.init, element.text)
            elif element.tag == "molecule":
                self.molecule = extract_1d_int(self.molecule, element.text)
            elif element.tag == "image":
                self.image = extract_3d_int(self.image, element.text)
            elif element.tag == "type":
                self.type = extract_1d_str(self.type, element.text)
            elif element.tag == "mass":
                self.mass = extract_1d_float(self.mass, element.text)
            elif element.tag == "charge":
                self.charge = extract_1d_float(self.charge, element.text)
            elif element.tag == "diameter":
                self.diameter = extract_1d_float(self.diameter, element.text)
            elif element.tag == "bond":
                self.bond = extract_3d_bond(self.bond, element.text)
            elif element.tag == "pair":
                self.pair = extract_3d_bond(self.pair, element.text)
            elif element.tag == "constraint":
                self.constraint = extract_3d_bond(self.constraint, element.text)
            elif element.tag == "angle":
                self.angle = extract_4d_angle(self.angle, element.text)
            elif element.tag == "vsite":
                self.vsite = extract_5d_dihedral(self.vsite, element.text)
            elif element.tag == "dihedral":
                self.dihedral = extract_5d_dihedral(self.dihedral, element.text)

    def __assign_record_info(self, assign):
        for element in self.__get_config(self.tree):
            if element.tag == "box" and element.tag in assign:
                self.box = [float(element.attrib["lx"]), float(element.attrib["ly"]), float(element.attrib["lz"])]
            elif element.tag == "position" and element.tag in assign:
                self.position = extract_3d(self.position, element.text)
            elif element.tag == "velocity":
                self.velocity = extract_3d(self.velocity, element.text)
            elif element.tag == "image" and element.tag in assign:
                self.image = extract_3d_int(self.image, element.text)
            elif element.tag == "type" and element.tag in assign:
                self.type = extract_1d_str(self.type, element.text)
            elif element.tag == "bond" and element.tag in assign:
                self.bond = extract_3d_bond(self.bond, element.text)
            elif element.tag == "mass" and element.tag in assign:
                self.mass = extract_1d_float(self.mass, element.text)

    def return_info(self):
        self.__record_all_info()
        info = {"box": self.box,
                "position": self.position,
                "type": self.type,
                "mass": self.mass,
                "diameter": self.diameter,
                "charge": self.charge,
                "molecule": self.molecule,
                "init": self.init,
                "cris": self.cris,
                "bond": self.bond,
                "pair": self.pair,
                "constraint": self.constraint,
                "vsite": self.vsite,
                "body": self.body,
                "image": self.image,
                "velocity": self.velocity,
                "angle": self.angle,
                "dihedral": self.dihedral,
                "time_step": self.__get_config(self.tree).attrib['time_step'],
                }

        return info

    def assign_info(self, assign=None):
        if assign is None:
            assign = []

        self.__assign_record_info(assign)
        info = {}
        for assi in assign:
            match assi:
                case 'box':
                    info[assi] = self.box
                case 'position':
                    info[assi] = self.position
                case 'image':
                    info[assi] = self.image
                case 'velocity':
                    info[assi] = self.velocity
                case 'type':
                    info[assi] = self.type
                case 'bond':
                    info[assi] = self.bond
                case 'mass':
                    info[assi] = self.mass
                case _:
                    pass

        info["time_step"] = self.__get_config(self.tree).attrib['time_step']
        return info

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
        if not self.read_indicator['pair']:
            self.pair = []
        if not self.read_indicator['constraint']:
            self.constraint = []


def extract_3d(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()

    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        element_line = element_line.replace('\t', ' ').split(' ')
        tem = []
        for element_i in element_line:
            if element_i:
                tem.append(float(element_i))
        selfAim.append(tem)
    return selfAim

def extract_3d_int(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()

    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        element_line = element_line.replace('\t', ' ').split(' ')
        tem = []
        for element_i in element_line:
            if element_i:
                tem.append(int(element_i))
        selfAim.append(tem)
    return selfAim


def extract_1d_str(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()
    
    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        selfAim.append(element_line)
    return selfAim


def extract_1d_float(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()
    
    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        selfAim.append(float(element_line))
    return selfAim


def extract_1d_int(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()
    
    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        selfAim.append(int(element_line))
    return selfAim


def extract_3d_bond(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()
    
    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        element_line = element_line.replace('\t', ' ').split(' ')
        tem = []
        for element_i in element_line:
            if element_i:
                tem.append(element_i)
        selfAim.append([tem[0], int(tem[1]), int(tem[2])])
    
    return selfAim


def extract_4d_angle(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()
    
    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        element_line = element_line.replace('\t', ' ').split(' ')
        tem = []
        for element_i in element_line:
            if element_i:
                tem.append(element_i)
        selfAim.append([tem[0], int(tem[1]), int(tem[2]), int(tem[3])])
    return selfAim


def extract_5d_dihedral(selfAim, elementText):
    element_IO = StringIO(elementText)
    element_IO.readline()

    while True:
        element_line = element_IO.readline()
        element_line = element_line.strip()
        if element_line == "":
            break
        element_line = element_line.replace('\t', ' ').split(' ')
        tem = []
        for element_i in element_line:
            if element_i:
                tem.append(element_i)
        selfAim.append([tem[0], int(tem[1]), int(tem[2]), int(tem[3]), int(tem[4])])
    return selfAim


"""
输出xml文件
"""

def backimage(posi, imagei, box):
    posi[0] += imagei[0] * box[0]
    posi[1] += imagei[1] * box[1]
    posi[2] += imagei[2] * box[2]


def reimage(posi, box):
    imagei = [0, 0, 0]
    dx = posi[0]
    dy = posi[1]
    dz = posi[2]
    roundx = round(dx / box[0])
    roundy = round(dy / box[1])
    roundz = round(dz / box[2])
    dx -= box[0] * roundx
    dy -= box[1] * roundy
    dz -= box[2] * roundz
    posi[0] = dx
    posi[1] = dy
    posi[2] = dz
    imagei[0] = int(roundx)
    imagei[1] = int(roundy)
    imagei[2] = int(roundz)
    return imagei


def View_S(filename, pos, types=None, box=None, bond=None, pair=None, constraint=None, velocity=None, body=None, image=None,
           init=None, cris=None, mass=None, charge=None, diameter=None, angle=None, dihedral=None, vsite=None,
           dimensions=3, time_step=0):
    if box is None:
        box = [1, 1, 1]

    if len(pos) != len(types):
        print(len(pos), len(types), "pos_len != type_len")
        return -1

    opd = open(filename, "w")
    opd.write('<?xml version ="1.0" encoding ="UTF-8" ?>\n<galamost_xml version="1.3">\n')
    opd.write('<configuration time_step="' + str(time_step) + '" dimensions="' + str(dimensions) + '" natoms="' +
              str(len(pos)) + '" >\n')
    opd.write('<box lx="' + str(box[0]) + '" ly="' + str(box[1]) + '" lz="' + str(box[2]) + '"/>\n')

    opd.write('<position num="' + str(len(pos)) + '">\n')
    for i in pos:
        opd.write(str(format(i[0], '17.10f')) + " " + str(format(i[1], '17.10f')) + " " + str(
            format(i[2], '17.10f')) + "\n")
    opd.write('</position>\n')

    opd.write('<type num="' + str(len(types)) + '">\n')
    for i in types:
        opd.write(i + '\n')
    opd.write('</type>\n')

    if velocity is not None:
        if len(velocity) > 0:
            opd.write('<velocity num="' + str(len(velocity)) + '">\n')
            for i in velocity:
                opd.write(str(format(i[0], '17.10f')) + " " + str(format(i[1], '17.10f')) + " " + str(
                    format(i[2], '17.10f')) + "\n")
            opd.write('</velocity>\n')

    if init is not None:
        if len(init) > 0:
            opd.write('<h_init num="' + str(len(init)) + '">\n')
            for i in init:
                opd.write(str(i) + "\n")
            opd.write('</h_init>\n')

    if cris is not None:
        if len(cris) > 0:
            opd.write('<h_cris num="' + str(len(cris)) + '">\n')
            for i in cris:
                opd.write(str(i) + "\n")
            opd.write('</h_cris>\n')

    if mass is not None:
        if len(mass) > 0:
            opd.write('<mass num="' + str(len(mass)) + '">\n')
            for i in mass:
                opd.write(str(i) + "\n")
            opd.write('</mass>\n')

    if charge is not None:
        if len(charge) > 0:
            opd.write('<charge num="' + str(len(charge)) + '">\n')
            for i in charge:
                opd.write(str(i) + "\n")
            opd.write('</charge>\n')
            
    if diameter is not None:
        if len(diameter) > 0:
            opd.write('<diameter num="' + str(len(diameter)) + '">\n')
            for i in diameter:
                opd.write(str(i) + "\n")
            opd.write('</diameter>\n')

    if body is not None:
        if len(body) > 0:
            opd.write('<body num="' + str(len(body)) + '">\n')
            for i in body:
                opd.write(str(i) + "\n")
            opd.write('</body>\n')

    if image is not None:
        if len(image) > 0:
            opd.write('<image num="' + str(len(image)) + '">\n')
            for i in image:
                opd.write(str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
            opd.write('</image>\n')

    if bond is not None and len(bond):
        opd.write('<bond num="' + str(len(bond)) + '">\n')
        for i in bond:
            opd.write(
                str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
        opd.write('</bond>\n')
    
    if pair is not None and len(pair):
        opd.write('<pair num="' + str(len(pair)) + '">\n')
        for i in pair:
            opd.write(
                str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
        opd.write('</pair>\n')
    
    if constraint is not None and len(constraint):
        opd.write('<constraint num="' + str(len(constraint)) + '">\n')
        for i in constraint:
            opd.write(
                str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
        opd.write('</constraint>\n')

    if angle is not None:
        if len(angle) > 0:
            opd.write('<angle num="' + str(len(angle)) + '">\n')
            for i in angle:
                opd.write(
                    str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + " " + str(i[3]) + "\n")
            opd.write('</angle>\n')

    if dihedral is not None:
        if len(dihedral) > 0:
            opd.write('<dihedral num="' + str(len(dihedral)) + '">\n')
            for i in dihedral:
                opd.write(
                    str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + " " + str(i[3]) + " " + str(i[4]) + "\n")
            opd.write('</dihedral>\n')

    if vsite is not None:
        if len(vsite) > 0:
            opd.write('<vsite num="' + str(len(vsite)) + '">\n')
            for i in vsite:
                opd.write(
                    str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + " " + str(i[3]) + " " + str(i[4]) + "\n")
            opd.write('</vsite>\n')

    opd.write('</configuration>\n</galamost_xml>\n')
    opd.close()


def View_(filename, infos):
    oxml = {"position": None, "type": None, "box": None, "bond": None, "pair": None, "constraint": None, "velocity": None, "body": None, "image": None,
            "init": None, "cris": None, "mass": None, "charge": None, "diameter": None, "angle": None, "dihedral": None, "vsite": None,
            "dimensions": 3, "time_step": 0}

    for key in infos:
        oxml[key] = infos[key]
        
    View_S(filename,
           pos=oxml["position"],
           types=oxml["type"],
           box=oxml["box"],
           bond=oxml["bond"],
           pair=oxml["pair"],
           constraint=oxml["constraint"],
           velocity=oxml["velocity"],
           body=oxml["body"],
           image=oxml["image"],
           init=oxml["init"],
           cris=oxml["cris"],
           mass=oxml["mass"],
           charge=oxml["charge"],
           diameter=oxml["diameter"],
           angle=oxml["angle"],
           dihedral=oxml["dihedral"],
           vsite=oxml["vsite"],
           dimensions=oxml["dimensions"],
           time_step=oxml["time_step"])


def newName(f_i, head_n):
    f_i = f_i.split('\\')
    f_i[-1] = head_n + f_i[-1]
    new_f_i = '\\'.join(f_i)
    return new_f_i
