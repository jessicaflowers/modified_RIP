import os
import shutil
import glob
import time
import random
import math
import copy
import sys
import re

parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent not in sys.path:
    sys.path.insert(0, parent)

import pdbstruct
from pdbstruct import vector3d

# Debugging: enable with environment variable RIP_DEBUG=1
DEBUG = os.environ.get('RIP_DEBUG', '0') == '1'
# limit printed debug lines to avoid flooding
_ATOM_DEBUG_LIMIT = int(os.environ.get('RIP_DEBUG_LIMIT', '100'))
_ATOM_DEBUG_COUNT = 0

# Defaults 
TOP = "../input/topol.top"
# MIN_OUTPUT = "seg_01_MIN"   # base name (without .gro)
MIN_OUTPUT = "seg_11_NPT_PROD_01"  # starting from production phase, rather than minimized strucutre
TEMP_PERTURB = 300.0        # K
OUT_NAME = "seg_11_NPT_PROD_01_PERT.gro"
BOX_LINE = None

def maxwell_velocity(temp, mass):
    # temp in K, mass in a.m.u., output in angstroms/ps
    factor = 8314.47148 * temp / mass  # in (m/s)^2
    convert_to_ang_per_ps = 0.01  # from m/s to angstroms/ps
    return random.gauss(0.0, math.sqrt(factor)) * convert_to_ang_per_ps

def mean_energy(temp, n_degree_of_freedom):
    boltzmann = 8314.47148  # Da (m/s)^2 K^-1
    convert_to_ang_per_ps = 1.0E-4  # (m/s)^2 to (angstroms/ps)^2
    return n_degree_of_freedom * 0.5 * boltzmann * temp * convert_to_ang_per_ps

def random_energy(temp, n_degree_of_freedom):
    average = mean_energy(temp, n_degree_of_freedom)
    std_dev = math.sqrt(average)
    return random.gauss(average, std_dev)

def kinetic_energy(atoms):
    en = 0.0
    for a in atoms:
        v2 = a.vel.length() ** 2
        en += 0.5 * a.mass * v2
    return en

def gas_randomize(atoms, temp):
    # set Maxwellian velocities for each atom (Å/ps)
    for a in atoms:
        a.vel.set(
            maxwell_velocity(temp, a.mass),
            maxwell_velocity(temp, a.mass),
            maxwell_velocity(temp, a.mass)
        )

def anderson_velocity_scale(atoms, temp, n_degree_of_freedom):
    target_energy = mean_energy(temp, n_degree_of_freedom)
    kin = kinetic_energy(atoms)
    if vector3d.is_near_zero(kin):
        # no kinetic energy -> randomize then scale to target
        gas_randomize(atoms, temp)
        kin = kinetic_energy(atoms)
        if vector3d.is_near_zero(kin):
            return
    scale = math.sqrt(target_energy / kin)
    for a in atoms:
        a.vel.scale(scale)
    # remove center-of-mass velocity
    total_mass = sum(a.mass for a in atoms)
    if total_mass > 0.0:
        com = vector3d.Vector3d()
        for a in atoms:
            tmp = a.vel.copy()
            tmp.scale(a.mass)
            com += tmp
        com.scale(1.0 / total_mass)
        for a in atoms:
            a.vel -= com

# --- rotational / chi utilities ---
def moment_of_inertia(atom, axis, anchor):
    r = atom.pos - anchor
    r_perp = r.perpendicular_vec(axis)
    r_len = r_perp.length()
    return atom.mass * (r_len ** 2)

def total_moment_of_inertia(atoms, axis, anchor):
    return sum(moment_of_inertia(a, axis, anchor) for a in atoms)

def rotational_velocity(atom, axis, anchor):
    # Estimate local angular speed omega for this atom about axis through anchor.
    r = atom.pos - anchor
    r_perp = r.perpendicular_vec(axis)
    r_len = r_perp.length()
    if vector3d.is_near_zero(r_len):
        return 0.0
    u = axis.copy()
    u.normalize()
    pos_ref = vector3d.CrossProductVec(u, r_perp)  # magnitude = r_len
    # omega = dot(v, pos_ref) / (r_len^2)
    denom = r_len * r_len
    if vector3d.is_near_zero(denom):
        return 0.0
    return vector3d.dot(atom.vel, pos_ref) / denom

def weighted_rotational_velocity(atoms, axis, anchor):
    moments = [moment_of_inertia(a, axis, anchor) for a in atoms]
    total_moment = sum(moments)
    if vector3d.is_near_zero(total_moment):
        return 0.0
    weights = [m / total_moment for m in moments]
    rot_vels = [rotational_velocity(a, axis, anchor) for a in atoms]
    return sum(rv * w for rv, w in zip(rot_vels, weights))

def add_rotational_velocity(atoms, rot_vel, axis, anchor):
    # Adds linear velocity v = Omega x r to each atom (Omega = rot_vel * axis_unit)
    u = axis.copy()
    u.normalize()
    for a in atoms:
        r = a.pos - anchor
        # r_perp = component of r perpendicular to axis
        r_perp = r.perpendicular_vec(u)
        # v_add = rot_vel * (u x r_perp)
        v_add = vector3d.CrossProductVec(u, r_perp)
        v_add.scale(rot_vel)
        a.vel += v_add

def get_axis_anchor(res, i):
    # returns (axis_vector, anchor_pos_vector)
    chi_topology = pdbstruct.get_res_chi_topology(res.type)
    # chi_topology[i] is a list of atom names defining that chi
    p = [res.atom(atom_type).pos for atom_type in chi_topology[i]]
    axis = p[2] - p[1]
    anchor = p[2].copy()
    return axis, anchor

def atoms_affected_by_chi(atoms, i):
    return [atom for atom in atoms if pdbstruct.get_atom_sidechain_nesting(atom.type) >= i]

def get_rot_vel_chi(res, i):
    axis, anchor = get_axis_anchor(res, i)
    atoms = atoms_affected_by_chi(res.atoms(), i)
    return weighted_rotational_velocity(atoms, axis, anchor)

def get_random_chi_rot_vel(res, i, temp):
    axis, anchor = get_axis_anchor(res, i)
    atoms = atoms_affected_by_chi(res.atoms(), i)
    moment = total_moment_of_inertia(atoms, axis, anchor)
    n_atom = len(atoms)
    if vector3d.is_near_zero(moment) or n_atom == 0:
        return 0.0
    energy = random_energy(temp, 3 * n_atom)
    return math.sqrt(max(0.0, 2.0 * energy / moment))

def add_rot_vel_to_chi(res, i, target_rot_vel):
    axis, anchor = get_axis_anchor(res, i)
    atoms = atoms_affected_by_chi(res.atoms(), i)
    add_rotational_velocity(atoms, target_rot_vel, axis, anchor)

def randomize_clean_high_chi_with_attractor(
    res, temp, mean_chi_values,
    max_delta_chi=vector3d.DEG2RAD * 45.0):
    "Chi-pulses the res at temp. Vel: angstroms per picosecond."
    n_chi = pdbstruct.get_n_chi(res)
    rot_vels = [get_rot_vel_chi(res, i) for i in range(n_chi)]

    # clear all velocities of this residue
    for atom in res.atoms():
        atom.vel.set(0.0, 0.0, 0.0)

    for i_chi in reversed(list(range(n_chi))):
        # decide sign/direction using attractor logic
        chi = pdbstruct.calculate_chi(res, i_chi)
        delta_chi = vector3d.normalize_angle(chi - mean_chi_values[i_chi])
        if delta_chi > max_delta_chi:
            sign = -1.0
        elif delta_chi < -max_delta_chi:
            sign = 1.0
        else:
            sign = 1.0 if rot_vels[i_chi] > 0.0 else -1.0

        target_rot_vel = sign * get_random_chi_rot_vel(res, i_chi, temp)
        add_rot_vel_to_chi(res, i_chi, target_rot_vel)

    anderson_velocity_scale(res.atoms(), temp, 3 * len(res.atoms()))

def read_top(top):
    # Returns: (masses_list, total_qtot, atomtype_mass_map)
    lines = open(top).readlines()
    masses = []
    total_qtot = 0.0
    atomtype_map = {}
    top_dir = os.path.dirname(top)

    i = 0
    nlines = len(lines)
    while i < nlines:
        l = lines[i]
        # process includes immediately
        if l.strip().startswith('#include'):
            print("Processing include:", l.strip())
            parts = l.split()
            if len(parts) > 1:
                itp = parts[1].strip()
                # remove quotes if present
                if itp.startswith('"') and itp.endswith('"'):
                    itp = itp[1:-1]
                itp_path = os.path.join(top_dir, itp)
                if os.path.isfile(itp_path):
                    inc_masses, inc_qtot, inc_atommap = read_top(itp_path)
                    masses.extend(inc_masses)
                    total_qtot += inc_qtot
                    # merge atomtype maps (later entries override earlier)
                    atomtype_map.update(inc_atommap)
        # parse [ atomtypes ] blocks
        if l.strip().lower().startswith('[ atomtypes ]'):
            # skip header line(s) until non-comment/blank
            i += 1
            while i < nlines:
                l2 = lines[i]
                if l2.strip().startswith('['):
                    i -= 1
                    break
                if l2.strip() and not l2.strip().startswith(';'):
                    words = l2.split()
                    # expected format: name at.num mass charge ptype sigma epsilon
                    if len(words) >= 3:
                        name = words[0]
                        try:
                            mass = float(words[2])
                            atomtype_map[name] = mass
                        except Exception:
                            pass
                i += 1
        # parse [ atoms ] blocks for explicit per-atom masses
        if l.strip().lower().startswith('[ atoms ]'):
            i += 1
            qtot = None
            while i < nlines:
                l2 = lines[i]
                if l2.strip().startswith('['):
                    i -= 1
                    break
                if l2.strip() and not l2.strip().startswith(';'):
                    words = l2.split()
                    try:
                        n = int(words[0])
                    except Exception:
                        i += 1
                        continue
                    mass = None
                    if len(words) > 7:
                        try:
                            mass = float(words[7])
                        except Exception:
                            mass = None
                    if mass is None:
                        for w in reversed(words):
                            try:
                                mass = float(w)
                                break
                            except Exception:
                                continue
                    if mass is not None:
                        masses.append(mass)
                    # check for per-residue qtot annotation
                    if ';' in l2 and 'qtot' in l2:
                        try:
                            tail = l2.split('qtot')[-1].strip()
                            q_line = float(tail.split()[0])
                            if q_line is not None:
                                total_qtot += q_line
                        except Exception:
                            pass
                i += 1
        i += 1

    return masses, total_qtot, atomtype_map

def AtomFromGroLine(line):
    atom = pdbstruct.Atom()
    atom.res_num = int(line[0:5])
    atom.res_type = line[5:10].strip()
    atom.type = line[10:15].strip()
    # print(f'atom.res_num: {atom.res_num}, atom.res_type: {atom.res_type}, atom.type: {atom.type}')
    if atom.res_type == "ILE" and atom.type == "CD":
        atom.type = "CD1"
    element = ''
    for c in line[12:15]:
        if not c.isdigit() and c != " ":
            element += c
    if element[:2] in pdbstruct.two_char_elements:
        atom.element = element[:2]
    else:
        atom.element = element[0] if element else ''
    atom.num = int(line[15:20])
    x = 10.0 * float(line[20:28])
    y = 10.0 * float(line[28:36])
    z = 10.0 * float(line[36:44])
    atom.pos.set(x, y, z)
    if len(line) > 62:
        x = 10.0 * float(line[44:52])
        y = 10.0 * float(line[52:60])
        z = 10.0 * float(line[60:68])
        atom.vel.set(x, y, z)
    else:
        atom.vel.set(0.0, 0.0, 0.0)
    # debug printing: show raw slices and parsed fields (limited)
    global _ATOM_DEBUG_COUNT
    if DEBUG and _ATOM_DEBUG_COUNT < _ATOM_DEBUG_LIMIT:
        _ATOM_DEBUG_COUNT += 1
        # show a short repr of the raw line (first 80 chars)
        short = line.rstrip('\n')
        if len(short) > 160:
            short = short[:160] + '...'
        print("[RIP_DEBUG] AtomFromGroLine:")
        print(f"  raw[{len(line)}]: {short!r}")
        print(f"  slices: resnum_slice={line[0:5]!r}, resname_slice={line[5:10]!r}, atomname_slice={line[10:15]!r}, atomnum_slice={line[15:20]!r}")
        print(f"  parsed -> res_num={atom.res_num!r}, res_type={atom.res_type!r}, type={atom.type!r}, num={atom.num!r}")
    return atom


def parse_included_molecule_masses(top):
    """Parse included .itp files for per-molecule [ atoms ] blocks.
    Returns dict {(moleculetype_name, atom_name): mass}
    """
    base = os.path.dirname(top)
    mol_atom_map = {}
    seen = set()
    try:
        lines = open(top).readlines()
    except Exception:
        return mol_atom_map
    includes = []
    for l in lines:
        l = l.strip()
        if l.startswith('#include'):
            parts = l.split()
            if len(parts) > 1:
                path = parts[1].strip().strip('"')
                includes.append(os.path.join(base, path))
    # recursively scan included files
    queue = list(includes)
    while queue:
        f = queue.pop(0)
        if not os.path.isfile(f) or f in seen:
            continue
        seen.add(f)
        try:
            flines = open(f).readlines()
        except Exception:
            continue
        n = len(flines)
        i = 0
        cur_mol = None
        while i < n:
            l = flines[i].strip()
            if l.lower().startswith('#include'):
                parts = l.split()
                if len(parts) > 1:
                    path = parts[1].strip().strip('"')
                    queue.append(os.path.join(os.path.dirname(f), path))
            if l.lower().startswith('[ moleculetype ]'):
                # next non-comment non-empty line is moleculetype name
                i += 1
                while i < n and (not flines[i].strip() or flines[i].strip().startswith(';')):
                    i += 1
                if i < n:
                    cur_mol = flines[i].split()[0]
            if l.lower().startswith('[ atoms ]'):
                i += 1
                while i < n:
                    l2 = flines[i]
                    if l2.strip().startswith('['):
                        i -= 1
                        break
                    if l2.strip() and not l2.strip().startswith(';'):
                        words = l2.split()
                        atom_name = None
                        mass = None
                        if len(words) >= 5:
                            atom_name = words[4]
                        for w in reversed(words):
                            try:
                                mass = float(w)
                                break
                            except Exception:
                                continue
                        if atom_name and mass is not None and cur_mol:
                            # normalize names: uppercase and strip whitespace
                            atom_key = atom_name.strip().upper()
                            mol_key = cur_mol.strip().upper()
                            mol_atom_map[(mol_key, atom_key)] = mass
                            # also store a base-atom name without trailing digits (HB1 -> HB)
                            base_atom = re.sub(r"[^A-Z0-9]", "", atom_key).rstrip('0123456789')
                            if base_atom and base_atom != atom_key:
                                mol_atom_map[(mol_key, base_atom)] = mass
                            # many CHARMM-GUI itps use a single moleculetype (e.g. PROC)
                            # but include a `residu` column giving the residue name (ALA, LEU, ...)
                            # when available, also index by that residue name so lookups by
                            # the GRO's residue (e.g. 'LEU') will succeed.
                            # words layout: nr type resnr residu atom cgnr charge mass
                            residu = None
                            if len(words) >= 5:
                                # words[3] is residu in typical GROMACS itp [ atoms ] format
                                try:
                                    residu = words[3]
                                except Exception:
                                    residu = None
                            if residu:
                                residu_key = residu.strip().upper()
                                mol_atom_map[(residu_key, atom_key)] = mass
                                if base_atom and base_atom != atom_key:
                                    mol_atom_map[(residu_key, base_atom)] = mass
                    i += 1
            i += 1
    return mol_atom_map

def SoupFromGromacs(top, gro, skip_solvent=True):
    atoms = []
    lines = open(gro, 'r').readlines()
    remaining_text = ""
    n_remaining_text = 0
    for i_line, line in enumerate(lines[2:-1]):
        # print('\n')
        # print(f'line: {line}')
        atom = AtomFromGroLine(line)
        # print(f'atom: {atom}')
        # print('\n')
        if skip_solvent and atom.res_type == "SOL":
            remaining_text = "".join(lines[i_line+2:-1])
            n_remaining_text = len(lines[i_line+2:-1])
            break
        atoms.append(atom)
    box = [float(w) for w in lines[-1].split()]
    masses, q_tot, atomtype_map = read_top(top)
    # also parse per-molecule [ atoms ] blocks from included itps
    mol_atom_map = parse_included_molecule_masses(top)
    # assign per-atom masses when available (explicit per-atom list)
    for a, mass in zip(atoms, masses):
        a.mass = mass
    # if per-atom masses incomplete, try to use per-molecule atom mass map first
    if len(masses) < len(atoms) and mol_atom_map:
        n_assigned_mol = 0
        for a in atoms:
            if getattr(a, 'mass', 0.0):
                continue
            res = getattr(a, 'res_type', '') or ''
            atype = getattr(a, 'type', '') or ''
            # normalize to uppercase and strip non-alphanumeric chars
            res_u = re.sub(r"[^A-Z0-9]", "", res.upper())
            at_u = re.sub(r"[^A-Z0-9]", "", atype.upper())
            # try several residue-name variants to match moleculetype keys
            res_variants = []
            res_variants.append(res_u)
            res_variants.append(res_u.rstrip('0123456789'))
            # try common alternate like TIP -> TIP3 and TIP3 -> TIP
            if not res_u.endswith('3'):
                res_variants.append(res_u + '3')
            if res_u.endswith('3'):
                res_variants.append(res_u[:-1])
            # unique
            res_variants = [r for r in dict.fromkeys(res_variants) if r]
            mass = None
            for rv in res_variants:
                rv_norm = re.sub(r"[^A-Z0-9]", "", rv)
                if (rv_norm, at_u) in mol_atom_map:
                    mass = mol_atom_map[(rv_norm, at_u)]
                    break
            # also try stripping numeric suffix from atom name (HB1 -> HB)
            if mass is None:
                base_at = re.sub(r"[^A-Z0-9]", "", at_u).rstrip('0123456789')
                for rv in res_variants:
                    rv_norm = re.sub(r"[^A-Z0-9]", "", rv)
                    if (rv_norm, base_at) in mol_atom_map:
                        mass = mol_atom_map[(rv_norm, base_at)]
                        break
            # water/ion quick aliases
            if mass is None:
                # map common water atom names
                if at_u in ('OH2', 'OW'):
                    for rv in res_variants:
                        rv_norm = re.sub(r"[^A-Z0-9]", "", rv)
                        if (rv_norm, 'OT') in mol_atom_map:
                            mass = mol_atom_map[(rv_norm, 'OT')]
                            break
                if at_u in ('H1', 'H2', 'HW'):
                    for rv in res_variants:
                        rv_norm = re.sub(r"[^A-Z0-9]", "", rv)
                        if (rv_norm, 'HT1') in mol_atom_map:
                            mass = mol_atom_map[(rv_norm, 'HT1')]
                            break
            # last-resort: try any mol_atom_map entry matching the base atom name
            if mass is None:
                base_at = re.sub(r"[^A-Z0-9]", "", at_u).rstrip('0123456789')
                if base_at and base_at != at_u:
                    for (molname, aname), m in mol_atom_map.items():
                        if aname == base_at:
                            mass = m
                            break
            if mass is not None:
                a.mass = mass
                n_assigned_mol += 1
        if n_assigned_mol:
            print(f"Info: assigned masses from molecule [ atoms ] blocks to {n_assigned_mol} atoms")
    # if per-atom masses still incomplete, try to use atom type -> mass map
    if len(masses) < len(atoms) and atomtype_map:
        n_assigned = 0
        for a in atoms:
            if getattr(a, 'mass', 0.0):
                continue
            key = getattr(a, 'type', '')
            # try direct lookup, then try common variants
            mass = atomtype_map.get(key)
            if mass is None:
                # uppercase lookup
                mass = atomtype_map.get(key.upper())
            if mass is None:
                # strip trailing digits (HB1 -> HB)
                k2 = key.rstrip('0123456789')
                if k2 and k2 != key:
                    mass = atomtype_map.get(k2) or atomtype_map.get(k2.upper())
            if mass is None:
                # sometimes atomtypes use alternate naming (CD1 vs CD)
                if key.endswith('1'):
                    mass = atomtype_map.get(key[:-1]) or atomtype_map.get((key[:-1]).upper())
            if mass is not None:
                a.mass = mass
                n_assigned += 1
        if n_assigned:
            print(f"Info: assigned masses from atomtypes to {n_assigned} atoms")
    # assign reasonable default masses for any atoms not covered by the topology parse
    default_mass_table = {
        'H': 1.008,
        'HE': 4.0026,
        'LI': 6.94,
        'C': 12.011,
        'N': 14.007,
        'O': 15.999,
        'F': 18.998,
        'P': 30.974,
        'S': 32.06,
        'CL': 35.45,
        'BR': 79.904,
        'I': 126.90
    }
    n_defaulted = 0
    for a in atoms:
        if not getattr(a, 'mass', 0.0):
            el = (getattr(a, 'element', '') or '').upper()
            mass = default_mass_table.get(el)
            if mass is None:
                # fallback generic mass (carbon-like)
                mass = 12.011
            a.mass = mass
            n_defaulted += 1
    if n_defaulted:
        print(f"Warning: assigned default masses to {n_defaulted} atoms (missing in topology).")
    soup = pdbstruct.Polymer()
    curr_res_num = -1
    for a in atoms:
        if curr_res_num != a.res_num:
            res = pdbstruct.Residue(a.res_type, a.chain_id, a.res_num)
            soup.append_residue_no_renum(res.copy())
            curr_res_num = a.res_num
        soup.insert_atom(-1, a)
    soup.box = box
    soup.remaining_text = remaining_text
    soup.n_remaining_text = n_remaining_text
    return soup

# --- writer: writes gro with velocities columns (positions in nm, velocities in nm/ps) ---
def write_gro(soup, filename):
    atoms = []
    try:
        residues = list(soup.residues())
        for r in residues:
            for a in r.atoms():
                atoms.append(a)
    except Exception:
        # fallback: try numeric indices
        atoms = []
        i = 0
        while True:
            try:
                r = soup.residue(i)
                for a in r.atoms():
                    atoms.append(a)
                i += 1
            except Exception:
                break
    with open(filename, 'w') as fh:
        fh.write("Perturbed structure\n")
        fh.write(f"{len(atoms):5d}\n")
        for a in atoms:
            resnum = a.res_num
            resname = a.res_type[:5].ljust(5)
            atomname = a.type[:5].rjust(5)
            atomnum = a.num
            x = a.pos.x / 10.0
            y = a.pos.y / 10.0
            z = a.pos.z / 10.0
            vx = a.vel.x / 10.0
            vy = a.vel.y / 10.0
            vz = a.vel.z / 10.0
            line = f"{resnum:5d}{resname}{atomname}{atomnum:5d}{x:8.3f}{y:8.3f}{z:8.3f}{vx:8.4f}{vy:8.4f}{vz:8.4f}\n"
            fh.write(line)
        # box
        box = getattr(soup, "box", None)
        if box is None:
            box = [1.0, 1.0, 1.0]
        fh.write(" ".join(f"{b:8.5f}" for b in box) + "\n")

# --- main per-residue perturbation loop ---
def ensure_dir(d):
    if not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)

def load_soup(top, gro):
    return SoupFromGromacs(top, gro, skip_solvent=True)

def get_residue_list(soup):
    try:
        return list(soup.residues())
    except Exception:
        # fallback to soup.residue(i) iteration
        res_list = []
        i = 0
        while True:
            try:
                r = soup.residue(i)
                res_list.append(r)
                i += 1
            except Exception:
                break
        return res_list

def main(top=TOP, min_base=MIN_OUTPUT, temp=TEMP_PERTURB, target_resnums=None):
    gro = min_base + ".gro"
    if not os.path.isfile(gro):
        print("Minimized gro not found:", gro)
        return
    if not os.path.isfile(top):
        print("Topology not found:", top)
        return

    # initial soup used only to find list of rotatable residues
    master_soup = load_soup(top, gro)
    residues = get_residue_list(master_soup)

    # normalize target_resnums into a set of ints (or None to process all)
    target_set = None
    if target_resnums is not None:
        if isinstance(target_resnums, (list, set)):
            try:
                target_set = set(int(x) for x in target_resnums)
            except Exception:
                target_set = None
        else:
            try:
                target_set = set(int(x) for x in str(target_resnums).split(',') if x.strip())
            except Exception:
                target_set = None

    for res in residues:
        nchi = pdbstruct.get_n_chi(res)
        if nchi <= 0:
            continue
        if res.type == "PRO":
            continue
        # use residue number for directory naming (prefer res.num if present)
        resnum = getattr(res, "num", None)
        if resnum is None:
            # fallback to position in list
            resnum = residues.index(res)
        # if the caller requested a specific residue(s), skip others
        if target_set is not None and resnum not in target_set:
            continue
        out_dir = f"pulse_res_{resnum}"
        if os.path.isdir(out_dir):
            print("Skipping existing:", out_dir)
            continue
        print("Creating perturbation for residue", resnum, "->", out_dir)
        ensure_dir(out_dir)
        shutil.copy(top, os.path.join(out_dir, os.path.basename(top)))

        # reload a fresh copy of the minimized structure
        soup = load_soup(top, gro)
        # find the corresponding residue object in this soup (match by residue num)
        target_res = None
        try:
            for r in soup.residues():
                if getattr(r, "num", None) == resnum:
                    target_res = r
                    break
        except Exception:
            # fallback iteration
            i = 0
            while True:
                try:
                    r = soup.residue(i)
                    if getattr(r, "num", None) == resnum:
                        target_res = r
                        break
                    i += 1
                except Exception:
                    break
        if target_res is None:
            print("Could not find residue", resnum, "in reloaded soup; skipping")
            continue

        # compute mean chis from minimized structure (use values from master_soup res)
        master_res = res
        n_chi = pdbstruct.get_n_chi(master_res)
        mean_chis = [pdbstruct.calculate_chi(master_res, j) for j in range(n_chi)]

        # apply single pulse
        randomize_clean_high_chi_with_attractor(target_res, temp, mean_chis,
                                               max_delta_chi=60.0 * vector3d.DEG2RAD)

        # write perturbed gro
        out_gro = os.path.join(out_dir, OUT_NAME)
        write_gro(soup, out_gro)
        print("Wrote:", out_gro)

if __name__ == "__main__":
    # allow overrides from CLI: top gro temp
    top = TOP
    min_base = MIN_OUTPUT
    temp = TEMP_PERTURB
    if len(sys.argv) > 1:
        min_base = sys.argv[1]
    if len(sys.argv) > 2:
        top = sys.argv[2]
    if len(sys.argv) > 3:
        temp = float(sys.argv[3])
    # optional fourth arg: single residue number or comma-separated list (e.g. "203" or "203,205")
    target_resnums = None
    if len(sys.argv) > 4:
        try:
            target_resnums = [int(x) for x in sys.argv[4].split(',') if x.strip()]
        except Exception:
            target_resnums = None

    main(top=top, min_base=min_base, temp=temp, target_resnums=target_resnums)
