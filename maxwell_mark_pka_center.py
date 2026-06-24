#!/usr/bin/env python3
import re
import math
import random
import sys

def maxwell_velocity_sigma(temperature, mass_amu):
    """
    计算麦克斯韦-玻尔兹曼速度分布的标准差
    sigma = sqrt(kB * T / m)
    单位: Å/fs
    """
    k_B = 8.617333262145e-5  # eV/K
    # 1 amu = 1.66053906660e-27 kg
    # 1 eV = 1.602176634e-19 J
    # 转换因子: 从 (eV/amu) 到 (m/s)^2
    # v^2 = 2 * E / m, 但这里只需要 sigma
    # sigma = sqrt(k_B * T / m) in m/s
    mass_kg = mass_amu * 1.66053906660e-27
    sigma_mps = math.sqrt(k_B * temperature / (mass_kg / 1.602176634e-19))
    # 转换到 Å/fs: 1 m/s = 1e-5 Å/fs
    sigma_afs = sigma_mps * 1e-5
    return sigma_afs

def generate_thermal_velocity(sigma):
    """生成高斯分布的随机速度 (Å/fs)"""
    return random.gauss(0.0, sigma)

def process_xyz(filename, temperature=300, pka_energy_kev=2, pka_direction=(1, 0, 0)):
    """
    处理 xyz 文件，标记 PKA 原子并初始化麦克斯韦速度
    
    参数:
        filename: 输入 xyz 文件名
        temperature: 温度 (K)，默认 300
        pka_energy_kev: PKA 能量 (keV)，默认 2
        pka_direction: PKA 方向向量，默认 (1,0,0) 即 x 方向
    """
    # 设置随机种子以便结果可重复
    random.seed(42)
    
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 找到第一帧
    n_atoms_line = 0
    for i, line in enumerate(lines):
        if line.strip().isdigit():
            n_atoms_line = i
            break
    n_atoms = int(lines[n_atoms_line].strip())
    header = lines[n_atoms_line + 1]
    atom_lines = lines[n_atoms_line+2 : n_atoms_line+2+n_atoms]

    # 解析 Lattice 得到盒子中心
    match = re.search(r'Lattice="([^"]*)"', header)
    if not match:
        raise ValueError("Cannot find Lattice in header")
    lattice_vals = list(map(float, match.group(1).split()))
    Lx = lattice_vals[0]
    Ly = lattice_vals[4] if len(lattice_vals) > 4 else Lx
    Lz = lattice_vals[8] if len(lattice_vals) > 8 else Lx
    center = (Lx/2.0, Ly/2.0, Lz/2.0)

    print(f"Box dimensions: {Lx:.6f} x {Ly:.6f} x {Lz:.6f} Å")
    print(f"Box center: ({center[0]:.6f}, {center[1]:.6f}, {center[2]:.6f}) Å")
    print(f"Temperature: {temperature} K")
    print(f"PKA energy: {pka_energy_kev} keV")

    # 铁原子质量 (amu)
    mass_Fe_amu = 55.845
    
    # 计算麦克斯韦速度分布的标准差
    sigma_afs = maxwell_velocity_sigma(temperature, mass_Fe_amu)
    print(f"Maxwell velocity sigma: {sigma_afs:.6f} Å/fs")
    
    # 计算 PKA 定向速度
    # E = 1/2 * m * v^2
    # v = sqrt(2E/m)
    pka_energy_J = pka_energy_kev * 1000 * 1.602176634e-19
    mass_kg = mass_Fe_amu * 1.66053906660e-27
    pka_speed_mps = math.sqrt(2 * pka_energy_J / mass_kg)
    pka_speed_afs = pka_speed_mps * 1e-5  # 转换到 Å/fs
    # 归一化方向向量
    dir_norm = math.sqrt(sum(d**2 for d in pka_direction))
    pka_vx = pka_speed_afs * pka_direction[0] / dir_norm
    pka_vy = pka_speed_afs * pka_direction[1] / dir_norm
    pka_vz = pka_speed_afs * pka_direction[2] / dir_norm
    print(f"PKA speed: {pka_speed_afs:.6f} Å/fs")
    print(f"PKA direction: {pka_direction}")
    print(f"PKA velocity: ({pka_vx:.6f}, {pka_vy:.6f}, {pka_vz:.6f}) Å/fs")

    # 找到离中心最近的原子作为 PKA
    best_idx = -1
    best_dist2 = float('inf')
    best_coord = None
    for i, line in enumerate(atom_lines):
        parts = line.split()
        if len(parts) < 4 or parts[0] != 'Fe':
            continue
        try:
            x, y, z = map(float, parts[1:4])
        except:
            continue
        dx = x - center[0]
        dy = y - center[1]
        dz = z - center[2]
        dist2 = dx*dx + dy*dy + dz*dz
        if dist2 < best_dist2:
            best_dist2 = dist2
            best_idx = i
            best_coord = (x, y, z)

    if best_idx == -1:
        raise RuntimeError("No valid Fe atoms found")

    print(f"\nClosest atom to center: atom index {best_idx} (0-based in atom list)")
    print(f"  Coordinates: ({best_coord[0]:.6f}, {best_coord[1]:.6f}, {best_coord[2]:.6f})")
    print(f"  Distance from center: {math.sqrt(best_dist2):.6f} Å")

    # 修改 header 添加 vel 和 group 属性
    # 确保有 vel:R:3 和 group:I:2 (group 可支持 0,1,2)
    prop_pat = re.compile(r'([Pp]roperties=)([^\s]+)')
    def modify_properties(m):
        prefix = m.group(1)
        props = m.group(2)
        if 'vel:R:3' not in props:
            props = props + ':vel:R:3'
        if 'group:I:2' not in props:
            # 如果已有 group:I:1，替换为 group:I:2
            if 'group:I:1' in props:
                props = props.replace('group:I:1', 'group:I:2')
            else:
                props = props + ':group:I:2'
        return prefix + props
    new_header = prop_pat.sub(modify_properties, header)
    if 'properties=' not in new_header.lower():
        new_header = new_header.rstrip() + ' properties=species:S:1:pos:R:3:vel:R:3:group:I:2\n'

    # 生成新的原子行（包含速度和 group）
    new_atom_lines = []
    for i, line in enumerate(atom_lines):
        stripped = line.rstrip()
        parts = stripped.split()
        if len(parts) >= 4 and parts[0] == 'Fe':
            # 生成热速度
            vx_th = generate_thermal_velocity(sigma_afs)
            vy_th = generate_thermal_velocity(sigma_afs)
            vz_th = generate_thermal_velocity(sigma_afs)
            
            if i == best_idx:
                # PKA 原子：热速度 + 定向速度
                vx = vx_th + pka_vx
                vy = vy_th + pka_vy
                vz = vz_th + pka_vz
                group = 1
                new_atom_lines.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]} {vx:.8f} {vy:.8f} {vz:.8f} {group}\n")
            else:
                # 普通原子：只有热速度
                vx = vx_th
                vy = vy_th
                vz = vz_th
                group = 0
                new_atom_lines.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]} {vx:.8f} {vy:.8f} {vz:.8f} {group}\n")
        else:
            new_atom_lines.append(stripped + '\n')

    # 写入新文件
    outfile = 'model_pka.xyz'
    with open(outfile, 'w') as f:
        f.write(f"{n_atoms}\n")
        f.write(new_header)
        f.writelines(new_atom_lines)

    print(f"\nWritten {outfile}")
    print(f"All atoms: Maxwell-Boltzmann velocity distribution at {temperature} K")
    print(f"PKA atom (group 1): thermal + directional velocity ({pka_energy_kev} keV)")
    print(f"Other atoms (group 0): thermal velocity only")
    print("\nNow you can use this file for PKA simulation.")

if __name__ == "__main__":
    # ========== 可修改参数 ==========
    input_file = 'equilibrated.xyz'      # 输入文件名
    temperature = 30                  # 温度 (K)
    pka_energy_kev = 0.018            # PKA 能量 (keV)
    pka_direction = (1, 0, 0)      # PKA 方向，例如 (1,0,0) 为 x 方向
    # ==============================
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        temperature = float(sys.argv[2])
    if len(sys.argv) > 3:
        pka_energy_kev = float(sys.argv[3])
    
    process_xyz(input_file, temperature, pka_energy_kev, pka_direction)