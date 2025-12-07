from fractions import Fraction
import re

atomic_weights = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867,
    "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98,
    "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91,
    "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209, "At": 210,
    "Rn": 222, "Fr": 223, "Ra": 226, "Ac": 227, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237,
    "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258,
    "No": 259, "Lr": 266, "Rf": 267, "Db": 268, "Sg": 269, "Bh": 270, "Hs": 269, "Mt": 278,
    "Ds": 281, "Rg": 282, "Cn": 285, "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293, "Ts": 294, "Og": 294
}


# ---- Formula Parser ----

def parse_formula(formula):
    tokens = re.findall(r'([A-Z][a-z]?|\(|\)|\d+)', formula)
    stack = [{}]

    i = 0
    while i < len(tokens):
        t = tokens[i]

        if t == "(":
            stack.append({})
        elif t == ")":
            group = stack.pop()
            i += 1
            mult = int(tokens[i]) if i < len(tokens) and tokens[i].isdigit() else 1
            for el, count in group.items():
                stack[-1][el] = stack[-1].get(el, 0) + count * mult
        elif t.isdigit():
            for el in stack[-1]:
                if stack[-1][el] == "_last":
                    stack[-1][el] = int(t)
                    break
        else:
            stack[-1][t] = "_last"

        i += 1

    for el in list(stack[0].keys()):
        if stack[0][el] == "_last":
            stack[0][el] = 1

    return stack[0]


# ---- Molecular mass ----

def molar_mass(formula):
    parsed = parse_formula(formula)
    total = 0
    for el, count in parsed.items():
        total += atomic_weights[el] * count
    return total


def percent_composition(formula):
    parsed = parse_formula(formula)
    total = molar_mass(formula)
    result = {}
    for el, count in parsed.items():
        result[el] = (atomic_weights[el] * count / total) * 100
    return result


# ---- Balancing Chemical Equations ----

def build_matrix(reactants, products):
    all_molecules = reactants + products
    elements = set()

    parsed_list = []
    for mol in all_molecules:
        p = parse_formula(mol)
        parsed_list.append(p)
        elements.update(p.keys())

    elements = list(elements)

    rows = []
    for el in elements:
        row = []
        for i, p in enumerate(parsed_list):
            row.append(Fraction(p.get(el, 0)) * (1 if i < len(reactants) else -1))
        rows.append(row)

    return rows


def rref(matrix):
    mat = [row[:] for row in matrix]
    r, c = 0, 0
    rows, cols = len(mat), len(mat[0])

    while r < rows and c < cols:
        pivot = None
        for i in range(r, rows):
            if mat[i][c] != 0:
                pivot = i
                break

        if pivot is None:
            c += 1
            continue

        mat[r], mat[pivot] = mat[pivot], mat[r]

        pivot_val = mat[r][c]
        mat[r] = [val / pivot_val for val in mat[r]]

        for i in range(rows):
            if i != r and mat[i][c] != 0:
                factor = mat[i][c]
                mat[i] = [mat[i][j] - factor * mat[r][j] for j in range(cols)]

        r += 1
        c += 1

    return mat


def null_space_vector(matrix):
    rows = rref(matrix)
    cols = len(rows[0])
    solution = [Fraction(0)] * cols
    solution[-1] = 1

    for i in reversed(range(rows.__len__())):
        pivot_col = None
        for j in range(cols):
            if rows[i][j] == 1:
                pivot_col = j
                break
        if pivot_col is None:
            continue

        s = 0
        for j in range(pivot_col + 1, cols):
            s += rows[i][j] * solution[j]
        solution[pivot_col] = -s

    lcm = 1
    for v in solution:
        if v.denominator != 1:
            lcm = (lcm * v.denominator) // gcd(lcm, v.denominator)

    return [int(v * lcm) for v in solution]


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def balance_equation(reactants, products):
    matrix = build_matrix(reactants, products)
    coeffs = null_space_vector(matrix)

    r = coeffs[:len(reactants)]
    p = coeffs[len(reactants):]

    return r, p


# ---- Main Menu ----

def main():
    while True:
        print("\n--- Computational Chemistry Tool ---")
        print("1. Molecular Mass")
        print("2. Percent Composition")
        print("3. Balance Chemical Equation")
        print("4. Exit")

        choice = input("Enter choice: ")

        if choice == "1":
            f = input("Enter chemical formula: ")
            print("Molar mass:", molar_mass(f))

        elif choice == "2":
            f = input("Enter chemical formula: ")
            pc = percent_composition(f)
            for el, pct in pc.items():
                print(f"{el}: {pct:.2f}%")

        elif choice == "3":
            n = int(input("Number of reactants: "))
            react = [input(f"Reactant {i+1}: ") for i in range(n)]

            m = int(input("Number of products: "))
            prod = [input(f"Product {i+1}: ") for i in range(m)]

            r, p = balance_equation(react, prod)
            out = ""
            for i in range(n):
                out += f"{r[i]} {react[i]}"
                if i < n - 1:
                    out += " + "
            out += " -> "
            for i in range(m):
                out += f"{p[i]} {prod[i]}"
                if i < m - 1:
                    out += " + "

            print("Balanced:", out)

        else:
            break


if __name__ == "__main__":
    main()
