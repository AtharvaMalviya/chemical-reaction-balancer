# Computational Chemistry Python Project

This project is a small command-line tool written in Python that performs three basic computational-chemistry tasks:

1. **Calculate the molecular mass of a chemical formula**
2. **Find percent composition of each element in a compound**
3. **Balance chemical equations using matrix methods**

It uses a formula parser that supports nested parentheses (for example: `Al2(SO4)3`) and a balancing method based on solving a homogeneous system of linear equations using exact rational arithmetic.

---

## Features

### **1. Molecular Mass Calculator**

* Accepts any chemical formula made of valid element symbols
* Handles parentheses
* Uses the full periodic table (118 elements)
* Outputs the molar mass in g/mol

### **2. Percent Composition**

* Shows the mass percentage of every element in the compound
* Uses atomic masses from a built-in periodic table
* Percentages are calculated relative to the total molar mass

### **3. Chemical Equation Balancer**

* Guided input: user enters reactants and products separately
* Automatically identifies all elements involved
* Builds and solves a linear algebra system using `fractions.Fraction`
* Produces integer coefficients for the balanced equation

---

## Example Usage

### **Molar Mass**

```
Enter chemical formula: H2O
Molar mass: 18.015
```

### **Percent Composition**

```
Enter chemical formula: C6H12O6
C: 40.00%
H: 6.71%
O: 53.29%
```

### **Balancing Equations**

Input:

```
Number of reactants: 2
Reactant 1: C2H6
Reactant 2: O2
Number of products: 2
Product 1: CO2
Product 2: H2O
```

Output:

```
Balanced: 2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O
```

---

## Project Structure

```
main.py
README.md
```

Everything is contained inside a single Python file for simplicity.

---

## How to Run

1. Install Python 3.8+
2. Download or copy `main.py`
3. Run:

```
python main.py
```

4. Choose an option from the menu and follow the instructions.

---

## Requirements

No external packages are needed.

The program only uses:

* `fractions` (built-in)
* `re` (built-in)

---

## Purpose

This project was created as part of an undergraduate assignment in **computational chemistry**.
It focuses on applying programming concepts like:

* Parsing
* Dictionaries
* Fractions
* Lists
* Matrix operations
* Linear algebra

While keeping the implementation clear and beginner-friendly.

---

## Notes

* Atomic masses are taken from standard values (IUPAC averages).
* The equation-balancing algorithm avoids floating point errors by using exact fractions.
* The program has been tested on multiple formulas and reactions.

---
