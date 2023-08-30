# Solving-Partial-Differential-Equations-using-Python

This repository is dedicated to solving Partial Differential Equations (PDEs) using Python. We specifically focus on the Heat Equation, the Wave Equation, and the uncoupled Thermo-elasticity problem. Solutions are obtained using both Forward Time Central Space (FTCS) and Backward Time Central Space (BTCS) numerical schemes.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [License](#license)

## Introduction

Partial Differential Equations (PDEs) are fundamental in modeling numerous physical phenomena. This repository provides a set of Python scripts that solve some common PDEs:

- **Heat Equation**: Describes the distribution of heat (or variation in temperature) in a given region over time.
  
- **Wave Equation**: Represents the behavior of waves, be it sound, light, or water waves.
  
- **Uncoupled Thermo-elasticity**: Represents the relation between thermal and elastic deformations, without them being directly coupled.

## Installation

**Prerequisites**: Ensure you have Python 3.x installed.

1. Clone the repository:

   ```bash
   git clone https://github.com/shadzzz90/Solving-Partial-Differential-Equations-using-Python.git
   ```

2. Navigate into the directory:

   ```bash
   cd Solving-Partial-Differential-Equations-using-Python
   ```

3. (Optional) It's recommended to create a virtual environment:

   ```bash
   python -m venv env
   source env/bin/activate  # On Windows, use `env\Scripts\activate`
   ```

4. Install the required packages:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

Each PDE and its respective solution method is encapsulated in a separate Python script. To run any of the scripts:

```bash
python [script_name].py
```

Replace `[script_name]` with the desired script (e.g., `heat_equation_FTCS.py`).

## Features

- **Numerical Schemes**:
  - FTCS: Forward Time Central Space
  - BTCS: Backward Time Central Space

- **Visualizations**: The solutions are visualized using matplotlib, showing how the quantities evolve over time.


## License

This project is under the MIT License. See the [LICENSE.md](LICENSE.md) file for details.

---

Happy solving!# Solving-Partial-Differential-Equations-using-Python
